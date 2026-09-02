# get_survival(): the analytical core of the BRCA survival tool.
#
# Takes ANY per-sample score (a gene's expression, a VIPER NES, a signature score)
# as a named numeric vector (names = sample_id, spanning any cohorts) and returns:
#   - per-cohort continuous Cox model (score as a per-SD covariate), stratified by
#     that cohort's registry stratifier when it has one (pam50 for breast, none for
#     ovarian). A stratified Cox assumes no PH across strata, which is what we want
#     where the stratum is a strong confounder (PAM50 subtype in BRCA);
#   - per-cohort median-split KM + log-rank (for visualization / an intuitive HR);
#   - an inverse-variance random-effects meta-pool of the continuous log(HR) across
#     cohorts (metafor::rma, REML), with I^2 heterogeneity.
#
# Design: scores are network-agnostic, so the same engine serves expression, VIPER,
# and signature queries. Because expr/VIPER matrices are z-scored within cohort, the
# continuous coefficient is a per-SD log(HR) and is directly poolable across cohorts.

suppressPackageStartupMessages({ library(survival); library(metafor) })
# Requires get_clinical()/SURV_COHORTS from data_access.R (source it before this file).

# Proportional-hazards check for the score term (Grambsch-Therneau: correlation of
# scaled Schoenfeld residuals with time). A SMALL p means the coefficient moves over
# follow-up, i.e. PH is VIOLATED and the reported HR is a time-average of a changing
# effect rather than one constant ratio. Returns NA rather than erroring on degenerate
# fits (too few events, zero-variance score) -- a missing diagnostic must not take down
# a query that otherwise succeeded.
PH_ALPHA <- 0.05

# Below this many cohorts, tau^2 (the between-study variance the random-effects CI depends
# on) is not reliably identified, so the pooled CI reports a heterogeneity estimate carrying
# almost no information. At k=2 the same inputs give REML [0.663, 0.915] p=0.0024 and
# Knapp-Hartung [0.273, 2.221] p=0.20 -- a stable point estimate wrapped in an arbitrary
# interval. The whole breast panel sits here (k=2 DFS, k=3 OS).
#
# The engine does NOT switch estimators below the threshold. Fixed-effect would be
# affirmatively wrong for breast OS, where 64% of scanned TFs show I2 > 50% -- it assumes
# the homogeneity the data rejects -- and Knapp-Hartung at k=2 is a t with 1 df. All three
# rules are differently unreliable; picking a different one would silently re-rank ~1200
# breast hits on an assumption the data does not support. So the pooled row DISCLOSES that
# tau^2 is unidentified and carries the fixed-effect result beside it, letting a reader see
# whether a conclusion depends on the pooling choice. Annotate, don't substitute -- the same
# rule as the "!" PH marker and the truncation labels.
#
# 5 is the conventional floor for estimating tau^2 with any reliability. Like tau=48 it is
# declared from study design, not tuned against the resulting p-values.
POOL_MIN_K <- 5

# REML tau^2 by Fisher scoring, with metafor's own documented remedy on non-convergence.
#
# rma(method = "REML") ERRORS -- not warns -- with "Fisher scoring algorithm did not
# converge" on a minority of inputs. It propagated out of get_survival(), and
# batch_regulators()'s per-TF tryCatch then discarded the row, so the regulator simply
# was not in the scan CSV and nothing in the file said it had been tested. 110 of 1446
# ovarian TFs were missing that way, and NOT at random: the failure needs a set of yi
# that is tight, consistently signed, and precise, which is exactly the profile of a
# strong effect (HIC1, NR2F1, ZFHX4, ZNF521 all reach |z| > 6 once they fit). The scan
# was dropping its best candidates and calling it silence.
#
# What this does and does not change. Damping the Fisher step (control$stepadj) is a
# property of the OPTIMIZER, not of the model: same REML likelihood, same tau^2, same
# estimand. It is the first remedy help(rma) lists. The default call is attempted FIRST
# and returned untouched whenever it succeeds, so every estimate this project has ever
# published is bit-identical -- the retry can only produce rows that used to be absent.
# (Damping everything unconditionally would have been simpler and was rejected: it
# shifts converging fits by ~1e-5 in tau^2, which is numerically nothing but would still
# have re-ranked already-reported results for no gain.)
#
# Two damping levels, then give up and let the error propagate. Silently returning some
# other estimator -- DL, or the fixed-effect fit -- would be a substitution, and this
# engine annotates rather than substitutes (same rule as the "!" PH marker, the tau^2
# disclosure, and the truncation labels). A row that needed damping SAYS so, via
# converged_by, all the way out to the scan CSV.
.pool_rma <- function(yi, vi) {
  fit <- tryCatch(rma(yi = yi, vi = vi, method = "REML"), error = function(e) NULL)
  if (!is.null(fit)) return(list(fit = fit, converged_by = "default"))
  for (ctl in list(list(stepadj = 0.5,  maxiter = 1000),
                   list(stepadj = 0.25, maxiter = 5000))) {
    fit <- tryCatch(rma(yi = yi, vi = vi, method = "REML", control = ctl),
                    error = function(e) NULL)
    if (!is.null(fit)) return(list(fit = fit, converged_by = "damped"))
  }
  stop("REML tau^2 did not converge for this set of cohorts, even with a damped ",
       "Fisher step (stepadj 0.5 then 0.25). k = ", length(yi))
}

.ph_p <- function(fit) {
  z <- tryCatch(cox.zph(fit), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(z) || !("score" %in% rownames(z$table))) return(NA_real_)
  p <- z$table["score", "p"]
  if (!is.finite(p)) return(NA_real_)
  unname(p)
}

get_survival <- function(score,
                         endpoint = c("OS", "DFS", "DSS"),
                         cohorts = NULL,
                         cancer_type = NULL,
                         adjust_strata = TRUE,
                         standardize = TRUE,
                         min_events = 10,
                         max_followup = NULL,
                         clinical = NULL) {
  endpoint <- match.arg(endpoint)
  # FAIL LOUD on a score that cannot produce a result. Before this the engine joined an
  # unusable score, matched nothing, and returned an EMPTY survresult with no error -- which
  # only surfaced later in forest_plot() as "no estimable cohort to plot", a message that
  # names the symptom, not the cause (a typo'd feature, an unnamed vector, the wrong
  # cohorts). A score MUST be a named numeric vector (names = sample_id); an unnamed or
  # non-numeric one is always a caller bug, never a legitimate "no data" outcome. (The one
  # place an empty score is legitimate -- a TF measurable in no cohort -- is guarded in
  # batch_regulators() before it reaches here, so it is counted, not raised.)
  if (!is.numeric(score))
    stop("`score` must be a named numeric vector (names = sample_id); got ",
         if (is.null(score)) "NULL" else paste(class(score), collapse = "/"))
  if (length(score) > 0 && is.null(names(score)))
    stop("`score` must be NAMED (names = sample_id) so it can be joined to clinical rows; ",
         "got an unnamed length-", length(score), " vector")
  if (!is.null(max_followup)) {
    if (!is.numeric(max_followup) || length(max_followup) != 1 ||
        !is.finite(max_followup) || max_followup <= 0)
      stop("max_followup must be a single positive number of months (or NULL)")
  }
  # Resolve cohorts: explicit cohorts win; else resolve from cancer_type (endpoint-
  # filtered via the registry); else fall back to the legacy breast default so existing
  # breast drivers (make_plots/validate_engine/signature_score) keep working unchanged.
  if (is.null(cohorts)) {
    if (!is.null(cancer_type)) {
      # Name the real cause: an unknown cancer_type otherwise resolves to character(0) and
      # dies downstream in validate_cohorts() as "no cohorts requested (empty cohort
      # vector)", blaming the caller for something they did not pass.
      if (!(cancer_type %in% COHORTS$cancer_type))
        stop("unknown cancer_type '", cancer_type, "'. Known: ",
             paste(unique(COHORTS$cancer_type), collapse = ", "))
      cohorts <- cohorts_for(cancer_type, endpoint)
      if (!length(cohorts))
        stop("cancer_type '", cancer_type, "' has no cohort carrying endpoint ", endpoint)
    } else cohorts <- SURV_COHORTS
  }
  # `clinical` lets a SEPARATE measurement layer reuse this engine without being a
  # registered cohort. RPPA is the case it exists for: those samples are not in
  # clinical.db and their cohorts are deliberately absent from config/cohorts.tsv,
  # because cohorts_for(cancer_type, endpoint) is called in a dozen places with no role
  # filter and every one of them means "the expression cohorts for this tissue" -- a
  # registered RPPA cohort would surface as an extra forest row drawn from patients
  # already on the plot. See R/rppa.R.
  #
  # Everything below this line is untouched by the choice: with explicit `cohorts` and
  # adjust_strata = FALSE, get_survival() reads the registry NOWHERE else, so an
  # injected frame flows through the same guards, standardization, Cox and KM as a
  # registry-backed one. The guards here are what keep that true.
  if (is.null(clinical)) {
    clin <- get_clinical(cohorts, endpoint)
  } else {
    if (!is.null(cancer_type))
      stop("`clinical` and `cancer_type` are mutually exclusive: an injected clinical ",
           "frame carries its own cohorts, so there is nothing for the registry to resolve")
    if (!is.data.frame(clinical))
      stop("`clinical` must be a data.frame; got ", paste(class(clinical), collapse = "/"))
    need <- c("sample_id", "dataset_id", "time", "event")
    miss <- setdiff(need, names(clinical))
    if (length(miss))
      stop("`clinical` is missing column(s): ", paste(miss, collapse = ", "),
           ". Required: ", paste(need, collapse = ", "))
    if (nrow(clinical) == 0)
      stop("`clinical` has zero rows -- nothing to model")
    # Without this, an injected frame whose dataset_id values do not appear in `cohorts`
    # is skipped one cohort at a time by the loop below and returns an EMPTY result with
    # no error -- the same silent-empty failure the score guards above exist to prevent.
    stray <- setdiff(unique(clinical$dataset_id), cohorts)
    if (length(stray))
      stop("`clinical` carries dataset_id(s) not in `cohorts`: ",
           paste(stray, collapse = ", "), ". Requested cohorts: ",
           paste(cohorts, collapse = ", "))
    if (adjust_strata)
      stop("`clinical` requires adjust_strata = FALSE: the per-cohort stratifier comes ",
           "from the cohort registry, which an injected frame is by definition not in")
    clin <- clinical
  }
  clin$score <- score[clin$sample_id]
  # Overlap captured BEFORE the validity filter so the fail-loud message below can tell
  # "the score's names matched no clinical row" (wrong cohorts / unnamed / typo) apart from
  # "matched, but every matched row was unusable" (all-non-finite score, non-positive time).
  n_clin <- nrow(clin); n_matched <- sum(!is.na(clin$score))
  clin <- clin[is.finite(clin$score) & is.finite(clin$time) &
               !is.na(clin$event) & clin$time > 0, ]
  # Zero usable rows is a caller error, not an empty result to return silently. It is the
  # distinct opposite of the LEGITIMATE thin case (rows matched but every cohort fell below
  # min_events): that keeps nrow > 0 and flows through as per-cohort skips with a NULL pool.
  # Only a total non-overlap lands here, so raising cannot swallow a real result.
  if (nrow(clin) == 0)
    stop(if (n_matched == 0)
           sprintf("`score` matched none of the %d clinical sample_ids for cohorts {%s} at endpoint %s (score carries %d name(s)). Check the names are sample_ids for these cohorts.",
                   n_clin, paste(cohorts, collapse = ", "), endpoint, length(score))
         else
           sprintf("`score` matched %d of %d clinical samples, but none had a finite score, a positive time, and a non-missing %s event -- no rows to model.",
                   n_matched, n_clin, endpoint))

  # Administrative censoring at a common horizon. Applied to ALL cohorts at once and
  # BEFORE anything downstream (event counting, min_events, standardization, Cox, KM,
  # cox.zph), so every number in the result -- including the PH verdict -- describes
  # the truncated window rather than being inherited from the full one.
  #
  # Censor, don't exclude: a patient still observed at tau keeps their follow-up up to
  # tau and becomes event-free there. Dropping them instead would delete exactly the
  # longest survivors and bias the sample toward early deaths.
  if (!is.null(max_followup)) {
    late <- clin$time > max_followup
    clin$event[late] <- 0
    clin$time[late]  <- max_followup
  }

  per <- list()
  for (coh in cohorts) {
    d <- clin[clin$dataset_id == coh, ]
    if (nrow(d) == 0) next
    # Standardize the score WITHIN cohort so the continuous Cox coefficient is a
    # per-SD log(HR) regardless of input type (expression, VIPER NES, signature)
    # -> HRs are interpretable and poolable on one scale. (Idempotent for the
    # already-z-scored expression matrices; monotonic, so median split unaffected.)
    # Reference population is the full cohort d, deliberately: the per-SD unit is a
    # property of the cohort, not of whichever subset survives stratifier filtering.
    if (standardize) {
      s <- d$score; sdev <- sd(s)
      if (is.finite(sdev) && sdev > 0) d$score <- (s - mean(s)) / sdev
    }
    # Per-cohort stratifier from the registry (pam50 for breast, none/stage for
    # others). strata() lets each stratum keep its own baseline hazard.
    # ENDPOINT-SCOPED since step 102: `stratify_by` may differ per endpoint, because a
    # cohort sits in more than one pool and its partners differ between them. Passing
    # `endpoint` is not optional here -- stratifier_for() stops rather than guess, and this
    # is the call whose answer IS the estimand.
    strat_var <- if (adjust_strata) stratifier_for(coh, endpoint) else NA_character_
    # FAIL LOUD when the registry declares a stratifier this frame cannot supply.
    #
    # The gate below is a conjunction, and `strat_var %in% names(d)` used to fail SILENTLY
    # inside it: a declared stratifier that is not a column made use_strat FALSE, the cohort
    # fitted UNSTRATIFIED, and the result was indistinguishable from a legitimately
    # unstratified one -- same fields, no warning, a normal-looking HR straight into the
    # meta-pool. A missing COLUMN is a schema or plumbing bug (the loader does not SELECT it;
    # the registry names something that does not exist), never a property of the patients.
    #
    # The two conditions AFTER it are the opposite kind of thing: a cohort with no values, or
    # with only one level, is a data fact. Those still fall back to unstratified on purpose --
    # refusing them would make a tissue unqueryable because one arm lacks the annotation --
    # and the pooled result now DISCLOSES the resulting mixture instead (see `strata_k`).
    #
    # The `clinical =` injection path cannot reach here: it forces adjust_strata = FALSE, so
    # strat_var is NA and this check is skipped, exactly as an unregistered frame requires.
    if (!is.na(strat_var) && !(strat_var %in% names(d)))
      stop(sprintf(paste0(
        "cohort '%s' declares stratify_by = '%s', but the clinical frame has no '%s' column, ",
        "so the Cox would have fitted UNSTRATIFIED with no warning. Columns present: %s. ",
        "Either the clinical loader does not SELECT '%s' yet, or the registry names a column ",
        "that does not exist."),
        coh, strat_var, strat_var, paste(names(d), collapse = ", "), strat_var))
    use_strat <- !is.na(strat_var) &&
      sum(!is.na(d[[strat_var]])) > 0 && length(unique(na.omit(d[[strat_var]]))) >= 2
    # d2 is the sample the Cox actually fits: under stratification, rows with a missing
    # stratifier are dropped. Derived BEFORE the min_events guard so the guard counts the
    # events that will be modelled, not a larger pre-filter set -- otherwise a cohort at
    # the boundary passes on events it never uses (TCGA_BRCA/DFS: 145 in d, 144 in d2).
    d2 <- if (use_strat) d[!is.na(d[[strat_var]]), ] else d

    n_ev <- sum(d2$event == 1)
    if (n_ev < min_events) {
      per[[coh]] <- list(n = nrow(d2), events = n_ev, skipped = TRUE,
                         reason = sprintf("only %d events (<%d)", n_ev, min_events))
      next
    }

    # A score with no within-cohort spread has no contrast for the Cox to fit, and collapses
    # the median split to a single group -- survdiff() would then abort the WHOLE query
    # ("There is only 1 group"), so one degenerate cohort takes down every other cohort's
    # pooled result. Treat it like too-few-events: a named per-cohort skip, not a crash. Only
    # an EXACTLY-constant score lands here; a merely tiny variance still yields a valid per-SD
    # coefficient (the reference SD is a real number) and is left to fit -- thresholding the
    # variance would be a shoppable cutoff, which this engine refuses (declarable, not tuned).
    if (length(unique(d2$score)) < 2L) {
      per[[coh]] <- list(n = nrow(d2), events = n_ev, skipped = TRUE,
                         reason = "score is constant in this cohort (no contrast to fit)")
      next
    }

    # --- continuous Cox (per-SD score) ---
    fmla <- if (use_strat)
      as.formula(sprintf("Surv(time, event) ~ score + strata(%s)", strat_var))
    else Surv(time, event) ~ score
    # coxph does NOT error on a non-convergent fit. On quasi-complete separation it returns a
    # finite-but-diverging coefficient (logHR ~403, HR ~1e175) with only a WARNING; on an
    # all-censored arm it returns coef = NA. Either would otherwise be poured into the meta-pool
    # as a real per-SD log(HR) and, in a scan, top the results on a numerical artifact. Capture
    # coxph's OWN convergence verdict (its warning) and reject a non-finite estimate, then skip
    # with a named reason -- trusting the fitter's verdict rather than a magnitude threshold we
    # invent, the same way .pool_rma() trusts REML's convergence. Non-convergence warnings are
    # consumed here; any other coxph warning still propagates.
    conv_warn <- NULL
    fit <- withCallingHandlers(
      coxph(fmla, data = d2),
      warning = function(w) {
        m <- conditionMessage(w)
        if (grepl("did not converge|infinite|Loglik converged before", m)) {
          conv_warn <<- m; invokeRestart("muffleWarning")
        }
      })
    co  <- summary(fit)$coefficients["score", ]
    logHR <- unname(co["coef"]); se <- unname(co["se(coef)"]); p <- unname(co["Pr(>|z|)"])
    if (!is.null(conv_warn) || !is.finite(logHR) || !is.finite(se) || se <= 0) {
      per[[coh]] <- list(n = nrow(d2), events = n_ev, skipped = TRUE,
                         reason = if (!is.null(conv_warn))
                           sprintf("Cox did not converge (%s)", conv_warn)
                         else "Cox produced a non-finite coefficient (all-censored arm or separation)")
      next
    }
    ph_p <- .ph_p(fit)

    # --- median-split KM + log-rank (visualization / intuitive HR) ---
    # Built on d2 -- the SAME rows the Cox above fitted -- not on d. The KM is a view of
    # this result, so it has to describe this result's patients: km_plot() draws these
    # curves and annotates them with the Cox HR and the "<strat>-stratified" label, and a
    # figure whose log-rank, HR, and label come from three different samples cannot be
    # reconciled against itself, against the forest slab, or against survtable()'s CSV.
    # (Under an unstratified fit d2 IS d, so this only ever differs where a stratifier is
    # missing on some rows: TCGA_BRCA/DFS, 12 rows and 1 event. METABRIC happens to have
    # no gap, which is why the divergence stayed invisible.)
    grp <- factor(ifelse(d2$score >= median(d2$score), "high", "low"),
                  levels = c("low", "high"))
    sd_ <- survdiff(Surv(time, event) ~ grp, data = d2)
    lr_p <- pchisq(sd_$chisq, df = length(sd_$n) - 1, lower.tail = FALSE)

    per[[coh]] <- list(
      n = nrow(d2), events = sum(d2$event == 1),
      logHR = logHR, se = se, HR = exp(logHR), p = p,
      ph_p = ph_p, ph_violated = !is.na(ph_p) && ph_p < PH_ALPHA,
      strata_adjusted = use_strat, strata_var = if (use_strat) strat_var else NA_character_, logrank_p = lr_p,
      km = survfit(Surv(time, event) ~ grp, data = d2),
      data = d2, skipped = FALSE)
  }

  # --- meta-pool the continuous log(HR) ---
  est <- Filter(function(x) isFALSE(x$skipped), per)
  # How many of the cohorts that ACTUALLY ENTER the pool were fitted with a strata() term,
  # and under which variable. Counted over `est`, not over `cohorts`: a cohort skipped for
  # too few events is not in the average and must not be in its denominator either.
  #
  # This exists because a pooled HR can silently average per-cohort estimates made under
  # DIFFERENT models. The per-cohort `strata_adjusted` flag has always been recorded, but the
  # pooled row is where the number a reader quotes lives, and nothing there said whether the
  # cohorts behind it were fitted the same way. With one stratifier (breast/pam50) that gap
  # was invisible because the three cohorts all carry PAM50; with a stratifier only some
  # cohorts can supply it is the default case, not the edge case.
  .k_strat <- sum(vapply(est, function(e) isTRUE(e$strata_adjusted), logical(1)))
  .s_vars  <- unique(na.omit(vapply(est, function(e) e$strata_var, character(1))))
  pooled <- NULL
  if (length(est) >= 2) {
    yi <- vapply(est, `[[`, numeric(1), "logHR")
    vi <- vapply(est, `[[`, numeric(1), "se")^2
    pr <- .pool_rma(yi, vi)
    m  <- pr$fit
    # Fixed-effect fit on the SAME inputs, carried as a sensitivity rather than as the
    # answer. Where tau^2 = 0 the two coincide exactly, which is itself the useful signal:
    # it says the pooling choice does not matter for this row.
    fe <- rma(yi = yi, vi = vi, method = "FE")
    pooled <- list(HR = exp(as.numeric(m$b)), ci_lb = exp(m$ci.lb), ci_ub = exp(m$ci.ub),
                   p = m$pval, I2 = m$I2, tau2 = m$tau2, k = m$k,
                   tau2_identified = m$k >= POOL_MIN_K,
                   converged_by = pr$converged_by,
                   HR_fe = exp(as.numeric(fe$b)), ci_fe_lb = exp(fe$ci.lb),
                   ci_fe_ub = exp(fe$ci.ub), p_fe = fe$pval,
                   strata_k = .k_strat, strata_vars = .s_vars,
                   # MIXED = some but not all. Neither 0 nor k is a mixture: an all-unstratified
                   # pool is one model applied throughout, and so is an all-stratified one.
                   strata_mixed = .k_strat > 0 && .k_strat < m$k)
  } else if (length(est) == 1) {
    # One cohort is not a pooled estimate: there is no between-study variance to identify,
    # and the fixed-effect "sensitivity" is the same single Wald interval. Flagged for the
    # same reason, so a k=1 row never reads as a meta-analytic result.
    e <- est[[1]]
    pooled <- list(HR = e$HR, ci_lb = exp(e$logHR - 1.96 * e$se),
                   ci_ub = exp(e$logHR + 1.96 * e$se), p = e$p, I2 = NA, tau2 = NA, k = 1,
                   tau2_identified = FALSE,
                   # No tau^2 search ran at k=1, so there was nothing to fail to converge.
                   converged_by = "default",
                   HR_fe = e$HR, ci_fe_lb = exp(e$logHR - 1.96 * e$se),
                   ci_fe_ub = exp(e$logHR + 1.96 * e$se), p_fe = e$p,
                   # k=1 cannot be mixed: there is one cohort and one model.
                   strata_k = .k_strat, strata_vars = .s_vars, strata_mixed = FALSE)
  }

  structure(list(per_cohort = per, pooled = pooled, endpoint = endpoint,
                 cohorts = cohorts, max_followup = max_followup), class = "survresult")
}

print.survresult <- function(x, ...) {
  cat(sprintf("get_survival — endpoint = %s (continuous per-SD Cox, HR>1 = worse)\n", x$endpoint))
  # Never let a truncated fit read like a full-follow-up one.
  if (!is.null(x$max_followup))
    cat(sprintf("follow-up truncated at %g months (patients still at risk are censored there)\n",
                x$max_followup))
  cat(strrep("-", 72), "\n")
  cat(sprintf("%-10s %6s %7s %8s %10s %10s %9s  %s\n",
              "cohort", "n", "events", "HR/SD", "p", "logrank", "PH p", "strat-by"))
  for (coh in names(x$per_cohort)) {
    r <- x$per_cohort[[coh]]
    if (isTRUE(r$skipped)) {
      cat(sprintf("%-10s %6d %7d   -- skipped: %s\n", coh, r$n, r$events, r$reason)); next
    }
    # "!" marks a PH violation: the HR on that row is a time-average of an effect that
    # changed over follow-up, not one constant hazard ratio.
    ph <- if (is.na(r$ph_p)) "      n/a" else
      sprintf("%8.2e%s", r$ph_p, if (isTRUE(r$ph_violated)) "!" else " ")
    cat(sprintf("%-10s %6d %7d %8.2f %10.2e %10.2e %9s  %s\n",
                coh, r$n, r$events, r$HR, r$p, r$logrank_p, ph,
                ifelse(r$strata_adjusted, r$strata_var, "-")))
  }
  cat(strrep("-", 72), "\n")
  if (any(vapply(x$per_cohort, function(r) isTRUE(r$ph_violated), logical(1))))
    cat("! = proportional hazards rejected (p <", PH_ALPHA,
        "): that cohort's HR averages a time-varying effect.\n")
  if (!is.null(x$pooled)) {
    p <- x$pooled
    cat(sprintf("POOLED (k=%d): HR/SD = %.2f  [%.2f, %.2f]   p = %.2e   I^2 = %s%%\n",
                p$k, p$HR, p$ci_lb, p$ci_ub, p$p,
                ifelse(is.na(p$I2), "NA", sprintf("%.0f", p$I2))))
    # Below POOL_MIN_K the interval above rests on a between-study variance the data cannot
    # identify. Say so, and show the fixed-effect fit on the same inputs: if the two agree
    # the conclusion is robust to the pooling choice, and if they disagree the reader has
    # been told which of the two they are looking at.
    # A damped fit is the same REML estimate, but the reader should know the default
    # optimizer path could not reach it -- the tau^2 surface here is flat or at a boundary.
    if (identical(p$converged_by, "damped"))
      cat("  ! tau^2 required a damped Fisher step (metafor's default scoring did not\n",
          "   converge). Same REML estimator and estimand; the likelihood surface is flat\n",
          "   or at a boundary for these cohorts.\n", sep = "")
    # A pooled HR that averages stratified and unstratified per-cohort estimates is an
    # average over two different models, and the unstratified arms are not adjusted at all.
    # Said in words here, where the pooled number is printed; carried as `pool_strata_k` by
    # survtable() so it survives into the CSV.
    if (isTRUE(p$strata_mixed))
      cat(sprintf(paste0(
        "  ! MIXED stratification: %d of the %d pooled cohorts fitted with strata(%s); the\n",
        "    other %d: no usable stratifier, fitted UNADJUSTED. The pooled HR averages\n",
        "    estimates made under two different models.\n"),
        p$strata_k, p$k, paste(p$strata_vars, collapse = "/"), p$k - p$strata_k))
    if (isFALSE(p$tau2_identified)) {
      cat(sprintf("  ! tau^2 not identified at k=%d (< %d): the CI above assumes a between-study\n",
                  p$k, POOL_MIN_K))
      cat(sprintf("    variance this many cohorts cannot estimate. fixed-effect, same inputs:\n"))
      cat(sprintf("    HR/SD = %.2f  [%.2f, %.2f]   p = %.2e%s\n",
                  p$HR_fe, p$ci_fe_lb, p$ci_fe_ub, p$p_fe,
                  if (isTRUE(abs(p$HR - p$HR_fe) < 1e-10))
                    "   (identical: tau^2 = 0, choice is moot here)" else ""))
    }
  } else cat("POOLED: n/a (no cohort met the event threshold)\n")
  invisible(x)
}

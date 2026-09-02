# R/brand.R -- the two brand colours, and the one function that derives every shade of them.
# Added 2026-08-30 (step 113).
#
# WHY THIS IS A FILE AND NOT TWO LINES IN app.R
# ---------------------------------------------
# It started in app.R (step 112). Then the Guide's callout and step numbers needed the same
# colours, and R/guide.R is sourced as its own unit with its own test -- so either the hexes
# were typed a second time in a second file, which is exactly what tests/test_brand_colors.R
# exists to prevent, or app.R grew a substitution pass over another file's stylesheet. One
# owner is the third option and the only one where "the brand colour" has a single definition.
#
# THE COLOURS ARE NOT CHOSEN, THEY ARE MEASURED
# ----------------------------------------------
# Each is the mean of the largest quantised ink cluster in logo.png: orange is 48.7% of the
# mark's ink, purple 31.9%. `python3 scripts/make_logo.py` prints both on every run, and
# tests/test_brand_colors.R re-derives them from the image in R -- by an independent path, in
# a different language -- and fails if either constant drifts from what the logo actually is.
# That is what makes "the brand colour" a fact here rather than a preference.
#
# WHAT IS DELIBERATELY *NOT* BRANDED
# -----------------------------------
# The three severity colours in app.R's PAGE_CSS -- .note-good #2a6b2a, .note-warn #a05a00,
# .note-err #b00 -- stay exactly as they are. They are green/amber/red because green, amber
# and red MEAN something to a reader before any word is read, and that meaning is worth more
# than palette consistency. Repainting a warning in the house purple would make the page
# tidier and the warning quieter. Do not "finish the job" by pulling them in.

BRAND_PURPLE <- "#4D127D"
BRAND_ORANGE <- "#F76F08"

# Mixes toward white (p > 0) or toward black (p < 0). Every tint in this project is a call to
# this rather than a typed hex, so that a logo redraw moves the whole palette and not just the
# two colours somebody remembered to update.
brand_tint <- function(hex, p) {
  stopifnot(is.character(hex), length(hex) == 1L, is.numeric(p), length(p) == 1L,
            p >= -1, p <= 1)
  v <- grDevices::col2rgb(hex)[, 1]
  v <- round(if (p >= 0) v * (1 - p) + 255 * p else v * (1 + p))
  sprintf("#%02X%02X%02X", v[[1]], v[[2]], v[[3]])
}

# R/ here is this project's plain function library (loaded via source()/sys.source()
# as needed), not shiny-module support code. Several files in it are CLI scripts with
# a commandArgs(FALSE)-based guard that only works under `Rscript file.R` - don't let
# Shiny's automatic R/ autoload source them blindly on app startup.
options(shiny.autoload.r = FALSE)

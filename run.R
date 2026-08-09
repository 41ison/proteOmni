## proteOmni launcher entry point.
##
## Installs every dependency first (base R only), then starts the Shiny app.
## This is the file the .command / .bat launchers invoke — never call
## shiny::runApp() directly on a fresh R installation, because `shiny::` is
## resolved before proteOmni.r is sourced.

app_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
if (!file.exists(file.path(app_dir, "proteOmni.r"))) app_dir <- getwd()

source(file.path(app_dir, "bootstrap.R"))

shiny::runApp(file.path(app_dir, "proteOmni.r"), launch.browser = TRUE)

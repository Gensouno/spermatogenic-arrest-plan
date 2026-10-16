# Install R packages into the project-local .Rlib.
# The system R library (/Library/Frameworks/R.framework/...) is NOT writable,
# so every install MUST pass lib= pointing inside the project.
# Usage: Rscript envs/install_pkgs.R pkg1 pkg2 ...
args <- commandArgs(trailingOnly = TRUE)
proj <- normalizePath(".", mustWork = TRUE)
lib  <- file.path(proj, ".Rlib")
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"),
        Ncpus = 1, timeout = 3600)

cat("== libPaths ==\n"); print(.libPaths())
for (p in args) {
  t0 <- Sys.time()
  already <- p %in% rownames(installed.packages(lib.loc = lib))
  if (already) {
    cat(sprintf("[RESULT] %-22s SKIP(exists)  %.1fs\n", p, 0)); next
  }
  # Route: try CRAN first, fall back to Bioconductor via BiocManager.
  ok <- tryCatch({
    install.packages(p, lib = lib, quiet = TRUE)
    p %in% rownames(installed.packages(lib.loc = lib))
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!isTRUE(ok)) {
    ok <- tryCatch({
      BiocManager::install(p, lib = lib, ask = FALSE, update = FALSE, quiet = TRUE)
      p %in% rownames(installed.packages(lib.loc = lib))
    }, error = function(e) { cat(sprintf("[ERROR] %s: %s\n", p, conditionMessage(e))); FALSE })
  }
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[RESULT] %-22s %-12s %.1fs\n", p, if (isTRUE(ok)) "OK" else "FAIL", el))
}

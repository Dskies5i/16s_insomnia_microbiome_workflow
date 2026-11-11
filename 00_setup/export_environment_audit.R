# ==============================================================================
# NAME: export_environment_audit.R
# DESCRIPTION: Scans project scripts and active system environment to export 
#              a SINGLE unified audit report inside the 00_setup/ folder.
# USAGE: Rscript 00_setup/export_environment_audit.R
# ==============================================================================

cat("====================================================\n")
cat("  UNIFIED ENVIRONMENT & DEPENDENCY AUDIT EXPORTER\n")
cat("====================================================\n\n")

# Ensure 00_setup directory exists
dir.create("00_setup", showWarnings = FALSE, recursive = TRUE)

# Define the output file path inside 00_setup/
report_file <- "00_setup/environment_audit_report.txt"
con <- file(report_file, open = "write")

# Write Report Header
writeLines("====================================================", con)
writeLines("         UNIFIED ENVIRONMENT AUDIT REPORT           ", con)
writeLines(sprintf("         Generated on: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), con)
writeLines("====================================================\n", con)

# 1. System, Python, and Conda Diagnostics
cat("🔍 Capturing system, Python, and Conda diagnostics...\n")
writeLines("=== 1. SYSTEM & PYTHON DIAGNOSTICS ===", con)
writeLines(sprintf("OS: %s", tryCatch(system("uname -a", intern = TRUE), error = function(e) "N/A")), con)
writeLines(sprintf("Python Version: %s", tryCatch(system("python --version", intern = TRUE), error = function(e) "N/A")), con)
writeLines(sprintf("Python Path (which python): %s", tryCatch(system("which python", intern = TRUE), error = function(e) "N/A")), con)
writeLines("\nConda Environment Info:", con)
conda_info <- tryCatch(system("conda info", intern = TRUE), error = function(e) "conda info unavailable")
writeLines(conda_info, con)
writeLines("\n----------------------------------------------------\n", con)

# 2. Conda Package List
cat("🔍 Exporting active Conda environment packages...\n")
conda_env_name <- Sys.getenv("CONDA_DEFAULT_ENV")
writeLines(sprintf("=== 2. CONDA ENVIRONMENT PACKAGES (%s) ===", ifelse(nzchar(conda_env_name), conda_env_name, "Unknown")), con)
conda_pkgs <- tryCatch(system("conda list", intern = TRUE), error = function(e) "conda list unavailable")
writeLines(conda_pkgs, con)
writeLines("\n----------------------------------------------------\n", con)

# 3. Project R & Bioconductor Packages
cat("🔍 Scanning project scripts for R / Bioconductor dependencies...\n")
files <- list.files(path = c("1_meta_analysis", "utils"), pattern = "\\.R$", full.names = TRUE)

if (length(files) == 0) {
  files <- list.files(pattern = "\\.R$", full.names = TRUE)
}

libs <- unique(unlist(lapply(files, function(f) {
  lines <- readLines(f, warn = FALSE)
  matches <- regmatches(lines, gregexpr("(library|requireNamespace)\\s*\\(\\s*[\"\x27]?([a-zA-Z0-9\\.]+)", lines))
  sub(".*\\(", "", unlist(matches))
})))

libs <- gsub("[\"\x27)]", "", libs)
libs <- intersect(libs, rownames(installed.packages()))

res <- data.frame(
  Package = libs, 
  Version = sapply(libs, function(p) as.character(packageVersion(p)))
)
rownames(res) <- NULL

writeLines("=== 3. PROJECT R & BIOCONDUCTOR PACKAGES ===", con)
writeLines(sprintf("%-35s %-15s", "Package", "Version"), con)
writeLines(paste(rep("-", 50), collapse=""), con)
for (i in seq_len(nrow(res))) {
  writeLines(sprintf("%-35s %-15s", res$Package[i], res$Version[i]), con)
}

close(con)

cat(sprintf("\n  ✔ Unified audit report successfully exported to: %s\n", report_file))
cat("====================================================\n")
cat("  🎉 AUDIT COMPLETE!\n")
cat("====================================================\n")
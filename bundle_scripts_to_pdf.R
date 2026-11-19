# bundle_scripts_to_pdf.R
# Bundles all scripts in 1_meta_analysis into HTML and converts to PDF via pandoc.

input_dir <- "1_meta_analysis"
out_html  <- "1_meta_analysis_scripts.html"
out_pdf   <- "1_meta_analysis_scripts.pdf"
tmp_md    <- tempfile(fileext = ".md")

if (!dir.exists(input_dir)) {
  stop("Directory not found: ", input_dir)
}

if (Sys.which("pandoc") == "") {
  stop("pandoc was not found on PATH. Install pandoc first.")
}

files <- list.files(input_dir, full.names = TRUE, recursive = TRUE)
files <- files[file.info(files)$isdir %in% FALSE]
files <- sort(files)

get_lang <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(
    ext,
    "r"   = "r",
    "sh"  = "bash",
    "rmd" = "rmarkdown",
    "qmd" = "markdown",
    "py"  = "python",
    ""
  )
}

md_lines <- c(
  "# 1_meta_analysis Scripts Bundle",
  "",
  paste0("_Generated: ", Sys.time(), "_"),
  ""
)

for (f in files) {
  rel_path <- sub(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", input_dir), "/?"), "", f)
  lang <- get_lang(f)
  code_lines <- readLines(f, warn = FALSE)
  
  md_lines <- c(
    md_lines,
    paste0("## `", rel_path, "`"),
    "",
    paste0("```", lang),
    code_lines,
    "```",
    "",
    '<div style="page-break-after: always;"></div>',
    ""
  )
}

writeLines(md_lines, tmp_md, useBytes = TRUE)

# 1. Render Markdown -> standalone HTML with syntax highlighting
status_html <- system2(
  "pandoc",
  args = c(
    tmp_md,
    "-f", "markdown",
    "-t", "html",
    "-s",
    "--toc",
    "--highlight-style=tango",
    "-o", out_html
  )
)

if (status_html != 0 || !file.exists(out_html)) {
  stop("Failed to generate HTML with pandoc.")
}
message("Created HTML: ", normalizePath(out_html))

# 2. Convert HTML -> PDF using an available HTML-to-PDF engine (no pdflatex needed)
pdf_converted <- FALSE

if (Sys.which("weasyprint") != "") {
  if (system2("weasyprint", args = c(out_html, out_pdf)) == 0 && file.exists(out_pdf)) {
    pdf_converted <- TRUE
  }
} else if (Sys.which("wkhtmltopdf") != "") {
  if (system2("wkhtmltopdf", args = c(out_html, out_pdf)) == 0 && file.exists(out_pdf)) {
    pdf_converted <- TRUE
  }
} else if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] == "Darwin") {
  # Fallback on macOS: print HTML to PDF via headless Chrome/Edge if installed
  browsers <- c(
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
    "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser"
  )
  browser_bin <- browsers[file.exists(browsers)][1]
  if (!is.na(browser_bin)) {
    abs_html <- paste0("file://", normalizePath(out_html))
    abs_pdf  <- normalizePath(out_pdf, mustWork = FALSE)
    status_pdf <- system2(
      browser_bin,
      args = c("--headless", "--disable-gpu", paste0("--print-to-pdf=", abs_pdf), abs_html),
      stdout = FALSE,
      stderr = FALSE
    )
    if (status_pdf == 0 && file.exists(out_pdf)) {
      pdf_converted <- TRUE
    }
  }
}

if (pdf_converted) {
  message("Created PDF: ", normalizePath(out_pdf))
} else {
  message(
    "No headless PDF engine (weasyprint, wkhtmltopdf, or Chrome/Edge) was found.\n",
    "Open '", out_html, "' in your browser and press Cmd+P -> Save as PDF."
  )
}
#!/usr/bin/env Rscript

# Download the selected GEO metadata and raw IDAT archive.
# Raw files are written under data/raw and are ignored by Git.

options(stringsAsFactors = FALSE)
set.seed(20261003)

project_dir <- normalizePath(".", mustWork = TRUE)
raw_dir <- file.path(project_dir, "data", "raw")
metadata_dir <- file.path(project_dir, "data", "metadata")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

accession <- "GSE212682"
message("Downloading GEO metadata and raw IDAT archive for ", accession)

series_url <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", accession,
  "&targ=self&form=text&view=quick"
)
series_file <- file.path(metadata_dir, paste0(accession, "_series.txt"))
download.file(series_url, series_file, mode = "wb", quiet = FALSE)

archive_file <- file.path(raw_dir, accession, paste0(accession, "_RAW.tar"))
dir.create(dirname(archive_file), recursive = TRUE, showWarnings = FALSE)
archive_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE212nnn/", accession,
  "/suppl/", accession, "_RAW.tar"
)
archive_part <- paste0(archive_file, ".part")
if (!file.exists(archive_file)) {
  status <- system2(
    "curl",
    c("--fail", "--location", "--retry", "10", "--retry-delay", "5",
      "--continue-at", "-", "--output", shQuote(archive_part), archive_url)
  )
  if (!identical(status, 0L)) {
    stop("The raw archive download failed. Re-run this script to resume it.")
  }
  file.rename(archive_part, archive_file)
}

archive_listing <- system2("tar", c("-tf", shQuote(archive_file)), stdout = TRUE, stderr = TRUE)
if (length(archive_listing) == 0L || any(grepl("Truncated input file|Unexpected EOF", archive_listing))) {
  stop("The raw archive is incomplete. Re-run this script after completing a resumable download.")
}

extract_dir <- file.path(raw_dir, accession, "idats")
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
utils::untar(archive_file, exdir = extract_dir)

idat_files <- list.files(extract_dir, pattern = "\\.idat\\.gz$", full.names = TRUE,
                         ignore.case = TRUE)
if (length(idat_files) == 0L) {
  stop("No IDAT files were found in the downloaded archive.")
}
message("Extracted ", length(idat_files), " IDAT files.")
message("Next step: retrieve sample metadata and build the clean sample sheet.")
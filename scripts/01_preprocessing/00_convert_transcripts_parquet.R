# Convert Xenium transcripts.parquet to CSV / CSV.GZ
# --------------------------------------------------
# Consolidated from the three region-specific conversion scripts used in the
# original analysis. This utility does not isolate TMA cores; it only converts
# the region-level Xenium transcript table into a CSV-compatible format.
#
# Usage from repository root:
#   Rscript scripts/01_preprocessing/00_convert_transcripts_parquet.R \
#     /path/to/Region_1/transcripts.parquet
#
# Multiple parquet files may be supplied in one call.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  stop(
    'Supply one or more Xenium transcripts.parquet files.\n',
    'Example: Rscript scripts/01_preprocessing/00_convert_transcripts_parquet.R ',
    '/data/Region_1/transcripts.parquet'
  )
}

chunk_size <- 1e6

convert_transcripts <- function(path, chunk_size = 1e6, gzip_output = TRUE) {
  if (!file.exists(path)) stop('File not found: ', path)
  if (!grepl('\\.parquet$', path, ignore.case = TRUE)) {
    stop('Expected a .parquet file: ', path)
  }

  output_csv <- sub('\\.parquet$', '.csv', path, ignore.case = TRUE)
  if (file.exists(output_csv)) file.remove(output_csv)

  message('Reading: ', path)
  parquet_file <- arrow::read_parquet(path, as_data_frame = FALSE)

  start <- 0
  while (start < parquet_file$num_rows) {
    end <- min(start + chunk_size, parquet_file$num_rows)
    chunk <- as.data.frame(parquet_file$Slice(start, end - start))
    data.table::fwrite(chunk, output_csv, append = start != 0)
    start <- end
    message('  wrote ', format(start, big.mark = ','), ' rows')
  }

  if (gzip_output && requireNamespace('R.utils', quietly = TRUE)) {
    gz_file <- paste0(output_csv, '.gz')
    if (file.exists(gz_file)) file.remove(gz_file)
    R.utils::gzip(output_csv, destname = gz_file, remove = TRUE, overwrite = TRUE)
    message('Saved: ', gz_file)
    return(invisible(gz_file))
  }

  message('Saved: ', output_csv)
  invisible(output_csv)
}

invisible(lapply(args, convert_transcripts, chunk_size = chunk_size))

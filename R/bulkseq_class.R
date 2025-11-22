#' Create a bulkseq object
#'
#' This constructor function bundles count data and metadata into a single
#' S3 object of class 'bulkseq'. It strictly enforces that sample names
#' in the count matrix columns match the metadata row names.
#'
#' @param counts A matrix or data frame of raw counts.
#'   Rows: Genes (Ensembl IDs or Symbols).
#'   Columns: Sample names.
#' @param metadata A data frame of sample information.
#'   Rows: Sample names (must match columns of counts).
#'   Columns: Variables (e.g., treatment, condition, batch).
#'
#' @return An object of class \code{bulkseq} containing sorted \code{counts} and \code{metadata}.
#' @export
#'
#' @examples
#' # 1. Create Dummy Data
#' cnts <- matrix(sample(1:100, 100), ncol = 5)
#' colnames(cnts) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
#' rownames(cnts) <- paste0("Gene", 1:20)
#'
#' # 2. Create Metadata (rownames match counts columns)
#' meta <- data.frame(condition = c("Ctl", "Ctl", "Trt", "Trt", "Trt"))
#' rownames(meta) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
#'
#' # 3. Create Object
#' bs_obj <- create_bulkseqvis_object(cnts, meta)
create_bulkseqvis_object <- function(counts, metadata) {

  # Ensure inputs are in the correct base format
  counts <- as.matrix(counts)
  metadata <- as.data.frame(metadata)

  # 1. Check for missing names
  if (is.null(colnames(counts))) stop("Error: 'counts' must have column names (sample IDs).")
  if (is.null(rownames(metadata))) stop("Error: 'metadata' must have row names (sample IDs).")

  # 2. Check if sample IDs match exactly
  common_samples <- intersect(colnames(counts), rownames(metadata))

  if (length(common_samples) == 0) {
    stop("Error: No common sample names found between counts columns and metadata rownames.")
  }

  if (length(common_samples) != ncol(counts) || length(common_samples) != nrow(metadata)) {
    warning("Mismatch detected: Subsetting to keep only samples present in BOTH counts and metadata.")
  }

  # 3. Align and Sort Data
  # This ensures column 1 of counts corresponds exactly to row 1 of metadata
  counts <- counts[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  # 4. Create the list structure
  obj <- list(
    counts = counts,
    metadata = metadata
  )

  # 5. Assign the class name
  class(obj) <- "bulkseq"

  return(obj)
}

#' Print method for bulkseq objects
#'
#' @param x A bulkseq object.
#' @param ... Additional arguments.
#' @export
print.bulkseq <- function(x, ...) {
  cat("--- BulkSeqVis Object ---\n")
  cat("Genes:", nrow(x$counts), "\n")
  cat("Samples:", ncol(x$counts), "\n")
  cat("Metadata variables:", paste(names(x$metadata), collapse=", "), "\n")
  cat("Sample names:", paste(head(rownames(x$metadata)), collapse=", "), "...\n")
}

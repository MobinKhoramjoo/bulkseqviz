#' Sample-to-Sample Distance Heatmap
#'
#' Generates a heatmap of sample-to-sample distances using VST-transformed counts.
#' Samples are sorted by the color_by variable and clustering is disabled to ensure
#' groups appear adjacent.
#'
#' @param bs_obj A \code{bulkseq} object created by \code{create_bulkseqvis_object}.
#' @param color_by Character string. The name of the column in metadata to use for coloring.
#' @param min_gene_counts Integer. Genes with total counts below this threshold are excluded. Default 10.
#' @param ... Additional arguments passed to \code{pheatmap::pheatmap}.
#'   Useful for changing fonts (e.g. \code{fontsize}), titles (\code{main}), or sizes (\code{cellwidth}).
#'
#' @return A pheatmap object (invisibly).
#' @export
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst
#' @importFrom SummarizedExperiment assay
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats dist setNames
#' @importFrom grDevices colorRampPalette rainbow
#'
#' @examples
#' \dontrun{
#'   # Standard run
#'   sample_distance_heatmap(my_obj, color_by = "condition")
#'
#'   # Customize using standard pheatmap arguments (via ...)
#'   sample_distance_heatmap(my_obj, color_by = "condition",
#'                           fontsize = 14,
#'                           main = "Custom Title")
#' }
sample_distance_heatmap <- function(bs_obj,
                                    color_by,
                                    min_gene_counts = 10,
                                    ...) {

  # 1. Validate Object
  if (!inherits(bs_obj, "bulkseq")) {
    stop("Input must be a 'bulkseq' object.")
  }

  count_mat <- bs_obj$counts
  metadata <- bs_obj$metadata

  # Ensure count matrix is integer for DESeq2
  storage.mode(count_mat) <- "integer"

  # 2. Validate Metadata Column
  if (!color_by %in% names(metadata)) {
    stop(paste("Error: Column", color_by, "not found in metadata."))
  }

  # 2b. Sort Data by Group (to ensure adjacency)
  # We order metadata and counts based on the factor level of the color_by column
  ord <- order(factor(metadata[[color_by]]))
  metadata <- metadata[ord, , drop = FALSE]
  count_mat <- count_mat[, ord, drop = FALSE]

  # 3. Filter low-count genes
  keep_genes <- rowSums(count_mat) >= min_gene_counts
  if (!any(keep_genes)) stop("No genes pass min_gene_counts filter.")
  count_mat <- count_mat[keep_genes, , drop = FALSE]

  # 4. DESeq2 VST Transformation
  # We suppress messages to keep the console clean
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_mat,
                                        colData = metadata,
                                        design = ~ 1)
  dds_vst <- DESeq2::vst(dds, blind = TRUE)

  # 5. Calculate Distance Matrix
  mat_vst <- SummarizedExperiment::assay(dds_vst)
  sampleDists <- stats::dist(t(mat_vst))
  sampleDistMatrix <- as.matrix(sampleDists)
  storage.mode(sampleDistMatrix) <- "double"
  colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)

  # 6. Prepare Annotation
  annotation_df <- data.frame(
    anno = factor(metadata[[color_by]]),
    row.names = rownames(metadata),
    check.names = FALSE
  )
  colnames(annotation_df) <- color_by

  # 7. Prepare Colors
  group_levels <- levels(annotation_df[[color_by]])
  base_cols <- c('#de425b', '#488f31', '#329db3')

  if (length(group_levels) > length(base_cols)) {
    base_cols <- grDevices::rainbow(length(group_levels))
  }

  group_colors <- stats::setNames(base_cols[seq_along(group_levels)], group_levels)
  ann_colors <- list()
  ann_colors[[color_by]] <- group_colors

  # 8. Generate Heatmap
  # Note: clustering is disabled to respect the sorted group order
  # The '...' argument allows users to pass extras like fontsize, main, etc.
  pheatmap::pheatmap(sampleDistMatrix,
                     annotation_row = annotation_df,
                     annotation_col = annotation_df,
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     col = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(255),
                     annotation_colors = ann_colors,
                     main = "Sample-to-sample distance (Sorted by Group)",
                     ...)
}

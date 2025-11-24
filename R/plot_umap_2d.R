#' 2D UMAP Plot
#'
#' Generates a 2D UMAP plot using VST-transformed counts.
#'
#' @param bs_obj A \code{bulkseq} object created by \code{create_bulkseqvis_object}.
#' @param color_by Character string. Metadata column for point color.
#' @param shape_by Character string. Optional. Metadata column for point shape.
#' @param subset_samples Character vector. Optional. Subset of sample IDs to include.
#' @param min_gene_counts Integer. Keep genes with >= this many total reads. Default 100.
#' @param ntop Integer. Number of most variable genes to use for UMAP. Default 1000.
#' @param umap_seed Integer. Seed for reproducibility. Default 123.
#' @param group_colors Character vector. Optional. Custom colors. Defaults to ggsci::nrc_npg.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst
#' @importFrom SummarizedExperiment assay
#' @importFrom matrixStats rowVars
#' @importFrom umap umap
#' @importFrom paletteer paletteer_d
#'
#' @examples
#' \dontrun{
#'   plot_umap_2d(my_obj, color_by = "condition")
#' }
plot_umap_2d <- function(bs_obj,
                         color_by,
                         shape_by = NULL,
                         subset_samples = NULL,
                         min_gene_counts = 100,
                         ntop = 1000,
                         umap_seed = 123,
                         group_colors = NULL) {

  # 1. Validation
  if (!inherits(bs_obj, "bulkseq")) stop("Input must be a 'bulkseq' object.")

  counts <- bs_obj$counts
  metadata <- bs_obj$metadata

  # 2. Subsetting
  if (!is.null(subset_samples)) {
    keep <- intersect(subset_samples, colnames(counts))
    if(length(keep) < 3) stop("Need at least 3 samples for UMAP after subsetting.")
    counts <- counts[, keep, drop=FALSE]
    metadata <- metadata[keep, , drop=FALSE]
  }

  if (!color_by %in% names(metadata)) stop(paste("Column", color_by, "not found in metadata."))
  if (!is.null(shape_by) && !shape_by %in% names(metadata)) stop(paste("Column", shape_by, "not found."))

  # 3. Filter Genes
  keep_genes <- rowSums(counts) >= min_gene_counts
  counts <- counts[keep_genes, , drop=FALSE]

  # 4. VST Transform
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~1)
  dds_vst <- DESeq2::vst(dds, blind = TRUE)
  norm_counts <- SummarizedExperiment::assay(dds_vst)

  # 5. Select most variable genes
  gene_vars <- matrixStats::rowVars(norm_counts)
  top_genes <- order(gene_vars, decreasing = TRUE)[seq_len(min(ntop, nrow(norm_counts)))]
  expr_top  <- norm_counts[top_genes, , drop = FALSE]

  # 6. Run UMAP
  set.seed(umap_seed)
  umap_input <- t(expr_top)
  umap_result <- umap::umap(umap_input, verbose = FALSE) # Verbose FALSE to keep console clean
  umap_df <- as.data.frame(umap_result$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")

  # 7. Merge metadata for plotting
  # Note: metadata rows match expression columns, so we can cbind directly
  umap_df <- cbind(umap_df, metadata)

  # 8. Handle Colors
  n_groups <- length(unique(as.character(metadata[[color_by]])))
  if (is.null(group_colors)) {
    group_colors <- paletteer::paletteer_d("ggsci::nrc_npg")
    if (n_groups > length(group_colors)) {
      group_colors <- rep(group_colors, length.out = n_groups)
    }
  }

  # 9. Plot
  p <- ggplot2::ggplot(umap_df, ggplot2::aes(x = .data[["UMAP1"]],
                                             y = .data[["UMAP2"]],
                                             color = .data[[color_by]],
                                             shape = if(!is.null(shape_by)) .data[[shape_by]] else NULL)) +
    ggplot2::geom_point(size = 4, alpha = 0.7) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 14, face = "bold"),
                   axis.text = ggplot2::element_text(size = 14),
                   legend.title = ggplot2::element_blank(),
                   legend.position = "right") +
    ggplot2::labs(x = "UMAP1", y = "UMAP2")

  return(p)
}

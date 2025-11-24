#' 2D PCA Plot
#'
#' Generates a 2D PCA plot using VST-transformed counts. Supports custom PC selection
#' and ellipses around groups.
#'
#' @param bs_obj A \code{bulkseq} object created by \code{create_bulkseqvis_object}.
#' @param color_by Character string. Metadata column for point color.
#' @param ellipse_by Character string. Optional. Metadata column for grouping ellipses.
#' @param subset_samples Character vector. Optional. Subset of sample IDs to include.
#' @param min_gene_counts Integer. Keep genes with >= this many total reads. Default 10.
#' @param ntop Integer. Number of most variable genes to calculate PCA on. Default 500.
#' @param pcs Numeric vector of length 2. Which PCs to plot (e.g., c(1, 2)). Default c(1, 2).
#' @param ellipse_level Numeric. Confidence level for ellipses. Default 0.99.
#' @param group_colors Character vector. Optional. Custom colors. Defaults to ggsci::nrc_npg.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst
#' @importFrom SummarizedExperiment assay
#' @importFrom ggforce geom_mark_ellipse
#' @importFrom ggthemes theme_base
#' @importFrom paletteer paletteer_d
#' @importFrom stats prcomp var
#'
#' @examples
#' \dontrun{
#'   plot_pca_2d(my_obj, color_by = "condition")
#'   plot_pca_2d(my_obj, color_by = "condition", pcs = c(3, 4))
#' }
plot_pca_2d <- function(bs_obj,
                        color_by,
                        ellipse_by = NULL,
                        subset_samples = NULL,
                        min_gene_counts = 10,
                        ntop = 500,
                        pcs = c(1, 2),
                        ellipse_level = 0.99,
                        group_colors = NULL) {

  # 1. Validation
  if (!inherits(bs_obj, "bulkseq")) stop("Input must be a 'bulkseq' object.")

  counts <- bs_obj$counts
  metadata <- bs_obj$metadata

  # 2. Subsetting
  if (!is.null(subset_samples)) {
    keep <- intersect(subset_samples, colnames(counts))
    if(length(keep) < 3) stop("Need at least 3 samples for PCA after subsetting.")
    counts <- counts[, keep, drop=FALSE]
    metadata <- metadata[keep, , drop=FALSE]
  }

  if (!color_by %in% names(metadata)) stop(paste("Column", color_by, "not found in metadata."))
  if (!is.null(ellipse_by) && !ellipse_by %in% names(metadata)) stop(paste("Column", ellipse_by, "not found."))

  # 3. Filter Genes
  keep_genes <- rowSums(counts) >= min_gene_counts
  counts <- counts[keep_genes, , drop=FALSE]

  # 4. VST Transform
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~1)
  dds_vst <- DESeq2::vst(dds, blind = TRUE)

  # 5. Manual PCA Calculation
  # We do this manually to support arbitrary PCs (standard DESeq2::plotPCA only gives PC1/2)
  vst_mat <- SummarizedExperiment::assay(dds_vst)
  rv <- apply(vst_mat, 1, stats::var)

  # Select top genes
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- stats::prcomp(t(vst_mat[select, ]))
  percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  # Create plotting data frame
  pca_df <- data.frame(
    PCx = pca$x[, pcs[1]],
    PCy = pca$x[, pcs[2]],
    metadata
  )

  # 6. Colors
  group_var <- if(is.null(ellipse_by)) color_by else ellipse_by
  n_groups <- length(unique(as.character(metadata[[color_by]])))

  if (is.null(group_colors)) {
    group_colors <- paletteer::paletteer_d("ggsci::nrc_npg")
    if (n_groups > length(group_colors)) {
      group_colors <- rep(group_colors, length.out = n_groups)
    }
  }

  # 7. Plot
  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = .data[["PCx"]], y = .data[["PCy"]], color = .data[[color_by]])) +
    ggplot2::geom_point(size = 3, alpha = 0.75) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::labs(
      x = paste0("PC", pcs[1], " (", percentVar[pcs[1]], "%)"),
      y = paste0("PC", pcs[2], " (", percentVar[pcs[2]], "%)"),
      color = color_by
    ) +
    ggthemes::theme_base() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_blank()
    )

  # Add Ellipses
  if (!is.null(ellipse_by)) {
    p <- p + ggforce::geom_mark_ellipse(
      ggplot2::aes(group = .data[[ellipse_by]], fill = .data[[ellipse_by]]),
      alpha = 0.15, show.legend = FALSE
    )
  } else {
    p <- p + ggplot2::stat_ellipse(
      ggplot2::aes(group = .data[[color_by]]),
      type = "t", level = ellipse_level, linewidth = 0.6
    )
  }

  return(p)
}

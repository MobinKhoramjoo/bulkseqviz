#' 2D t-SNE Plot
#'
#' Generates a 2D t-SNE plot using VST-transformed counts.
#'
#' @param bs_obj A \code{bulkseq} object created by \code{create_bulkseqvis_object}.
#' @param color_by Character string. Metadata column for point color.
#' @param shape_by Character string. Optional. Metadata column for point shape.
#' @param subset_samples Character vector. Optional. Subset of sample IDs to include.
#' @param min_gene_counts Integer. Keep genes with >= this many total reads. Default 100.
#' @param ntop Integer. Number of most variable genes to use for t-SNE. Default 1000.
#' @param tsne_perplexity Integer. Perplexity parameter for t-SNE. Default 20.
#' @param tsne_seed Integer. Seed for reproducibility. Default 123.
#' @param group_colors Character vector. Optional. Custom colors. Defaults to ggsci::nrc_npg.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst
#' @importFrom SummarizedExperiment assay
#' @importFrom matrixStats rowVars
#' @importFrom Rtsne Rtsne
#' @importFrom paletteer paletteer_d
#'
#' @examples
#' \dontrun{
#'   plot_tsne_2d(my_obj, color_by = "condition")
#' }
plot_tsne_2d <- function(bs_obj,
                         color_by,
                         shape_by = NULL,
                         subset_samples = NULL,
                         min_gene_counts = 100,
                         ntop = 1000,
                         tsne_perplexity = 20,
                         tsne_seed = 123,
                         group_colors = NULL) {

  # 1. Validation
  if (!inherits(bs_obj, "bulkseq")) stop("Input must be a 'bulkseq' object.")

  counts <- bs_obj$counts
  metadata <- bs_obj$metadata

  # 2. Subsetting
  if (!is.null(subset_samples)) {
    keep <- intersect(subset_samples, colnames(counts))
    if(length(keep) < 3) stop("Need at least 3 samples for t-SNE after subsetting.")
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

  # 6. Run t-SNE
  set.seed(tsne_seed)
  tsne_input <- t(expr_top)
  # Rtsne generally requires no duplicates. We set check_duplicates to FALSE strictly if needed,
  # but usually it's safer to let it fail or filter. Here we assume uniqueness or let Rtsne handle.
  # Perplexity needs to be less than (N-1)/3
  max_perplexity <- floor((nrow(tsne_input) - 1) / 3)
  if (tsne_perplexity > max_perplexity) {
    warning(paste("Perplexity is too large for the number of samples. Reducing to", max_perplexity))
    tsne_perplexity <- max_perplexity
  }

  tsne_out <- Rtsne::Rtsne(tsne_input,
                           dims = 2,
                           perplexity = tsne_perplexity,
                           verbose = FALSE,
                           check_duplicates = FALSE)

  tsne_df <- as.data.frame(tsne_out$Y)
  colnames(tsne_df) <- c("tSNE1", "tSNE2")

  # 7. Merge metadata
  tsne_df <- cbind(tsne_df, metadata)

  # 8. Handle Colors
  n_groups <- length(unique(as.character(metadata[[color_by]])))
  if (is.null(group_colors)) {
    group_colors <- paletteer::paletteer_d("ggsci::nrc_npg")
    if (n_groups > length(group_colors)) {
      group_colors <- rep(group_colors, length.out = n_groups)
    }
  }

  # 9. Plot
  p <- ggplot2::ggplot(tsne_df, ggplot2::aes(x = .data[["tSNE1"]],
                                             y = .data[["tSNE2"]],
                                             color = .data[[color_by]],
                                             shape = if(!is.null(shape_by)) .data[[shape_by]] else NULL)) +
    ggplot2::geom_point(size = 4, alpha = 0.7) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 14, face = "bold"),
                   axis.text = ggplot2::element_text(size = 14),
                   legend.title = ggplot2::element_blank(),
                   legend.position = "right") +
    ggplot2::labs(x = "t-SNE 1", y = "t-SNE 2")

  return(p)
}

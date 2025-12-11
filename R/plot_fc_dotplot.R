#' Dotplot of Log2 Fold Changes
#'
#' Generates a dotplot visualizing log2 fold changes and significance for a specific set of genes
#' across multiple comparisons stored in the bulkseq object.
#'
#' @param bs_obj A \code{bulkseq} object containing DE results from \code{DEG()}.
#' @param genes Character vector. A list of gene symbols or Ensembl IDs to plot.
#' @param padj_cutoff Numeric. Filter to show only points with adjusted p-value < cutoff. Default 0.05.
#' @param id_col Character string. Optional. Explicitly specify column to search ("gene_name" or "gene_id").
#'   If NULL (default), searches both.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom dplyr filter select mutate bind_rows arrange case_when
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   my_genes <- c("GAPDH", "ACTB", "TP53")
#'   plot_fc_dotplot(bs_obj, genes = my_genes)
#' }
plot_fc_dotplot <- function(bs_obj, genes, padj_cutoff = 0.05, id_col = NULL) {

  # 1. Validation
  if (!inherits(bs_obj, "bulkseq")) stop("Input must be a 'bulkseq' object.")
  if (length(bs_obj$DE_results) == 0) stop("No DE results found. Run DEG() first.")

  # 2. Extract Data
  res_list <- lapply(names(bs_obj$DE_results), function(comp_name) {
    df <- bs_obj$DE_results[[comp_name]]

    # Determine which rows match the input genes
    if (!is.null(id_col)) {
      if (!id_col %in% colnames(df)) return(NULL)
      match_idx <- df[[id_col]] %in% genes
    } else {
      # Search both ID columns if available
      match_id <- if("gene_id" %in% colnames(df)) df$gene_id %in% genes else FALSE
      match_name <- if("gene_name" %in% colnames(df)) df$gene_name %in% genes else FALSE
      match_idx <- match_id | match_name
    }

    if (sum(match_idx) == 0) return(NULL)

    df_sub <- df[match_idx, ]

    # Prepare for plotting
    # Ensure gene_name column exists for the X-axis label
    label_col <- if("gene_name" %in% colnames(df_sub)) "gene_name" else "gene_id"

    df_sub %>%
      dplyr::select(gene_label = .data[[label_col]], .data$log2FoldChange, .data$padj) %>%
      dplyr::mutate(source = comp_name)
  })

  # Remove NULLs and combine
  res_list <- res_list[!sapply(res_list, is.null)]
  if (length(res_list) == 0) stop("None of the provided genes were found in the results.")

  plot_df <- dplyr::bind_rows(res_list)

  # 3. Transform Data
  plot_df <- plot_df %>%
    dplyr::filter(.data$padj <= padj_cutoff) %>%
    dplyr::mutate(
      padj_trans = -log10(.data$padj + 1e-10), # Avoid log(0)
      padj_trans = ifelse(is.infinite(.data$padj_trans), NA, .data$padj_trans)
    ) %>%
    dplyr::arrange(.data$gene_label)

  if (nrow(plot_df) == 0) stop(paste("No significant genes found with padj <", padj_cutoff))

  # 4. Plot
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$gene_label, y = .data$source)) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$padj_trans, color = .data$log2FoldChange)
    ) +
    ggplot2::scale_color_gradient2(
      low = "#5c6bc0", mid = "white", high = "#C05C6B", midpoint = 0,
      name = "log2FC"
    ) +
    ggplot2::scale_size_continuous(name = "-log10(padj)") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = 'bold', size = 12),
      axis.text.y = ggplot2::element_text(face = 'bold', size = 12),
      text = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(linewidth = 1.2)
    )

  return(p)
}

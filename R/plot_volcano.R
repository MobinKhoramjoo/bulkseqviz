#' Generate Volcano Plot
#'
#' Creates a publication-ready volcano plot from differential expression results.
#'
#' @param bs_obj A \code{bulkseq} object containing DE results from \code{DEG()}.
#' @param comparison_id Character string. The name of the comparison to plot.
#' @param fc_thresh Numeric. The Fold Change threshold (raw scale).
#'    Example: Enter 2 to filter for >2-fold change. Default 2.
#' @param padj_thresh Numeric. Adjusted p-value (FDR) threshold. Default 0.05.
#' @param repel_force Numeric. Force of repulsion for labels. Default 1.5.
#' @param repel_max_overlaps Integer. Max label overlaps allowed. Default 15.
#' @param top_labels Integer. Number of top significant genes to label. Default 0 (no labels).
#' @param plot_title Character string. Optional title. Defaults to comparison_id.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom dplyr select mutate arrange case_when if_else row_number
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggthemes theme_base
#' @importFrom rlang .data
plot_volcano <- function(bs_obj,
                         comparison_id,
                         fc_thresh = 2,
                         padj_thresh = 0.05,
                         repel_force = 1.5,
                         repel_max_overlaps = 15,
                         top_labels = 0,
                         plot_title = NULL) {

  # 1. Validation
  if (!inherits(bs_obj, "bulkseq")) stop("Input must be a 'bulkseq' object.")
  if (is.null(bs_obj$DE_results[[comparison_id]])) {
    stop(paste("Comparison", comparison_id, "not found in bs_obj$DE_results."))
  }

  # 2. Extract Data
  df <- bs_obj$DE_results[[comparison_id]]

  # Ensure gene_name exists for labeling, fallback to gene_id if missing
  if (!"gene_name" %in% colnames(df)) {
    df$gene_name <- df$gene_id
  } else {
    df$gene_name <- dplyr::if_else(is.na(df$gene_name), df$gene_id, df$gene_name)
  }

  # 3. Prepare Plot Data
  if (is.null(plot_title)) plot_title <- comparison_id

  # We still calculate this variable just for the vertical lines plot
  log2_limit <- log2(fc_thresh)

  df <- df %>%
    dplyr::select(.data$gene_id, .data$log2FoldChange, .data$padj, .data$gene_name) %>%
    dplyr::mutate(
      padj = ifelse(is.na(.data$padj), 1, .data$padj),
      log10p = -log10(.data$padj),
      category = dplyr::case_when(
        # EXACT LOGIC REQUESTED:
        .data$log2FoldChange >= log2(fc_thresh) & .data$padj < padj_thresh ~ "Up",
        .data$log2FoldChange <= -log2(fc_thresh) & .data$padj < padj_thresh ~ "Down",
        TRUE ~ "NS"
      ),
      label = ifelse(.data$category != "NS", .data$gene_name, NA_character_)
    )

  # 4. Filter Labels
  if (top_labels > 0) {
    df <- df %>%
      dplyr::arrange(.data$padj) %>%
      dplyr::mutate(label = dplyr::if_else(dplyr::row_number() <= top_labels & .data$category != "NS",
                                           .data$label, NA_character_))
  } else {
    df$label <- NA_character_
  }

  # 5. Summary Counts
  up_fdr <- sum(df$category == "Up", na.rm = TRUE)
  down_fdr <- sum(df$category == "Down", na.rm = TRUE)

  # 6. Axes Limits
  xmax <- max(abs(df$log2FoldChange), na.rm = TRUE)
  ymax <- max(df$log10p, na.rm = TRUE)

  # 7. Plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$log2FoldChange, y = .data$log10p)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$log2FoldChange, size = .data$log10p), alpha = 0.7) +

    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$label),
      max.overlaps = repel_max_overlaps,
      force = repel_force,
      size = 4,
      min.segment.length = 0,
      na.rm = TRUE
    ) +

    ggplot2::scale_color_gradient2(
      low = "#2676b8", mid = "grey80", high = "#d62728", midpoint = 0,
      name = expression(log[2]*"FC")
    ) +

    ggplot2::scale_size_continuous(
      range = c(1, 6),
      name = expression(-log[10]*"(FDR)")
    ) +

    ggplot2::geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = c(-log2(fc_thresh), log2(fc_thresh)), linetype = "dashed", size = 0.5) +
    ggplot2::annotate(
      "label", x = -Inf, y = Inf, label = paste0("Down: ", down_fdr),
      hjust = -0.1, vjust = 1.1, size = 4.2, fontface = "bold",
      fill = "#e8eef7", color = "#2676b8", label.size = 0.2
    ) +

    ggplot2::annotate(
      "label", x = Inf, y = Inf, label = paste0("Up: ", up_fdr),
      hjust = 1.1, vjust = 1.1, size = 4.2, fontface = "bold",
      fill = "#f8e6ea", color = "#d62728", label.size = 0.2
    ) +

    ggplot2::labs(
      title = plot_title,
      x = expression(log[2]*"(FC)"),
      y = expression(-log[10]*"(FDR)")
    ) +

    ggthemes::theme_base() +

    ggplot2::theme(
      text = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 14),
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(0.25, 0.75, 0.25, 0.25, "cm"),
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.grid = ggplot2::element_blank()
    ) +

    ggplot2::coord_cartesian(clip = "off",
                             xlim = c(-xmax*1.1, xmax*1.1),
                             ylim = c(0, ymax*1.15))

  return(p)
}

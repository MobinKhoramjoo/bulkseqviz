#' Compare Log2 Fold Changes
#'
#' Generates a scatter plot comparing log2 fold changes from two different comparisons.
#'
#' @param bs_obj A \code{bulkseq} object containing DE results.
#' @param name1 Character string. The name of the first comparison (must exist in \code{bs_obj$DE_results}).
#' @param name2 Character string. The name of the second comparison (must exist in \code{bs_obj$DE_results}).
#' @param fc_cutoff Numeric. Log2 Fold Change cutoff for quadrants (default 1.5).
#' @param padj_cutoff Numeric. P-adj cutoff for significance (default 0.05).
#'
#' @return A gtable object (plot grid) invisibly. The plot is drawn to the current device.
#' @export
#' @import ggplot2
#' @importFrom dplyr select mutate filter arrange left_join case_when
#' @importFrom tibble tibble deframe
#' @importFrom grid unit textGrob pointsGrob gpar convertWidth grobHeight grobWidth grid.newpage grid.draw
#' @importFrom gtable gtable gtable_add_grob
#' @importFrom stats cor.test complete.cases na.omit
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   plot_fcvsfc(my_obj, name1 = "TreatA_vs_Ctrl", name2 = "TreatB_vs_Ctrl")
#' }
plot_fcvsfc <- function(bs_obj, name1, name2, fc_cutoff = 1.5, padj_cutoff = 0.05) {

  # 1. Validation
  if (!inherits(bs_obj, "bulkseq")) stop("Input must be a 'bulkseq' object.")
  if (is.null(bs_obj$DE_results[[name1]])) stop(paste("Comparison", name1, "not found."))
  if (is.null(bs_obj$DE_results[[name2]])) stop(paste("Comparison", name2, "not found."))

  df1 <- bs_obj$DE_results[[name1]]
  df2 <- bs_obj$DE_results[[name2]]

  # Helpers (defined locally to avoid export)
  safe_count <- function(x, key) { v <- unname(x[key]); if (is.null(v) || is.na(v)) 0L else as.integer(v) }
  fmt_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 1e-4) formatC(p, format = "e", digits = 2)
    else if (p < 0.001) sprintf("%.4f", p)
    else if (p < 0.01) sprintf("%.3f", p)
    else sprintf("%.2f", p)
  }

  # Labels & Colors
  lab_sig_both <- "Sig in both"
  lab_sig_1    <- paste0("Sig only ", name1)
  lab_sig_2    <- paste0("Sig only ", name2)
  lab_nonsig   <- "non-Sig in both"
  legend_labels <- c(lab_sig_both, lab_sig_1, lab_sig_2)
  legend_cols   <- c("#3176CE", "#DD9BF2", "#8CE543")  # blue, lilac, green
  names(legend_cols) <- legend_labels

  # Process Inputs
  process_df <- function(df) {
    df %>%
      dplyr::select(.data$gene_id, .data$log2FoldChange, .data$pvalue, .data$padj) %>%
      dplyr::mutate(
        pvalue = ifelse(is.na(.data$pvalue), 1, .data$pvalue),
        padj   = ifelse(is.na(.data$padj),   1, .data$padj),
        category = dplyr::case_when(
          .data$log2FoldChange > log2(fc_cutoff) & .data$padj < padj_cutoff ~ "Upregulated",
          .data$log2FoldChange < -log2(fc_cutoff) & .data$padj < padj_cutoff ~ "Downregulated",
          TRUE ~ "Non-significant"
        )
      )
  }

  d1 <- process_df(df1)
  d2 <- process_df(df2)

  # Inner Join on shared genes
  FC <- dplyr::inner_join(d1, d2, by = "gene_id", suffix = c("_1", "_2"))

  # Cross-category
  FC <- FC %>%
    dplyr::mutate(
      Category = dplyr::case_when(
        .data$category_1 != "Non-significant" & .data$category_2 != "Non-significant" ~ lab_sig_both,
        .data$category_1 != "Non-significant" & .data$category_2 == "Non-significant" ~ lab_sig_1,
        .data$category_1 == "Non-significant" & .data$category_2 != "Non-significant" ~ lab_sig_2,
        TRUE ~ lab_nonsig
      )
    )

  # Quadrant counts
  q_counts <- FC %>%
    dplyr::filter(.data$padj_1 < padj_cutoff | .data$padj_2 < padj_cutoff) %>%
    dplyr::mutate(quadrant = dplyr::case_when(
      .data$log2FoldChange_1 > log2(fc_cutoff) & .data$log2FoldChange_2 > log2(fc_cutoff) ~ "Q1",
      .data$log2FoldChange_1 < -log2(fc_cutoff) & .data$log2FoldChange_2 > log2(fc_cutoff) ~ "Q2",
      .data$log2FoldChange_1 < -log2(fc_cutoff) & .data$log2FoldChange_2 < -log2(fc_cutoff) ~ "Q3",
      .data$log2FoldChange_1 > log2(fc_cutoff) & .data$log2FoldChange_2 < -log2(fc_cutoff) ~ "Q4",
      TRUE ~ "Other")) %>%
    dplyr::count(.data$quadrant) %>%
    tibble::deframe()

  # Axes
  M <- max(abs(c(FC$log2FoldChange_1, FC$log2FoldChange_2)), na.rm = TRUE)
  pad <- max(0.5, 0.05 * M)
  xlim <- c(-M - pad, M + pad); ylim <- c(-M - pad, M + pad)
  x_max <- 0.8 * M; y_max <- 0.8 * M

  # Correlation
  corr_df <- FC %>%
    dplyr::filter(.data$Category != lab_nonsig) %>%
    dplyr::transmute(x = .data$log2FoldChange_1, y = .data$log2FoldChange_2) %>%
    stats::na.omit()

  if (nrow(corr_df) >= 3) {
    ct <- suppressWarnings(stats::cor.test(corr_df$x, corr_df$y, method = "pearson"))
    r  <- unname(ct$estimate); r2 <- r^2; p <- ct$p.value
    corr_lines <- paste0("r = ", sprintf("%.3f", r), "\n",
                         "r\u00B2 = ", sprintf("%.3f", r2), "\n",
                         "p = ",  fmt_p(p))
  } else {
    corr_lines <- "r = NA\nr\u00B2 = NA\np = NA"
  }

  # Visibility
  vis_df <- FC %>%
    dplyr::filter(.data$Category %in% legend_labels) %>%
    dplyr::mutate(Alpha = ifelse(.data$Category == lab_sig_both, 0.75, 0.25),
                  Category = factor(.data$Category, levels = legend_labels))

  hid_df <- FC %>% dplyr::filter(.data$Category == lab_nonsig)

  # Main Plot
  p_main <- ggplot2::ggplot() +
    ggplot2::geom_point(data = hid_df,
                        ggplot2::aes(x = .data$log2FoldChange_1, y = .data$log2FoldChange_2),
                        alpha = 0, show.legend = FALSE) +
    ggplot2::geom_point(data = vis_df,
                        ggplot2::aes(x = .data$log2FoldChange_1, y = .data$log2FoldChange_2,
                                     color = .data$Category, alpha = .data$Alpha),
                        size = 3, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = legend_cols, breaks = legend_labels, limits = legend_labels, drop = FALSE) +
    ggplot2::scale_alpha_identity() +
    ggplot2::geom_smooth(data = corr_df, ggplot2::aes(x = .data$x, y = .data$y),
                         method = "lm", se = FALSE, linetype = 1, color = "black", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = c(-log2(fc_cutoff),log2(fc_cutoff)), linetype = 4, color = "black", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = c(-log2(fc_cutoff),log2(fc_cutoff)), linetype = 4, color = "black", linewidth = 0.8) +
    ggplot2::annotate("text", x =  x_max+0.5, y =  y_max+0.5, label = safe_count(q_counts,"Q1"), fontface = "bold", size = 5) +
    ggplot2::annotate("text", x = -x_max-0.5, y =  y_max+0.5, label = safe_count(q_counts,"Q2"), fontface = "bold", size = 5) +
    ggplot2::annotate("text", x = -x_max-0.5, y = -y_max-0.5, label = safe_count(q_counts,"Q3"), fontface = "bold", size = 5) +
    ggplot2::annotate("text", x =  x_max+0.5, y = -y_max-0.5, label = safe_count(q_counts,"Q4"), fontface = "bold", size = 5) +
    ggplot2::labs(
      x = paste0("log2FC (", name1, ")"),
      y = paste0("log2FC (", name2, ")"),
      title = paste0("(", name1, ") vs (", name2, ")")
    ) +
    ggplot2::coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "right",
      text = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(1.5, 1.5, 1.5, 1.5)
    )

  # Gtable Construction
  point_col_w <- grid::unit(6, "mm")
  txt_grobs <- lapply(legend_labels, function(lbl) grid::textGrob(lbl, x = 0, hjust = 0, gp = grid::gpar(fontface = "bold", cex = 0.9)))
  txt_w <- do.call(grid::unit.pmax, lapply(txt_grobs, grid::grobWidth))

  leg_tbl <- gtable::gtable(widths = grid::unit.c(point_col_w + grid::unit(3, "pt"), txt_w + grid::unit(3, "pt")),
                            heights = grid::unit(rep(1, length(legend_labels)), "lines"))

  for (i in seq_along(legend_labels)) {
    pt <- grid::pointsGrob(x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"),
                           pch = 16, size = grid::unit(4, "mm"),
                           gp = grid::gpar(col = legend_cols[i], fill = legend_cols[i]))
    leg_tbl <- gtable::gtable_add_grob(leg_tbl, pt, t = i, l = 1, b = i, r = 1, name = paste0("pt", i))
    leg_tbl <- gtable::gtable_add_grob(leg_tbl, txt_grobs[[i]], t = i, l = 2, b = i, r = 2, name = paste0("tx", i))
  }

  corr_grob <- grid::textGrob(corr_lines, x = 0, y = 1, hjust = 0, vjust = 1,
                              gp = grid::gpar(fontface = "bold", cex = 0.9))

  right_col <- gtable::gtable(widths = grid::unit(1, "null"),
                              heights = grid::unit.c(grid::grobHeight(leg_tbl), grid::grobHeight(corr_grob) + grid::unit(6, "pt")))
  right_col <- gtable::gtable_add_grob(right_col, leg_tbl,  1, 1, 1, 1, name = "legend")
  right_col <- gtable::gtable_add_grob(right_col, corr_grob, 2, 1, 2, 1, name = "corrtext")

  plot_grob <- ggplot2::ggplotGrob(p_main)
  right_w <- max(grid::convertWidth(grid::grobWidth(leg_tbl), "pt", valueOnly = TRUE),
                 grid::convertWidth(grid::grobWidth(corr_grob), "pt", valueOnly = TRUE)) + 2

  final <- gtable::gtable(widths = grid::unit.c(grid::unit(1, "null"), grid::unit(right_w, "pt")),
                          heights = grid::unit(1, "null"))
  final <- gtable::gtable_add_grob(final, plot_grob, 1, 1, 1, 1, name = "plot")
  final <- gtable::gtable_add_grob(final, right_col, 1, 2, 1, 2, name = "rightcol")

  # Draw the plot
  grid::grid.newpage()
  grid::grid.draw(final)

  return(invisible(final))
}

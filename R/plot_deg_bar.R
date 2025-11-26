#' Barplot of DEGs Summary
#'
#' Generates a bar plot summarizing the number of upregulated and downregulated genes
#' across different comparisons stored in the bulkseq object.
#'
#' @param bs_obj A \code{bulkseq} object containing DE results.
#' @param method Character string. P-value column to use: "padj" or "pvalue". Default "padj".
#' @param cutoff Numeric. Significance threshold. Default 0.05.
#' @param fc_cutoff Numeric. Log2 Fold Change threshold (absolute value). Default 1 (fold change 2).
#' @param custom_colors Named character vector. Colors for "Upregulated", "Downregulated", and "All".
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom dplyr mutate filter summarise n
#' @importFrom tidyr pivot_longer
#' @importFrom purrr imap_dfr
#' @importFrom rlang sym !! .data
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#'   plot_deg_bar(my_obj)
#' }
plot_deg_bar <- function(bs_obj,
                         method = "padj",
                         cutoff = 0.05,
                         fc_cutoff = 1,
                         custom_colors = NULL) {

  # Validation
  if (!inherits(bs_obj, "bulkseq")) stop("Input must be a 'bulkseq' object.")
  if (length(bs_obj$DE_results) == 0) stop("No DE results found. Run DEG() first.")
  if (!method %in% c("padj", "pvalue")) stop("Method must be 'padj' or 'pvalue'.")

  # Prepare Data
  deg_list <- bs_obj$DE_results

  # Summarize
  deg_summary <- purrr::imap_dfr(deg_list, function(df, name) {

    # Ensure column exists
    if(!method %in% names(df)) return(NULL)

    # Handle NAs
    p_vals <- df[[method]]
    p_vals[is.na(p_vals)] <- 1

    log2fc <- df[["log2FoldChange"]]

    lfc_thresh <- log2(2^fc_cutoff) # Ensure we treat input as log2fc threshold directly or linear?
    # Standard convention: fc_cutoff=1 usually means log2FC=1 (2-fold).
    # The original user code had log2(fc_cutoff), implying linear input.
    # Let's assume input is LOG2FC for consistency with other package functions.
    # If user wants 2-fold, they input 1.

    n_up <- sum(log2fc >= fc_cutoff & p_vals < cutoff, na.rm = TRUE)
    n_down <- sum(log2fc <= -fc_cutoff & p_vals < cutoff, na.rm = TRUE)
    n_all <- n_up + n_down

    tibble::tibble(
      Dataset = name,
      Upregulated = n_up,
      Downregulated = n_down,
      All = n_all
    )
  })

  if(nrow(deg_summary) == 0) stop("Could not calculate summary. Check DE result structures.")

  # Long format
  deg_long <- tidyr::pivot_longer(deg_summary,
                                  cols = c("Upregulated", "Downregulated", "All"),
                                  names_to = "Category",
                                  values_to = "Count")

  # Factor levels for plotting order
  deg_long$Category <- factor(deg_long$Category, levels = c("All", "Upregulated", "Downregulated"))

  # Colors
  if (is.null(custom_colors)) {
    custom_colors <- c(
      "Upregulated" = "#C05C6B",
      "Downregulated" = "#5c6bc0",
      "All" = "gray70"
    )
  }

  # Plot
  p <- ggplot2::ggplot(deg_long, ggplot2::aes(x = .data$Dataset, y = .data$Count, fill = .data$Category)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Count),
                       position = ggplot2::position_dodge(width = 0.9),
                       vjust = -0.5, size = 4) +
    ggplot2::scale_fill_manual(values = custom_colors) +
    ggplot2::labs(
      title = paste0("Differentially Expressed Genes\n(", method, " < ", cutoff, ", |log2FC| > ", fc_cutoff, ")"),
      y = "Number of Genes",
      x = ""
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = 'bold'),
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.grid = ggplot2::element_blank(),
      title = ggplot2::element_text(face = 'bold'),
      axis.text.y = ggplot2::element_text(face = 'bold')
    )

  return(p)
}

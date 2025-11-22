#' Overall Count Boxplot
#'
#' Generates a boxplot of log10(counts + 1) for all samples, colored by a metadata variable.
#'
#' @param bs_obj A \code{bulkseq} object created by \code{create_bulkseq}.
#' @param color_by Character string. The name of the column in metadata to use for coloring.
#' @param group_order Character vector. Optional. Specify the order of levels for the \code{color_by} variable.
#'   Defaults to alphabetical order if NULL.
#' @param group_colors Character vector. Optional. Custom colors for the groups.
#'   Defaults to \code{paletteer_d("ggsci::nrc_npg")} if NULL.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom dplyr left_join distinct arrange pull
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom paletteer paletteer_d
#'
#' @examples
#' \dontrun{
#'   overall_count_boxplot(my_obj, color_by = "condition")
#' }
overall_count_boxplot <- function(bs_obj,
                                  color_by,
                                  group_order = NULL,
                                  group_colors = NULL) {

  # 1. Validate Object
  if (!inherits(bs_obj, "bulkseq")) {
    stop("Input must be a 'bulkseq' object.")
  }

  count_mat <- bs_obj$counts
  metadata <- bs_obj$metadata

  # 2. Validate Metadata Column
  if (!color_by %in% names(metadata)) {
    stop(paste("Error: Column", color_by, "not found in metadata."))
  }

  # 3. Handle Defaults
  # Default Order: Alphabetical
  if (is.null(group_order)) {
    group_order <- sort(unique(as.character(metadata[[color_by]])))
  }

  # Default Colors: ggsci::nrc_npg
  if (is.null(group_colors)) {
    group_colors <- paletteer::paletteer_d("ggsci::nrc_npg")

    # Handle case where we have more groups than colors in the palette
    if (length(group_order) > length(group_colors)) {
      warning("More groups than colors in default palette. Recycling colors.")
      group_colors <- rep(group_colors, length.out = length(group_order))
    }
  }

  # 4. Prepare Data
  # Set factor levels for plotting order
  metadata[[color_by]] <- factor(metadata[[color_by]], levels = group_order)

  log_counts <- log10(count_mat + 1)

  # Reshape
  long_df <- as.data.frame(log_counts) %>%
    tibble::rownames_to_column(var = "gene_id") %>%
    tidyr::pivot_longer(-gene_id, names_to = "Sample", values_to = "log10_count") %>%
    # Join using Sample as key. We convert rownames to a column named "Sample" to match pivot_longer
    dplyr::left_join(metadata %>% tibble::rownames_to_column("Sample"), by = "Sample")

  # Reorder Samples on X-axis based on group_order
  long_df$Sample <- factor(long_df$Sample,
                           levels = long_df %>%
                             dplyr::distinct(Sample, .data[[color_by]]) %>%
                             dplyr::arrange(factor(.data[[color_by]], levels = group_order)) %>%
                             dplyr::pull(Sample))

  # 5. Plot
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = Sample, y = log10_count, fill = .data[[color_by]])) +
    ggplot2::geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    ggplot2::scale_fill_manual(values = group_colors, drop = FALSE) +
    ggplot2::labs(x = "Samples", y = "log10(Counts + 1)", fill = color_by) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 9),
                   axis.title = ggplot2::element_text(size = 13, face = "bold"),
                   legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 12),
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())

  return(p)
}

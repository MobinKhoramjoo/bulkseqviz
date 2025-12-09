#' Overall Count Boxplot
#'
#' Generates a boxplot of log10(counts + 1) for all samples, colored by a metadata variable.
#' Returns a ggplot object that can be further customized (e.g., adding theme layers).
#'
#' @param bs_obj A \code{bulkseq} object created by \code{create_bulkseqvis_object}.
#' @param color_by Character string. The name of the column in metadata to use for coloring.
#' @param group_order Character vector. Optional. Specify the order of levels for the \code{color_by} variable.
#'   Defaults to alphabetical order if NULL.
#' @param group_colors Character vector. Optional. Custom colors for the groups.
#'   Defaults to \code{paletteer_d("ggsci::nrc_npg")} if NULL.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom dplyr left_join distinct arrange pull %>%
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom paletteer paletteer_d
#'
#' @examples
#' \dontrun{
#'   # You can add ggplot themes to the result:
#'   overall_count_boxplot(my_obj, color_by = "condition") + ggplot2::theme_classic()
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
  if (is.null(group_order)) {
    group_order <- sort(unique(as.character(metadata[[color_by]])))
  }

  if (is.null(group_colors)) {
    group_colors <- paletteer::paletteer_d("ggsci::nrc_npg")
    if (length(group_order) > length(group_colors)) {
      warning("More groups than colors in default palette. Recycling colors.")
      group_colors <- rep(group_colors, length.out = length(group_order))
    }
  }

  # 4. Prepare Data
  metadata[[color_by]] <- factor(metadata[[color_by]], levels = group_order)
  log_counts <- log10(count_mat + 1)

  # We use cols = -1 to exclude the 'gene_id' column without naming it explicitly
  # to avoid R CMD check 'global variable' notes.
  long_df <- as.data.frame(log_counts) %>%
    tibble::rownames_to_column(var = "gene_id") %>%
    tidyr::pivot_longer(cols = -1, names_to = "Sample", values_to = "log10_count") %>%
    dplyr::left_join(metadata %>% tibble::rownames_to_column("Sample"), by = "Sample")

  # Fix factor levels for ordering on X-axis
  long_df$Sample <- factor(long_df$Sample,
                           levels = long_df %>%
                             dplyr::distinct(.data[["Sample"]], .data[[color_by]]) %>%
                             dplyr::arrange(factor(.data[[color_by]], levels = group_order)) %>%
                             dplyr::pull(.data[["Sample"]]))

  # 5. Plot
  # Use .data[["name"]] to avoid global variable notes in R CMD check
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = .data[["Sample"]],
                                             y = .data[["log10_count"]],
                                             fill = .data[[color_by]])) +
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

#' Gene Expression Boxplot
#'
#' Generates boxplots for specific genes across conditions with statistical comparisons.
#' Automatically fetches gene symbols via BioMart if needed.
#'
#' @param bs_obj A \code{bulkseq} object.
#' @param genes Character vector. Gene IDs or Symbols to plot (must match rownames of counts).
#' @param x_var Character string. Metadata column for x-axis.
#' @param region_var Character string. Metadata column for faceting.
#' @param x_order Character vector. Optional. Order of levels for x-axis.
#' @param region_order Character vector. Optional. Order of levels for facet.
#' @param x_colors Named character vector. Optional. Colors for x-axis groups.
#' @param jitter_points Logical. Add jitter points? Default TRUE.
#' @param normalize_to_baseline Logical. Normalize expression to the median of the first x-level? Default FALSE.
#' @param log2_transform Logical. Log2 transform counts? Default FALSE.
#' @param show_ns Logical. Show non-significant brackets? Default TRUE.
#' @param biomart_dataset Character string. BioMart dataset for symbol mapping (e.g., "hsapiens_gene_ensembl").
#'   Default "hsapiens_gene_ensembl". Set to NULL to skip fetching.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
#' @importFrom dplyr mutate select filter group_by summarise ungroup left_join all_of any_of
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr geom_pwc
#' @importFrom rlang sym !!
#' @importFrom paletteer paletteer_d
#' @importFrom stats median
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom utils head
#'
#' @examples
#' \dontrun{
#'   plot_gene_boxplot(bs_obj, genes = c("ENSG00000139618"),
#'                     x_var = "condition", region_var = "batch")
#' }
plot_gene_boxplot <- function(bs_obj,
                              genes,
                              x_var,
                              region_var,
                              x_order = NULL,
                              region_order = NULL,
                              x_colors = NULL,
                              jitter_points = TRUE,
                              normalize_to_baseline = FALSE,
                              log2_transform = FALSE,
                              show_ns = TRUE,
                              biomart_dataset = "hsapiens_gene_ensembl") {

  # 1. Validation
  if (!inherits(bs_obj, "bulkseq")) stop("Input must be a 'bulkseq' object.")

  count_data <- bs_obj$counts
  metadata <- bs_obj$metadata

  if (!all(genes %in% rownames(count_data))) stop("Some 'genes' are not found in count matrix rownames.")
  if (!x_var %in% names(metadata)) stop(paste("Column", x_var, "not found in metadata."))
  if (!region_var %in% names(metadata)) stop(paste("Column", region_var, "not found in metadata."))

  # Capture symbols
  x_sym <- rlang::sym(x_var)
  region_sym <- rlang::sym(region_var)

  # 2. Normalization (DESeq2 size factors)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~1)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  # 3. Subset and Format
  sub_counts <- norm_counts[genes, , drop = FALSE] %>% as.data.frame()
  sub_counts$gene_id <- rownames(sub_counts)
  sub_counts$plot_name <- sub_counts$gene_id # Default fallback

  # 4. BioMart Annotation (Fetch Symbol)
  if (!is.null(biomart_dataset)) {
    # Heuristic: Check if first gene looks like Ensembl ID (starts with ENS)
    first_id <- utils::head(sub_counts$gene_id, 1)
    is_ensembl <- grepl("^ENS", first_id)

    if (is_ensembl) {
      tryCatch({
        mart <- biomaRt::useEnsembl("ensembl", dataset = biomart_dataset)
        map <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                              filters = "ensembl_gene_id",
                              values = sub_counts$gene_id,
                              mart = mart)

        # Merge mapping
        if(nrow(map) > 0) {
          colnames(map) <- c("gene_id", "symbol")
          sub_counts <- dplyr::left_join(sub_counts, map, by = "gene_id")
          # Use symbol if available, otherwise stick to ID
          sub_counts$plot_name <- ifelse(!is.na(sub_counts$symbol) & sub_counts$symbol != "",
                                         sub_counts$symbol, sub_counts$gene_id)
        }
      }, error = function(e) {
        warning("BioMart fetch failed. Using IDs as names.")
      })
    }
    # If not Ensembl ID (e.g. user already has symbols as rownames), plot_name stays as gene_id.
  }

  # Fix: Use any_of() to safely exclude columns that might not exist (e.g. symbol)
  long_df <- tidyr::pivot_longer(sub_counts,
                                 cols = -dplyr::any_of(c("gene_id", "plot_name", "symbol")),
                                 names_to = "sample",
                                 values_to = "value")

  # Add metadata
  meta_sub <- metadata %>%
    dplyr::mutate(sample = rownames(metadata)) %>%
    dplyr::select(dplyr::all_of(c("sample", x_var, region_var)))

  dat <- dplyr::left_join(long_df, meta_sub, by = "sample") %>%
    dplyr::filter(!is.na(!!region_sym), !is.na(!!x_sym))

  # 5. Factor Ordering
  if (!is.null(x_order)) dat[[x_var]] <- factor(dat[[x_var]], levels = x_order)
  if (!is.null(region_order)) dat[[region_var]] <- factor(dat[[region_var]], levels = region_order)

  # 6. Baseline Normalization
  if (normalize_to_baseline) {
    base_level <- levels(factor(dat[[x_var]]))[1]
    dat <- dat %>%
      dplyr::group_by(!!region_sym, .data$plot_name) %>%
      dplyr::mutate(
        base_med = stats::median(.data$value[!!x_sym == base_level], na.rm = TRUE),
        y = .data$value / .data$base_med
      ) %>%
      dplyr::ungroup()
  } else {
    dat$y <- dat$value
  }

  # 7. Log Transform
  if (log2_transform) dat$y <- log2(dat$y + 1)

  # 8. Colors
  n_x <- length(unique(dat[[x_var]]))
  if (is.null(x_colors)) {
    x_colors <- paletteer::paletteer_d("ggsci::nrc_npg", n = n_x)
  }

  # 9. Plot
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = !!x_sym, y = .data$y, fill = !!x_sym)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.72, alpha = 0.25, colour = "black")

  if (jitter_points) {
    p <- p + ggplot2::geom_jitter(ggplot2::aes(color = !!x_sym), width = 0.14, size = 2.2, alpha = 0.85, show.legend = FALSE)
  }

  p <- p +
    ggpubr::geom_pwc(
      ggplot2::aes(group = !!x_sym),
      method = "wilcox.test",
      p.adjust.method = "fdr",
      label = "p.adj.signif",
      hide.ns = !show_ns,
      tip.length = 0.01,
      group.by = region_var,
      label.size = 6
    ) +
    # Use plot_name for facet header
    ggplot2::facet_grid(rows = ggplot2::vars(.data$plot_name), cols = ggplot2::vars(!!region_sym), scales = "free_y") +
    ggplot2::scale_fill_manual(values = x_colors, name = x_var) +
    ggplot2::scale_color_manual(values = x_colors, guide = "none") +
    ggplot2::labs(x = NULL, y = if(normalize_to_baseline) "Fold Change vs Baseline" else "Normalized Expression") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.y = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(size = 11, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )

  if(normalize_to_baseline) {
    p <- p + ggplot2::geom_hline(yintercept = 1, linetype = 2, linewidth = 0.4, colour = "grey30")
  }

  return(p)
}

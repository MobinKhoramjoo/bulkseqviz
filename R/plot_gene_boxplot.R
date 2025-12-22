#' Gene Expression Boxplot
#'
#' Generates boxplots for specific genes across conditions with statistical comparisons.
#' Automatically fetches gene symbols via BioMart if needed.
#'
#' @param bs_obj A \code{bulkseq} object.
#' @param genes Character vector. Gene IDs or Symbols to plot (must match rownames of counts).
#' @param x_var Character string. Metadata column for x-axis.
#' @param region_var Character string. Optional. Metadata column for faceting. Default NULL.
#' @param x_order Character vector. Optional. Order of levels for x-axis.
#' @param region_order Character vector. Optional. Order of levels for facet.
#' @param x_colors Named character vector. Optional. Colors for x-axis groups.
#' @param jitter_points Logical. Add jitter points? Default TRUE.
#' @param normalize_to_baseline Logical. Normalize expression to the median of the first x-level? Default FALSE.
#' @param log2_transform Logical. Log2 transform counts? Default FALSE.
#' @param show_ns Logical. Show non-significant brackets? Default TRUE.
#' @param biomart_dataset Character string. BioMart dataset for symbol mapping (e.g., "hsapiens_gene_ensembl").
#'    Default "hsapiens_gene_ensembl". Set to NULL to skip fetching.
#'
#' @return A ggplot object.
#' @export
#' @import ggplot2
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts
#' @importFrom dplyr mutate select filter group_by summarise ungroup left_join all_of across
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr geom_pwc
#' @importFrom rlang sym !!
#' @importFrom paletteer paletteer_d
#' @importFrom stats median
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom utils head
plot_gene_boxplot <- function(bs_obj,
                              genes,
                              x_var,
                              region_var = NULL,
                              x_order = NULL,
                              region_order = NULL,
                              x_colors = NULL,
                              jitter_points = TRUE,
                              normalize_to_baseline = FALSE,
                              log2_transform = FALSE,
                              show_ns = TRUE,
                              biomart_dataset = "hsapiens_gene_ensembl") {

  # 1. Input Validation
  if (!inherits(bs_obj, "bulkseq")) stop("Input must be a 'bulkseq' object.")
  count_data <- bs_obj$counts
  metadata <- bs_obj$metadata

  if (!x_var %in% names(metadata)) stop(paste("Column", x_var, "not found in metadata."))
  if (!is.null(region_var) && !region_var %in% names(metadata)) {
    stop(paste("Column", region_var, "not found in metadata."))
  }

  # Capture symbols for plotting
  x_sym <- rlang::sym(x_var)
  region_sym <- if (!is.null(region_var)) rlang::sym(region_var) else NULL

  # 2. Normalization
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~1)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)

  # 3. Gene Name Resolution
  gene_map <- data.frame(input = genes, matrix_id = NA, plot_name = NA, stringsAsFactors = FALSE)

  # Check direct matches first
  direct_match <- genes %in% rownames(norm_counts)
  gene_map$matrix_id[direct_match] <- genes[direct_match]

  # Try BioMart if needed
  if (!is.null(biomart_dataset)) {
    mart <- tryCatch({
      biomaRt::useEnsembl("ensembl", dataset = biomart_dataset)
    }, error = function(e) { warning("BioMart connection failed."); NULL })

    if (!is.null(mart)) {
      first_id <- rownames(norm_counts)[1]
      matrix_is_ens <- grepl("^ENS", first_id)

      missing_mask <- is.na(gene_map$matrix_id)
      if (any(missing_mask) && matrix_is_ens) {
        syms <- gene_map$input[missing_mask]
        map_ids <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                                  filters = "external_gene_name", values = syms, mart = mart)
        idx <- match(gene_map$input[missing_mask], map_ids$external_gene_name)
        gene_map$matrix_id[missing_mask] <- map_ids$ensembl_gene_id[idx]
        gene_map$plot_name[missing_mask] <- gene_map$input[missing_mask]
      }

      ens_mask <- !is.na(gene_map$matrix_id) & grepl("^ENS", gene_map$matrix_id)
      if (any(ens_mask)) {
        ids <- gene_map$matrix_id[ens_mask]
        map_sym <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                                  filters = "ensembl_gene_id", values = ids, mart = mart)
        idx <- match(gene_map$matrix_id[ens_mask], map_sym$ensembl_gene_id)
        found_sym <- map_sym$external_gene_name[idx]
        has_sym <- !is.na(found_sym) & found_sym != ""
        gene_map$plot_name[ens_mask][has_sym] <- found_sym[has_sym]
      }
    }
  }

  # Fallbacks
  found_mask <- !is.na(gene_map$matrix_id)
  gene_map$plot_name[found_mask & is.na(gene_map$plot_name)] <- gene_map$input[found_mask & is.na(gene_map$plot_name)]

  if (any(is.na(gene_map$matrix_id))) {
    stop(paste("Could not find gene(s) in count matrix:", paste(gene_map$input[is.na(gene_map$matrix_id)], collapse=", ")))
  }

  # 4. Subset and Plot prep
  sub_counts <- norm_counts[gene_map$matrix_id, , drop = FALSE] %>% as.data.frame()
  sub_counts$matrix_id <- rownames(sub_counts)
  sub_counts <- dplyr::left_join(sub_counts, unique(gene_map[,c("matrix_id", "plot_name")]), by="matrix_id")

  long_df <- tidyr::pivot_longer(sub_counts, cols = -c("matrix_id", "plot_name"),
                                 names_to = "sample", values_to = "value")

  # Select columns safely
  sel_cols <- c("sample", x_var)
  if(!is.null(region_var)) sel_cols <- c(sel_cols, region_var)

  meta_sub <- metadata %>%
    dplyr::mutate(sample = rownames(metadata)) %>%
    dplyr::select(dplyr::all_of(sel_cols))

  dat <- dplyr::left_join(long_df, meta_sub, by = "sample") %>%
    dplyr::filter(!is.na(!!x_sym))

  if (!is.null(region_var)) {
    dat <- dat %>% dplyr::filter(!is.na(!!region_sym))
  }

  # 5. Factor Ordering
  if (!is.null(x_order)) dat[[x_var]] <- factor(dat[[x_var]], levels = x_order)
  if (!is.null(region_var) && !is.null(region_order)) {
    dat[[region_var]] <- factor(dat[[region_var]], levels = region_order)
  }

  # 6. Baseline Normalization
  if (normalize_to_baseline) {
    base_level <- levels(factor(dat[[x_var]]))[1]

    # Group by region if it exists, otherwise just by gene (plot_name)
    group_cols <- "plot_name"
    if (!is.null(region_var)) group_cols <- c(region_var, "plot_name")

    dat <- dat %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
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

  # Add Stats (Conditionally handle grouping)
  pwc_args <- list(
    ggplot2::aes(group = !!x_sym),
    method = "wilcox.test",
    p.adjust.method = "fdr",
    label = "p.adj.signif",
    hide.ns = !show_ns,
    tip.length = 0.01,
    label.size = 6
  )

  # Only add group.by if region_var exists
  if (!is.null(region_var)) {
    pwc_args$group.by <- region_var
  }

  p <- p + do.call(ggpubr::geom_pwc, pwc_args)

  # Faceting logic
  if (!is.null(region_var)) {
    # Facet grid: Gene rows vs Region columns
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data$plot_name),
                                 cols = ggplot2::vars(!!region_sym),
                                 scales = "free_y")
  } else {
    # Facet wrap: Just genes
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$plot_name),
                                 scales = "free_y")
  }

  p <- p +
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

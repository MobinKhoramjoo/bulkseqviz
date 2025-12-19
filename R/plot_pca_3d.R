#' 3D PCA Plot (Interactive)
#'
#' Generates an interactive 3D PCA plot using plotly. Includes 3D ellipses for groups.
#'
#' @param bs_obj A \code{bulkseq} object created by \code{create_bulkseqvis_object}.
#' @param color_by Character string. Metadata column for point color.
#' @param ntop Integer. Number of most variable genes to calculate PCA on. Default 500.
#' @param colors Character vector. Optional. Custom colors (named or unnamed). Defaults to ggsci::nrc_npg.
#' @param size Numeric. Size of the scatter points. Default is 5.
#'
#' @return A plotly object.
#' @export
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst
#' @importFrom SummarizedExperiment assay
#' @importFrom matrixStats rowVars
#' @importFrom plotly plot_ly add_trace layout %>%
#' @importFrom rgl ellipse3d
#' @importFrom Morpho quad2trimesh
#' @importFrom stats prcomp cov
#' @importFrom paletteer paletteer_d
plot_pca_3d <- function(bs_obj,
                        color_by,
                        ntop = 500,
                        colors = NULL,
                        size = 5) {

  # 1. Validation
  if (!inherits(bs_obj, "bulkseq")) stop("Input must be a 'bulkseq' object.")

  counts <- bs_obj$counts
  metadata <- bs_obj$metadata

  if (!color_by %in% names(metadata)) {
    stop(paste("Error: Column", color_by, "not found in metadata."))
  }

  # 2. VST Transform
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData = metadata,
                                        design = ~1)
  dds_vst <- DESeq2::vst(dds, blind = TRUE)

  # 3. PCA Calculation
  rv <- matrixStats::rowVars(SummarizedExperiment::assay(dds_vst))
  top_genes <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]

  pca <- stats::prcomp(t(SummarizedExperiment::assay(dds_vst)[top_genes,]))
  percentVar <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)

  # 4. Prepare Data Frame
  pc_df <- data.frame(metadata[[color_by]], pca$x[, 1:3])
  colnames(pc_df)[1] <- "Group"
  pc_df$Group <- factor(pc_df$Group)

  # 5. Handle Colors (Robust)
  levels_list <- levels(pc_df$Group)
  n_groups <- length(levels_list)

  # Default colors if none provided
  if (is.null(colors)) {
    cols_gen <- paletteer::paletteer_d("ggsci::nrc_npg")
    if (n_groups > length(cols_gen)) {
      cols_gen <- rep(cols_gen, length.out = n_groups)
    }
    final_colors <- as.vector(cols_gen)[1:n_groups]
    names(final_colors) <- levels_list
  } else {
    # If user provided named vector, match it. If unnamed, assign order.
    if (!is.null(names(colors))) {
      final_colors <- colors[levels_list]
    } else {
      final_colors <- colors[1:n_groups]
      names(final_colors) <- levels_list
    }
  }

  # 6. Build Plotly Object (Points)
  fig <- plotly::plot_ly(data = pc_df,
                         x = ~PC1, y = ~PC2, z = ~PC3,
                         color = ~Group,
                         colors = final_colors,
                         type = "scatter3d",
                         mode = "markers",
                         marker = list(size = size, opacity = 1.0)) # Use user-defined size

  # 7. Add 3D Ellipsoids (Mesh)
  # We iterate specifically through the levels found in the data
  for (g in levels_list) {
    group_data <- pc_df[pc_df$Group == g, 2:4]

    # Ellipse requires at least 4 points
    if (nrow(group_data) > 3) {
      # Calculate covariance
      e <- rgl::ellipse3d(stats::cov(group_data), centre = colMeans(group_data))
      el <- Morpho::quad2trimesh(e)

      # Get the color specifically for this group
      this_color <- final_colors[[g]]

      # Add the mesh trace
      fig <- fig %>%
        plotly::add_trace(x = el$vb[1,], y = el$vb[2,], z = el$vb[3,],
                          i = el$it[1,]-1, j = el$it[2,]-1, k = el$it[3,]-1,
                          type = "mesh3d",
                          opacity = 0.15,
                          facecolor = rep(this_color, ncol(el$it)), # Apply specific color
                          name = paste0(g, " (ellipse)"),
                          showlegend = FALSE,
                          inherit = FALSE) # Prevent inheriting scatter props
    }
  }

  # 8. Layout
  fig <- fig %>%
    plotly::layout(scene = list(
      xaxis = list(title = paste0("PC1 (", percentVar[1], "%)")),
      yaxis = list(title = paste0("PC2 (", percentVar[2], "%)")),
      zaxis = list(title = paste0("PC3 (", percentVar[3], "%)")),
      bgcolor = "white"
    ),
    title = paste("3D PCA - Colored by", color_by))

  return(fig)
}

#' 3D PCA Plot (Interactive)
#'
#' Generates an interactive 3D PCA plot using plotly. Includes 3D ellipses for groups.
#'
#' @param bs_obj A \code{bulkseq} object created by \code{create_bulkseqvis_object}.
#' @param color_by Character string. Metadata column for point color.
#' @param ntop Integer. Number of most variable genes to calculate PCA on. Default 500.
#' @param group_colors Character vector. Optional. Custom colors. Defaults to ggsci::nrc_npg.
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
#'
#' @examples
#' \dontrun{
#'   plot_pca_3d(my_obj, color_by = "condition")
#' }
plot_pca_3d <- function(bs_obj,
                        color_by,
                        ntop = 500,
                        group_colors = NULL) {

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
  pc_df$Group <- factor(pc_df$Group) # Ensure factor for levels

  # 5. Handle Colors
  n_groups <- length(levels(pc_df$Group))
  if (is.null(group_colors)) {
    group_colors <- paletteer::paletteer_d("ggsci::nrc_npg")
    if (n_groups > length(group_colors)) {
      group_colors <- rep(group_colors, length.out = n_groups)
    }
  }
  # Ensure group_colors is a plain vector to access by index
  group_colors <- as.vector(group_colors)

  # 6. Calculate 3D Ellipsoids
  ellipsoids <- list()
  levels_list <- levels(pc_df$Group)

  for (g in levels_list) {
    group_data <- pc_df[pc_df$Group == g, 2:4]
    # Ellipse requires at least 4 points for 3D covariance stability
    if (nrow(group_data) > 3) {
      # Use colMeans directly (from base), not stats::colMeans
      e <- rgl::ellipse3d(stats::cov(group_data), centre = colMeans(group_data))
      ellipsoids[[g]] <- Morpho::quad2trimesh(e)
    }
  }

  # 7. Build Plotly Object
  fig <- plotly::plot_ly()

  # Add Scatter Points
  fig <- fig %>%
    plotly::add_trace(data = pc_df,
                      x = ~PC1, y = ~PC2, z = ~PC3,
                      type = "scatter3d", mode = "markers",
                      color = ~Group,
                      colors = group_colors[1:n_groups],
                      marker = list(size = 4, opacity = 0.9))

  # Add Ellipsoids (Mesh)
  # We iterate through levels to ensure colors match index
  for (i in seq_along(levels_list)) {
    g_name <- levels_list[i]
    if (g_name %in% names(ellipsoids)) {
      el <- ellipsoids[[g_name]]
      fig <- fig %>%
        plotly::add_trace(x = el$vb[1,], y = el$vb[2,], z = el$vb[3,],
                          i = el$it[1,]-1, j = el$it[2,]-1, k = el$it[3,]-1,
                          type = "mesh3d", opacity = 0.2,
                          facecolor = rep(group_colors[i], ncol(el$it)),
                          showlegend = FALSE)
    }
  }

  # Layout
  fig <- fig %>%
    plotly::layout(scene = list(
      xaxis = list(title = paste0("PC1 (", percentVar[1], "%)")),
      yaxis = list(title = paste0("PC2 (", percentVar[2], "%)")),
      zaxis = list(title = paste0("PC3 (", percentVar[3], "%)")),
      bgcolor = "white"
    ))

  return(fig)
}

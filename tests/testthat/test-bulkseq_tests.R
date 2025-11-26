# Tell rgl to not use a display window (prevents crashes on headless servers)
options(rgl.useNULL = TRUE)

test_that("create_bulkseqvis_object works correctly", {
  cnts <- matrix(1:10, ncol = 2)
  colnames(cnts) <- c("S1", "S2")
  rownames(cnts) <- paste0("G", 1:5)
  meta <- data.frame(cond = c("A", "B"))
  rownames(meta) <- c("S1", "S2")

  obj <- create_bulkseqvis_object(cnts, meta)
  expect_s3_class(obj, "bulkseq")
  expect_equal(ncol(obj$counts), nrow(obj$metadata))
})

test_that("boxplot runs without error", {
  cnts <- matrix(sample(1:100, 200, replace=TRUE), ncol = 4)
  colnames(cnts) <- c("S1", "S2", "S3", "S4")
  rownames(cnts) <- paste0("G", 1:50)
  meta <- data.frame(cond = c("A", "A", "B", "B"))
  rownames(meta) <- c("S1", "S2", "S3", "S4")
  obj <- create_bulkseqvis_object(cnts, meta)

  p <- overall_count_boxplot(obj, color_by = "cond")
  expect_s3_class(p, "ggplot")
})

test_that("heatmap runs without error", {
  # Setup bigger dummy data for DESeq2/heatmap
  n_genes <- 2000
  n_samples <- 10
  cnts <- matrix(sample(10:1000, n_genes * n_samples, replace=TRUE), ncol = n_samples)
  colnames(cnts) <- paste0("S", 1:n_samples)
  rownames(cnts) <- paste0("Gene", 1:n_genes)
  meta <- data.frame(cond = rep(c("A", "B"), each = 5))
  rownames(meta) <- colnames(cnts)

  obj <- create_bulkseqvis_object(cnts, meta)
  expect_error(sample_distance_heatmap(obj, color_by = "cond"), NA)
})

test_that("PCA 2D plot runs without error", {
  n_genes <- 2000
  n_samples <- 10
  cnts <- matrix(sample(10:1000, n_genes * n_samples, replace=TRUE), ncol = n_samples)
  colnames(cnts) <- paste0("S", 1:n_samples)
  rownames(cnts) <- paste0("Gene", 1:n_genes)
  meta <- data.frame(cond = rep(c("A", "B"), each = 5))
  rownames(meta) <- colnames(cnts)

  obj <- create_bulkseqvis_object(cnts, meta)
  p <- plot_pca_2d(obj, color_by = "cond")
  expect_s3_class(p, "ggplot")
})

test_that("PCA 3D plot runs without error", {
  n_genes <- 2000
  n_samples <- 10
  cnts <- matrix(sample(10:1000, n_genes * n_samples, replace=TRUE), ncol = n_samples)
  colnames(cnts) <- paste0("S", 1:n_samples)
  rownames(cnts) <- paste0("Gene", 1:n_genes)
  meta <- data.frame(cond = rep(c("A", "B"), each = 5))
  rownames(meta) <- colnames(cnts)

  obj <- create_bulkseqvis_object(cnts, meta)
  p <- plot_pca_3d(obj, color_by = "cond")
  expect_s3_class(p, "plotly")
})

test_that("UMAP plot runs without error", {
  # UMAP requires slightly more data to avoid errors about 'n_neighbors'
  n_genes <- 2000
  n_samples <- 20 # Increased samples for UMAP stability
  cnts <- matrix(sample(10:1000, n_genes * n_samples, replace=TRUE), ncol = n_samples)
  colnames(cnts) <- paste0("S", 1:n_samples)
  rownames(cnts) <- paste0("Gene", 1:n_genes)
  meta <- data.frame(cond = rep(c("A", "B"), each = 10))
  rownames(meta) <- colnames(cnts)

  obj <- create_bulkseqvis_object(cnts, meta)
  p <- plot_umap_2d(obj, color_by = "cond")
  expect_s3_class(p, "ggplot")
})

test_that("t-SNE plot runs without error", {
  # t-SNE needs N > perplexity * 3 + 1 generally
  n_genes <- 2000
  n_samples <- 30
  cnts <- matrix(sample(10:1000, n_genes * n_samples, replace=TRUE), ncol = n_samples)
  colnames(cnts) <- paste0("S", 1:n_samples)
  rownames(cnts) <- paste0("Gene", 1:n_genes)
  meta <- data.frame(cond = rep(c("A", "B", "C"), each = 10))
  rownames(meta) <- colnames(cnts)

  obj <- create_bulkseqvis_object(cnts, meta)

  # Note: perplexity automatically adjusts in our function if N is small
  p <- plot_tsne_2d(obj, color_by = "cond")
  expect_s3_class(p, "ggplot")
})

test_that("DE analysis, Volcano, Barplot, and Log2FC Comp run without error", {
  # Generate synthetic data
  n_genes <- 1000
  n_samples <- 9
  cnts <- matrix(as.integer(sample(10:1000, n_genes * n_samples, replace=TRUE)), ncol = n_samples)
  colnames(cnts) <- paste0("S", 1:n_samples)
  rownames(cnts) <- paste0("Gene", 1:n_genes)

  # 3 groups to allow 2 comparisons
  meta <- data.frame(condition = rep(c("Ctrl", "TrtA", "TrtB"), each = 3))
  rownames(meta) <- colnames(cnts)

  obj <- create_bulkseqvis_object(cnts, meta)

  # Run DE for two contrasts
  obj <- DEG(obj, design_col = "condition", compare_levels = c("TrtA", "Ctrl"), biomart_dataset = NULL)
  obj <- DEG(obj, design_col = "condition", compare_levels = c("TrtB", "Ctrl"), biomart_dataset = NULL)

  # Check results
  expect_true("TrtA_vs_Ctrl" %in% names(obj$DE_results))
  expect_true("TrtB_vs_Ctrl" %in% names(obj$DE_results))

  # Run Volcano
  p_vol <- plot_volcano(obj, comparison_id = "TrtA_vs_Ctrl")
  expect_s3_class(p_vol, "ggplot")

  # Run Barplot
  p_bar <- plot_deg_bar(obj)
  expect_s3_class(p_bar, "ggplot")

  # Run Log2FC Comparison
  p_fc <- plot_fcvsfc(obj, name1 = "TrtA_vs_Ctrl", name2 = "TrtB_vs_Ctrl")
  expect_s3_class(p_fc, "gtable")
})

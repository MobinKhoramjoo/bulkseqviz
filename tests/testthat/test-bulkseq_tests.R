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

test_that("DE analysis runs without error", {
  # Generate synthetic data
  # We need replicates to run DESeq2 successfully
  n_genes <- 1000
  n_samples <- 6
  # Ensure counts are integers for DESeq2
  cnts <- matrix(as.integer(sample(10:1000, n_genes * n_samples, replace=TRUE)), ncol = n_samples)
  colnames(cnts) <- paste0("S", 1:n_samples)
  rownames(cnts) <- paste0("Gene", 1:n_genes)

  meta <- data.frame(condition = rep(c("Control", "Treat"), each = 3))
  rownames(meta) <- colnames(cnts)

  obj <- create_bulkseqvis_object(cnts, meta)

  # Run DE (skip BioMart for testing to avoid internet deps)
  obj <- DEG(obj,
             design_col = "condition",
             compare_levels = c("Treat", "Control"),
             biomart_dataset = NULL) # Skip annotation

  # Check if results were stored
  expect_true("Treat_vs_Control" %in% names(obj$DE_results))

  res <- obj$DE_results[["Treat_vs_Control"]]
  expect_true("log2FoldChange" %in% colnames(res))
  expect_true("padj" %in% colnames(res))
})

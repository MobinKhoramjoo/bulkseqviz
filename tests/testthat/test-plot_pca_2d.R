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
  # DESeq2 VST requires > 1000 genes to estimate dispersion trends securely
  n_genes <- 2000
  n_samples <- 10

  cnts <- matrix(sample(10:1000, n_genes * n_samples, replace=TRUE), ncol = n_samples)
  colnames(cnts) <- paste0("S", 1:n_samples)
  rownames(cnts) <- paste0("Gene", 1:n_genes)

  meta <- data.frame(cond = rep(c("A", "B"), each = 5))
  rownames(meta) <- colnames(cnts)

  obj <- create_bulkseqvis_object(cnts, meta)

  # We expect pheatmap to return a list (it actually plots to device)
  expect_error(sample_distance_heatmap(obj, color_by = "cond"), NA)
})

test_that("PCA plot runs without error", {
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

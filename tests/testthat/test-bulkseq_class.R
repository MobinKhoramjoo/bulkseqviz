test_that("create_bulkseq works with correct input", {
  # Setup dummy data
  cnts <- matrix(1:10, ncol = 2)
  colnames(cnts) <- c("S1", "S2")
  rownames(cnts) <- paste0("G", 1:5)

  meta <- data.frame(cond = c("A", "B"))
  rownames(meta) <- c("S1", "S2")

  # Run function
  obj <- create_bulkseqvis_object(cnts, meta)

  # Check expectations
  expect_s3_class(obj, "bulkseq")
  expect_equal(ncol(obj$counts), nrow(obj$metadata))
  expect_equal(colnames(obj$counts), rownames(obj$metadata))
})

test_that("create_bulkseqvis_object fails when names do not match", {
  cnts <- matrix(1:10, ncol = 2)
  colnames(cnts) <- c("S1", "S2")

  meta <- data.frame(cond = c("A", "B"))
  rownames(meta) <- c("X1", "X2") # Totally different names

  expect_error(create_bulkseqvis_object(cnts, meta), "No common sample names")
})

test_that("create_bulkseqvis_object subsets correctly when there is partial overlap", {
  cnts <- matrix(1:15, ncol = 3)
  colnames(cnts) <- c("S1", "S2", "S3")

  meta <- data.frame(cond = c("A", "B"))
  rownames(meta) <- c("S1", "S2") # S3 is missing here

  expect_warning(obj <- create_bulkseqvis_object(cnts, meta), "Mismatch detected")

  # Should result in only 2 samples
  expect_equal(ncol(obj$counts), 2)
  expect_equal(colnames(obj$counts), c("S1", "S2"))
})

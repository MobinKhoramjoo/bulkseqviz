test_that("boxplot runs without error", {
  # Setup dummy data
  cnts <- matrix(sample(1:100, 200, replace=TRUE), ncol = 4)
  colnames(cnts) <- c("S1", "S2", "S3", "S4")
  rownames(cnts) <- paste0("G", 1:50)

  meta <- data.frame(cond = c("A", "A", "B", "B"))
  rownames(meta) <- c("S1", "S2", "S3", "S4")

  obj <- create_bulkseqvis_object(cnts, meta)

  # Test 1: Basic run with defaults
  p <- overall_count_boxplot(obj, color_by = "cond")
  expect_s3_class(p, "ggplot")

  # Test 2: Custom order and colors
  p2 <- overall_count_boxplot(obj, color_by = "cond",
                              group_order = c("B", "A"),
                              group_colors = c("blue", "red"))
  expect_s3_class(p2, "ggplot")
})

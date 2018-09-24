testthat::context("Testing plot_heatmap function")

testthat::test_that("Test if plot_heatmap runs", {
  set.seed(42)

  imported_tx <- dplyr::data_frame(
    gene = paste0("gene", stringr::str_pad(1:100, 3, pad = "0")),
    `sample 1` = stats::rnbinom(n = 100, size = 1, prob = 0.05),
    `sample X` = stats::rnbinom(n = 100, size = 5, prob = 0.05),
    `_31` = stats::rnbinom(n = 100, size = 1, prob = 0.05),
    `X12` = stats::rnbinom(n = 100, size = 10, prob = 0.05),
    `SAMPLE.FIVE` = stats::rnbinom(n = 100, size = 10, prob = 0.05)
  )
  sample_df <- dplyr::data_frame(
    sample = names(imported_tx)[2:length(imported_tx)],
    condition = c("A", "2", "|", "A", "A"),
    color = c("red", "green", "green", "red", "green")
  )


  res <- plot_heatmap(imported_tx, num = 30, show_rownames = TRUE)

  testthat::expect_equal(class(res), "pheatmap")

  testthat::expect_equal(is.na(res$kmeans), TRUE)

  testthat::expect_equal(res$tree_row$merge, structure(c(
    -16L, -5L, -21L, -9L, -11L, 3L, -2L, -7L, 4L, -8L,
    1L, -3L, 2L, 6L, 7L, -14L, 15L, 13L, 9L, 18L, 17L, 14L, 11L,
    10L, 23L, 19L, 25L, 26L, 27L, -23L, -19L, -26L, -15L, -20L, -27L,
    -24L, -25L, -12L, -18L, 5L, -30L, -6L, 8L, -4L, -28L, -29L, -17L,
    12L, -13L, -1L, -22L, 16L, -10L, 22L, 21L, 20L, 24L, 28L
  ), .Dim = c(
    29L,
    2L
  )))

  testthat::expect_error(plot_heatmap(imported_tx, num = 0))

  res2 <- plot_heatmap(head(imported_tx, 10))

  testthat::expect_equal(res2$tree_col$labels, c("sample 1", "sample X", "_31", "X12", "SAMPLE.FIVE"))

  res3 <- plot_heatmap(imported_tx,
    sample_table = sample_df, color_by = c("color", "condition"),
    num = -30, show_rownames = FALSE, scale = "column"
  )

  testthat::expect_equal(res3$tree_col$order, c(1L, 3L, 2L, 4L, 5L))

  res4 <- plot_heatmap(imported_tx,
    sample_table = sample_df, color_by = c("color", "condition"),
    num = -30, show_rownames = FALSE, scale = "none"
  )
})

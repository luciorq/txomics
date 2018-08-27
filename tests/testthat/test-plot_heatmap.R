testthat::context("Testing plot_heatmap function")

testthat::test_that("Test if plot_heatmap runs", {

  set.seed(42)

  imported_transcripts <- dplyr::data_frame(
    gene = paste0("gene",stringr::str_pad(1:100 ,3,pad = "0")),
    `sample 1` = stats::rnbinom(n = 100,size = 1,prob = 0.05),
    `sample X` = stats::rnbinom(n = 100,size = 5,prob = 0.05),
    `_31` = stats::rnbinom(n = 100,size = 1,prob = 0.05),
    `X12` = stats::rnbinom(n = 100,size = 10,prob = 0.05)
  )

  sample_df <- dplyr::data_frame(
    sample = names(imported_transcripts)[2:length(imported_transcripts)],
    condition = c("A","2", "|", "A"),
    color = c("red", "green", "green","red")
  )
  res <- plot_heatmap(imported_transcripts, num = 30,show_rownames = TRUE)

  testthat::expect_equal(class(res), "pheatmap")

  testthat::expect_equal(is.na(res$kmeans), TRUE)

  testthat::expect_equal(res$tree_row$merge,structure(c(-15L, -21L, -6L, -1L, -12L, 2L, 4L, -8L, -14L, -7L,
                                            3L, 1L, -3L, -22L, 11L, -10L, 6L, 12L, 13L, -11L, 17L, 10L, 8L,
                                            22L, 18L, 25L, 21L, 26L, 28L, -23L, -25L, -24L, -13L, -20L, -9L,
                                            -5L, -17L, -16L, -28L, -18L, 5L, -26L, -27L, -19L, -29L, 7L,
                                            15L, 14L, -30L, -4L, 16L, 9L, -2L, 19L, 20L, 24L, 23L, 27L), .Dim = c(29L,
                                                                                                                  2L)))

  testthat::expect_error(plot_heatmap(imported_transcripts, num = 0))

  res2 <- plot_heatmap(head(imported_transcripts, 10))

  testthat::expect_equal(res2$tree_col$labels, c("sample 1", "sample X", "_31", "X12"))

  res3 <- plot_heatmap(imported_transcripts, sample_table = sample_df, color_by = c("color","condition"),
                       num = -30,show_rownames = FALSE, scale = "column")

  testthat::expect_equal(res3$tree_col$order, c(1,3,4,2))

})



###
# library(txomics)
# tx <- imported_transcripts
# expression_matrix <- tx
#
# expression_matrix %>% as.matrix()
#
# all(!is.na(as.numeric(as.matrix(expression_matrix))))
#
# all(!is.na(as.numeric(as.matrix(expression_matrix[,2:length(expression_matrix)]))))
#
# temp_mat <- as.matrix(expression_matrix[,2:length(expression_matrix)])
# rownames(temp_mat) <- expression_matrix %>% dplyr::pull(1)
#
# if(!is.matrix(expression_matrix)){
#   #expression_matrix %>% class()
#   if(!all(!is.na(suppressWarnings(as.numeric(as.matrix(expression_matrix)))))){
#     temp_mat <- as.matrix(expression_matrix[,2:length(expression_matrix)])
#     rownames(temp_mat) <- expression_matrix %>% dplyr::pull(1)
#     expression_matrix <- temp_mat
#   }
# }
#

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
  res <- plot_heatmap(imported_transcripts, num = 30,show_rownames = TRUE)

  testthat::expect_equal(class(res), "pheatmap")

  testthat::expect_equal(is.na(res$kmeans), TRUE)

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

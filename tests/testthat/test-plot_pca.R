testthat::context("Testing plot_pca function")

testthat::test_that("Test if plot_pca runs", {
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

  plot_pca(imported_tx,
    sample_table = sample_df,
    color_by = "color", num = 30
  )

  plot_pca(imported_tx,
    sample_table = sample_df,
    color_by = "sample", num = -10
  )

  testthat::expect_error(plot_pca(imported_tx, num = 0))

  res2 <- plot_pca(head(imported_tx, 5))


  testthat::expect_equal(res2$data$value, c(
    -0.141589618299063, -0.180204758946487, 0.402524792746824,
    -0.933770590153989, 0.853040174652716, -0.262972431433991, 0.752433817191767,
    -0.30550222096098, -0.211692682937508, 0.0277335181407126
  ))

  res3 <- plot_pca(imported_tx)


  testthat::expect_equal(res3$data$PC1, c(
    -9.08978758902333, -1.36567534552286, -8.65240409661001, 10.1086525613534,
    8.99921446980285, -9.08978758902333, -1.36567534552286, -8.65240409661001,
    10.1086525613534, 8.99921446980285
  ))
  # imported_tx_df <- data.frame(
  #  gene = paste0("gene",stringr::str_pad(1:100 ,3,pad = "0")),
  #  `sample 1` = stats::rnbinom(n = 100,size = 1,prob = 0.05),
  #  `sample X` = stats::rnbinom(n = 100,size = 5,prob = 0.05),
  #  `_31` = stats::rnbinom(n = 100,size = 1,prob = 0.05),
  #  `X12` = stats::rnbinom(n = 100,size = 10,prob = 0.05)
  # )
})

# tx <- imported_transcripts; sample_table <- samples_metadata
# color_by = c("condition","population"); num = 500
# breaks_x <- base::pretty(df$log2FoldChange)
# if( isTRUE(max(-log10(df$pvalue)) <= 4) ) {
#  range_y <- c(0,4)
# } else{
#  range_y <- -log10(df$pvalue)
# }
# breaks_y <- base::pretty(range_y)
# ggplot2::scale_x_continuous( breaks = breaks_x, limits = base::range(breaks_x*1.1)) +
# ggplot2::scale_y_continuous( breaks = breaks_y, limits = base::range(breaks_y)) +

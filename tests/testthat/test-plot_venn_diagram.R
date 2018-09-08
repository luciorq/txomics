testthat::context("Testing plot_venn_diagram function")
#
testthat::test_that("Test if plot_venn_diagram runs",{
#
  set.seed(42)
#
  #res_list <- list(de_res_rio_2013,de_res_rio_2015,de_res_bot_2014,de_res_neo_2014)
  #filename<-"temp_test/deg_venn_diagram.pdf";width=9;height=NULL
  temp_de_res <- dplyr::data_frame(
    gene = paste0("gene",stringr::str_pad(1:400 ,6,pad = "0")),
    pvalue = stats::runif(400,0,0.5),
    padj = stats::runif(400,0,1)
  )
  temp_res_1 <- temp_de_res %>%
    dplyr::sample_n(200)
  temp_res_2 <- temp_de_res %>%
    dplyr::sample_n(200)
  temp_res_3 <- temp_de_res %>%
    dplyr::sample_n(200)
  temp_res_4 <- temp_de_res %>%
    dplyr::sample_n(200)

  fs::dir_create("temp_test/plots/")

  testthat::expect_error(plot_venn_diagram())
  testthat::expect_error(plot_venn_diagram(temp_res_1,filename = "temp_test/plots/temp_test_1.pdf"))
  testthat::expect_error(plot_venn_diagram(temp_res_1,filename = "temp_test/plots/temp_test_1.npg"))


#  plot_venn_diagram(de_res_neo_2014, de_res_rio_2013,
#                    de_res_rio_2015, de_res_bot_2014,
#                    names = c( "Neopolis 2014", "Rio 2013", "Rio 2015", "Botucatu 2014"),
#                    filename = "temp_test/plots/venn_diagram_deg.pdf",
#                    width = 9)


  res <- plot_venn_diagram(temp_res_1, temp_res_2,
                    temp_res_3, temp_res_4,
                    names = c("x","!zsad","     ","_123"),
                    filename = "temp_test/plots/venn_diagram_deg2.pdf",
                    width = 9)
  testthat::expect_equal(res,NULL)

  res_list <- list(temp_res_4,temp_res_3,temp_res_2,temp_res_1)
  #res <- plot_venn_diagram(res_list,filename = "temp_test/x.pdf")


  #testthat::expect_equal()


  #plot_venn_diagram( list(
  #  sample_A = temp_res_1,
  #  `Sample B`= temp_res_2,
  #  green = temp_res_3,
  #  `_Sample.?` = temp_res_4
  #), filename = "temp_test/plots/named_list_venn.pdf", fdr = TRUE)

  fs::dir_delete("temp_test/")

})
#

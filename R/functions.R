#' \code{txomics} package
#'
#' Functional Analysis and Visualization of Transcriptomic Data
#'
#' See the README on
#' \href{https://github.com/luciorq/txomics#readme}{GitHub}
#'
#' @docType package
#' @name txomics
#' @importFrom dplyr %>%
#'
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
## reference: https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))


.onAttach <- function(libname, pkgname) {
  ##suppressPackageStartupMessages()
  packageStartupMessage("Attaching txomics package")
}

#' retrieve metadata from NCBI SRA
#'
#' retrieve metadata from NCBI SRA
#'
retrive_metadata_sra <- function() {
  ## This function retrieve metadata from sequencing run
  ## Using NCBI SRA
  ## INCOMPLETE
}

#' Import Expression Data
#'
#' Function used to import transcripts count and abundance data
#'
#' @param dir path to files
#'
#' @param source method used
#'
#' @param names source of gene names identifier
#'
#' @return A \code{txi} object containing
#'
#' @export
#'
import_tx <- function(dir, source = "salmon", names = "vectorbase" ){
  ## Function used to import transcripts abundance data
  #library(tidyverse)
  #library(tximport)
  #message("using tximport")
  files <- file.path(dir, list.files(dir), "quant.sf")
  if (!all(file.exists(files))) {
    stop("Not all files exists")
    #files <- files[file.exists(files)]
  }
  names(files) <- list.files(dir) %>%
    stringr::str_remove("_quant$")
  tx_vector <- readr::read_delim(files[1],
                                 "\t", escape_double = FALSE, trim_ws = TRUE) %>%
    .$Name
  if(names == "vectorbase"){
    gene_vector <- tx_vector %>% stringr::str_remove("-R.*")
    tx2gene <- data.frame("TXNAME" = tx_vector, "GENEID" = gene_vector)
  }else{
    tx2gene <- data.frame("TXNAME" = tx_vector, "GENEID" = tx_vector)
  }
  txi <- tximport::tximport(files, type = source, tx2gene = tx2gene)
  txi
}

#' Infer Library Construction Type
#'
#' Retrieve optimal library preparation
#' protocol used for Salmon method of
#' mapping.
#'
#' @param dir path to Salmon results files
#'
#' @return A \code{tibble} object
#'
#' @export
#'
salmon_libtype <- function(dir) {
  ## Function to retrieve Salmon library type from files in Path
  salmon_libtype_ <- function(lib, dir) {
    ## Function to retrieve Salmon library from file
    lib_type <- paste0(dir,"/",lib,"/lib_format_counts.json") %>%
      jsonlite::read_json() %>%
      .$expected_format
    dplyr::data_frame(lib, lib_type)
  }
  files <- list.files(dir)
  purrr::map_df(files, salmon_libtype_, dir)
}

#' Retrieve Mapping Rate
#'
#' Retrieve total number of fragments,
#' total mapped reads and
#' mapping rate from a Salmon mapped
#' sample.
#'
#' @param dir path to Salmon results files
#'
#' @return A \code{tibble} object
#'
#' @export
#'
retrieve_mapping_rate <- function( dir ) {
  libs <- list.files(dir)
  purrr::map_df(libs, ~{
    salmon_log <- readr::read_lines(paste0(dir,"/",.x,"/logs/salmon_quant.log"))
    mapping_pct <- salmon_log %>%
      stringr::str_extract("info.*Mapping rate = .*%") %>%
      .[!is.na(.)] %>%
      stringr::str_extract("(\\d)+.(\\d)+%")
    mapped_reads <- salmon_log %>%
      stringr::str_extract("info.*Counted.*total reads") %>%
      .[!is.na(.)] %>%
      stringr::str_extract("(\\d)+")
    total_reads <- salmon_log %>%
      stringr::str_extract("Observed.*total fragments") %>%
      .[!is.na(.)] %>%
      stringr::str_extract("(\\d)+")
    dplyr::data_frame(sample = .x,mapped_reads,total_reads,mapping_pct) %>%
      dplyr::distinct()
  })
}

#' Filter Imported Transcripts
#'
#' Filter salmon quantification files
#' in a \code{txi} format
#' based on metadata column and levels
#'
#' @param txi txi object
#'
#' @param sample_table metadata data frame
#'
#' @param var_column variable to filter the input
#'
#' @param var_levels levels of var_column to filter the input
#'
#' @return A \code{txi} object
#'
#' @export
#'
#txi <- imported_transcripts
#sample_table <- metadata_df
#var_column <- "population"
#var_levels <- c("Botucatu2014", "Rio2013")
filter_txi <- function(txi, sample_table, var_column = NULL, var_levels = NULL) {
  first_column <- colnames(sample_table)[1]
  if(!is.null(var_levels)){
    extracted_cols <- sample_table %>%
      dplyr::filter((!!as.name(var_column)) %in% var_levels) %>%
      dplyr::select(1) %>%
      dplyr::filter((!!as.name(first_column)) %in% colnames(txi$counts) ) %>%
      unlist()
  } else{
    extracted_cols <- sample_table %>%
      dplyr::select(1) %>%
      unlist()
  }

  txi$counts <- txi$counts[, extracted_cols]
  txi$abundance <- txi$abundance[, extracted_cols]
  txi$length <- txi$length[, extracted_cols]
  txi
}

#' Perform Differential Expression Analysis
#'
#' Run DESeq2 pipeline for differentially expressed genes.
#'
#' @param txi tximport object containing reads counts
#'
#' @param sample_table table containing experimental design
#'                     and sample info
#'
#' @param contrast_var column name from the sample_table
#'
#' @param  numerator factor from contrast_var to be considered as treatment
#'
#' @param denominator factor from contrast_var to be considered as control
#'
#' @param batch_var column from sample_table to be considered for batch effect
#'
#' @param beta_prior default = TRUE, set to FALSE to use apeglm method
#'
#' @param force_rep default = FALSE, experimental
#'
#' @return A \code{tibble} object
#'
#' @importFrom stats as.formula
#'
#' @export
#'
#txi <- tx_rio_2013; sample_table = metadata_df = samples_metadata
#contrast_var <- "condition"; numerator <- "resistant"; denominator <- "susceptible"
de_analysis <- function(txi, sample_table, contrast_var,
                        numerator, denominator,
                        batch_var = NULL,
                        beta_prior = TRUE,
                        force_rep = FALSE){
  lib_names <- colnames( txi$counts )
  var_column <- colnames(sample_table[1])
  sample_table <- sample_table %>%
    dplyr::filter( !!as.name(var_column) %in% lib_names )
  ##library(DESeq2)
  if(is.null(batch_var)){
    design_formula <- as.formula(paste0("~",contrast_var))
  } else{
    design_formula <- as.formula(paste0("~",batch_var," + ",contrast_var))
  }
  dds1 <- DESeq2::DESeqDataSetFromTximport(txi, sample_table, design = design_formula)

  number_of_numerator_samples <- sample_table %>%
    dplyr::filter(!!as.name(contrast_var) == numerator) %>%
    base::nrow()
  number_of_denominator_samples <- sample_table %>%
    dplyr::filter(!!as.name(contrast_var) == denominator) %>%
    base::nrow()
  sample_replicates <- TRUE
  if(number_of_numerator_samples < 2) {
    sample_replicates <- FALSE
  }
  if(number_of_denominator_samples < 2) {
    sample_replicates <- FALSE
  }

  ## START Experimental
  if(isTRUE(force_rep)){
    dds <- DESeq2::DESeq(dds1, betaPrior = TRUE)
    res <- DESeq2::results(dds, alpha = 0.05,
                           contrast = c(contrast_var,numerator,denominator))
    if(!isTRUE("data.frame" %in% class(res))){
      res <- res %>%
        base::as.data.frame() %>%
        dplyr::as_tibble(rownames = "gene")
    }
    return(res)
  }
  ## END Experimental
  if(isTRUE(sample_replicates)){
    if(isTRUE(beta_prior)){
      dds <- DESeq2::DESeq(dds1, betaPrior = TRUE)
      res <- DESeq2::results(dds, alpha = 0.05,
                             contrast = c(contrast_var,numerator,denominator))
    } else{
      dds <- DESeq2::DESeq(dds1, betaPrior = FALSE)
      res <- DESeq2::lfcShrink(dds, parallel = TRUE,
                               coef=paste0(contrast_var,"_",numerator,"_vs_",denominator),
                               type="apeglm")
    }
  } else {
    #txi_row_names <- base::rownames(txi$abundance)
    numerator_name <- sample_table %>%
      dplyr::filter(!!as.name(contrast_var) == numerator)
    numerator_name <- as.character(numerator_name[,1])
    denominator_name <- sample_table %>%
      dplyr::filter(!!as.name(contrast_var) == denominator)
    denominator_name <- as.character(denominator_name[,1])

    res <- txi$abundance %>%
      dplyr::as_tibble(rownames = "gene") %>%
      dplyr::mutate(FC = !!as.name(numerator_name)/!!as.name(denominator_name)) %>%
      dplyr::arrange( -FC ) %>%
      dplyr::mutate(log2FoldChange = log2(FC)) %>%
      dplyr::mutate(log2FoldChange = dplyr::if_else(is.finite(log2FoldChange), log2FoldChange,NA_real_)) %>%
      dplyr::mutate(pvalue = NA_real_) %>%
      dplyr::mutate(padj = NA_real_) %>%
      dplyr::select( -(!!numerator_name)) %>%
      dplyr::select(-(!!denominator_name)) %>%
      dplyr::select( -FC ) %>%
      dplyr::arrange(gene)
  }
  #cat(summary(res))
  if(!isTRUE("data.frame" %in% class(res))){
    res <- res %>%
      base::as.data.frame() %>%
      dplyr::as_tibble(rownames = "gene")
  }
  res
}

#' Perform Gene Set Enrichment Analysis
#'
#' Perform GSEA analysis.
#'
#' @param de_res results from de_analysis function
#'
#' @param gene_sets tidy data frame of the gene sets or
#'                  named list of genes named by gene set
#'
#' @param file_ext File extension where the gene
#'                 sets should be loaded GMT is the default for GSEA
#'                 rds to load rds saved files. DEPRECATED
#'
#' @return A \code{tibble} containing GSEA results
#'
#' @export
#'
#de_res <- de_res_rio2013
#gene_sets <- aaegdata::kegg_gene_sets
gsea_analysis <- function(de_res, gene_sets, file_ext = "gmt") {
  de_res$log2FoldChange[is.na(de_res$log2FoldChange)] <- 0
  sorted_res <- de_res %>%
    dplyr::select(gene, log2FoldChange) %>%
    dplyr::arrange(log2FoldChange)
  ## ordered from strongest down-regulated to strongest upregulated
  ranked_genes <- sorted_res$log2FoldChange
  names(ranked_genes) <- sorted_res$gene
  if(isTRUE("data.frame" %in% class(de_res))) {
    gene_set_col <- base::colnames(gene_sets)[1]
    temp_list <- gene_sets %>%
      dplyr::group_by(!!as.name(gene_set_col)) %>%
      dplyr::summarise(gene = list(gene))
    gene_sets_list <- temp_list$gene
    names(gene_sets_list) <- temp_list %>%
      base::as.data.frame() %>%
      .[,gene_set_col]
    gene_sets <- gene_sets_list
    if(!is.list(gene_sets)){
      if(file_ext == "gmt") {
        gene_sets <- fgsea::gmtPathways(gene_sets)
      } else {
        if( file_ext == "rds") {
          gene_sets <- readr::read_rds(gene_sets)
        }
      }
    }
  }
  gene_sets <- purrr::map(gene_sets, unique)
  fgsea::fgsea(pathways = gene_sets,
               stats = ranked_genes,
               nperm = 100000)
}

#' Perform Leading Edge Analysis
#'
#' Extract Leading Edge genes from GSEA results.
#'
#' @param gsea_res result from gsea_analysis
#'
#' @param set gene set to be extracted
#'
#' @return A \code{vector} of leading edge genes for the set
#'
#' @export
#'
le_analysis <- function(gsea_res, set){
  gsea_res %>%
    tidyr::unnest() %>%
    dplyr::filter(pathway == set) %>%
    .$leadingEdge
}

#' Retrieve Differentially Expressed Genes
#'
#' Extract differentially expressed genes from
#' DESeq2 results based on p-value cuttof.
#'
#' @param de_res de_analysis result
#'
#' @param alpha default = 0.05
#'
#' @param fdr logical, should consider False discovery rate
#'            default = FALSE
#'
#' @return a \code{vector} of gene names
#'
#' @export
#'
retrieve_deg <- function(de_res, alpha = 0.05 ,fdr = FALSE){
  if (!isTRUE("data.frame" %in% class(de_res))) {
    de_res <- de_res %>%
      base::as.data.frame() %>%
      dplyr::as_tibble(rownames = "gene")
  }
  if( isTRUE(fdr)){
    deg_genes <- de_res %>%
      dplyr::filter(padj <= alpha) %>%
      .$gene
  } else {
    deg_genes <- de_res %>%
      dplyr::filter(pvalue <= alpha) %>%
      .$gene
  }
  deg_genes
}

#use MAP family functions instead of FOR loops in R
#
#Use map family of functions from purrr package.
#
#example:
#
#  map_dbl(mtcars, mean)
#  map_dbl(mtcars, median)
#
### multiple functions abstraction
#  funs <- list(mean, median, sd)
#  funs %>%
#    map(~ mtcars %>% map_dbl(.x))

###############################################
##                                            #
##             Plot Functions                 #
##                                            #
###############################################
#' Plot a Heatmap From Expression Data
#'
#' Plot a heatmap from expression data
#'
#' @param tx txi object, matrix or data frame with trancript abundance
#'
#' @param sample_table metadata data frame
#'
#' @param color_by variable from sample table to group samples in the plot
#'
#' @param num number of genes to plot, filtering by variance between samples
#'
#' @export
#'
#tx <- imported_transcripts; sample_table <- samples_metadata
#color_by = c("condition","population"); num = 500
plot_heatmap <- function(tx, sample_table, color_by = NULL, num = NULL) {
  RowVar <- function(x, ...) {
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  }
  if ( is.list(tx) ) {
    expression_matrix <- tx$abundance
  } else {
    expression_matrix <- tx
  }
  if ( is.null(num) ) {
    num <- nrow(expression_matrix)
  }
  mat <- expression_matrix[rowSums(expression_matrix) > 1,]
  filtered_mat <- mat %>%
    dplyr::as_data_frame(rownames = "gene") %>%
    dplyr::mutate( row_var = RowVar(.[2:length(.)])) %>%
    dplyr::arrange(-row_var) %>%
    utils::head(num) %>%
    dplyr::select( -row_var)

  mat <- filtered_mat %>%
    dplyr::select(-gene) %>%
    as.matrix()
  rownames(mat) <- filtered_mat$gene
  if ( is.null(color_by) ) {
    pheatmap::pheatmap( mat, scale = "row",
                        show_rownames = FALSE,
                        border_color = NA,
                        treeheight_row = 20
    )
  } else {
    annotation_df <- sample_table %>%
      dplyr::select( sample,!!color_by  )
    row_names <- annotation_df$sample
    annotation_df <- annotation_df %>%
      dplyr::select( -sample ) %>%
      as.data.frame()
    rownames(annotation_df) <-  row_names
    pheatmap::pheatmap( mat, scale = "row",
                        annotation_col = annotation_df,
                        show_rownames = FALSE,
                        border_color = NA,
                        treeheight_row = 20
    )
  }
}
#' Plot Two Dimensional PCA
#'
#' Plot the results of a two dimensional principal component
#' analysis.
#'
#' @param tx txi object, matrix or data frame with trancript abundance
#'
#' @param sample_table metadata data frame
#'
#' @param color_by variable from sample table to group samples in the plot
#'
#' @param num number of genes to plot, filtering by variance between samples
#'
#' @importFrom stats prcomp
#'
#' @importFrom utils head
#'
#' @export
#'
#tx <- imported_transcripts; sample_table <- samples_metadata
#color_by = c("condition","population"); num = 500
plot_pca <- function (tx, sample_table, color_by = NULL, num = NULL) {

  RowVar <- function(x, ...) {
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  }
  if(is.list(tx)) {
    expression_matrix <- tx$abundance
  } else{
    expression_matrix <- tx
  }
  if( is.null(num) ) {
    num <- nrow(expression_matrix)
  }
  mat <- expression_matrix[rowSums(expression_matrix) > 1,]
  filtered_mat <- mat %>%
    dplyr::as_data_frame(rownames = "gene") %>%
    dplyr::mutate( row_var = RowVar(.[2:length(.)])) %>%
    dplyr::arrange(-row_var) %>%
    utils::head(num) %>%
    dplyr::select( -row_var)

  mat <- filtered_mat %>%
    dplyr::select(-gene) %>%
    as.matrix()
  rownames(mat) <- filtered_mat$gene

  pca <- mat %>%
    t() %>%
    stats::prcomp()

  percentVar <- pca$sdev^2/sum(pca$sdev^2)

  pc_1 <- pca$x[, 1]
  pc_1_df <- dplyr::data_frame(
    sample = names(pc_1),
    PC1 = pc_1
  )
  pc_2 <- pca$x[, 2]
  pc_2_df <- data_frame(
    sample = names(pc_2),
    PC2 = pc_2
  )
  if (is.null(color_by)) {
    ## plot without colors
    pc_1_df %>%
      dplyr::left_join(pc_2_df, by = "sample") %>%
      dplyr::left_join(group_df, by = "sample") %>%
      ggplot2::ggplot( ggplot2::aes(x = PC1, y = PC2, color = sample)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ggplot2::ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      ggplot2::coord_fixed() +
      ggpubr::theme_pubr() %>%
      return()

  }
  if (!all(color_by %in% names(sample_table))) {
    stop("the values of 'color_by' should specify columns of sample_table metadata")
  }

  group_df <- sample_table %>%
    dplyr::select( sample, !!color_by )

  df <- pc_1_df %>%
    dplyr::left_join(pc_2_df, by = "sample") %>%
    dplyr::left_join(group_df, by = "sample") %>%
    tidyr::unite( group, !!color_by )

  df %>%
    ggplot2::ggplot( ggplot2::aes(x = PC1, y = PC2, color = group)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
    ggplot2::ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
    ggplot2::coord_fixed() +
    ggpubr::theme_pubr()
}#' Volcano Plot
#'
#' A volcano plot from DE genes analysis results.
#'
#' @param de_res de_analysis result
#'
#' @param alpha default = 0.05
#'
#' @param lfc_threshold absolute value of log2FoldChange to limit the significance
#'
#' @param fdr logical, should consider False discovery rate
#'            default = FALSE
#'
#' @return A \code{ggplot} object
#'
#' @export
#'
#de_res <- de_res_rio_2013; alpha <- 0.05;fdr <- TRUE; lfc_threshold <- 0.9
plot_volcano <- function(de_res, alpha = 0.05, lfc_threshold = NULL, fdr = FALSE ) {
  ##scatter_plot <- function(data, x, y) {
  ##  x <- enquo(x)
  ##  y <- enquo(y)
  ##
  ##  ggplot(data) + geom_point(aes(!!x, !!y))
  ##}
  ##scatter_plot(mtcars, disp, drat)
  df <- de_res %>%
    dplyr::select(gene, log2FoldChange, pvalue, padj)
  if (isTRUE(fdr)){
    df <- df %>%
      dplyr::filter(!is.na(padj)) %>%
      dplyr::mutate(sig_alpha = dplyr::if_else( padj < alpha, 1, 0 ) )
  } else {
    df <- df %>%
      dplyr::filter(!is.na(pvalue)) %>%
      dplyr::mutate(sig_alpha = dplyr::if_else( pvalue < alpha, 1, 0 ) )
  }
  if (is.null(lfc_threshold)) {
    if (base::nrow(dplyr::filter(df, sig_alpha == 1 )) == 0  ){
      lfc_cutoff <- 2
    } else{
      lfc_cutoff <- df %>%
        dplyr::filter(sig_alpha == 1) %>%
        .$log2FoldChange %>%
        abs() %>%
        max()
    }
  } else{
    lfc_cutoff <- lfc_threshold
  }
  df <- df %>%
    dplyr::mutate( sig_lfc = dplyr::if_else(abs(log2FoldChange) > lfc_cutoff, 2, 0 ) )
  df <- df %>%
    dplyr::mutate(sig = as.character(sig_alpha + sig_lfc)) %>%
    dplyr::select(-c(sig_alpha, sig_lfc)) %>%
    dplyr::filter(!is.na(log2FoldChange))
  df$sig <- factor(df$sig, levels = c("0", "1", "2", "3") )

  limits_x <- max(abs(df$log2FoldChange)) * 1.2
  #if( limits_x < 2 ) {
  #  limits_x <- 2
  #}

  p <- df %>%
    ggplot2::ggplot( ggplot2::aes(log2FoldChange, -log10(pvalue)) ) +
    ggplot2::geom_point( ggplot2::aes( color = sig ) ) +
    ggplot2::scale_color_manual(
      values=c(`0` = "black",`1` = "red",`2` = "red",`3` = "orange")
    ) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(df, sig %in% "3" ),
      ggplot2::aes(label = gene)
    ) +
    ggplot2::geom_hline(yintercept = -log10(alpha), alpha = 0.5 ) +
    ggplot2::geom_vline(xintercept = lfc_cutoff, alpha = 0.5 ) +
    ggplot2::geom_vline(xintercept = -lfc_cutoff, alpha = 0.5 ) +
    ggplot2::expand_limits(x = c(-limits_x,limits_x)) +
    ggplot2::expand_limits(y = c(0,5)) +
    #ggplot2::coord_fixed() +
    ggpubr::theme_pubr()
  p
}
#' Plot GSEA Enrichment Curve and Table
#'
#' Plots enrichment table and individual enrichment graph for
#' GSEA analysis results
#'
#' @param gsea_res result from gsea_analysis
#'
#' @param de_res results from de_analysis function
#'
#' @param gene_sets named list of genes named by gene set
#'
#' @param top_sets Number of sets to be plotted in the gsea table,
#'                 the number of sets to be plotted is the double of the value
#'                 of top_sets
#'                 or a string for individual gene set for enrichment plot
#'
#' @param enrichment If set to TRUE, enrichment for individual
#'                   gene set will be plotted.
#'                   The default value is FALSE, for plotting GSEA table.
#'
#' @param title title label to be shown for enchment plot,
#'              should be a character vector of the text to be plotted
#'              or TRUE for using the .default = NULL, no text plotted.
#'
#' @return graphs for gsea table or gsea individual set
#'
#' @export
#'
#gsea_res <- gsea_res_rio_2013; de_res <- de_res_rio_2013; gene_sets <- aaegdata::kegg_gene_sets;
#top_sets <- 10; enrichment <- FALSE; title = TRUE
plot_gsea <- function(gsea_res, de_res, gene_sets,
                      top_sets = 10, enrichment = FALSE,
                      title = NULL) {
  de_res$log2FoldChange[is.na(de_res$log2FoldChange)] <- 0
  sorted_res <- de_res %>%
    dplyr::select(gene, log2FoldChange) %>%
    dplyr::arrange(log2FoldChange)
  ## ordered from strongest down-regulated to strongest upregulated
  ranked_genes <- sorted_res$log2FoldChange
  names(ranked_genes) <- sorted_res$gene
  if(isTRUE("data.frame" %in% class(de_res))) {
    gene_set_col <- base::colnames(gene_sets)[1]
    temp_list <- gene_sets %>%
      dplyr::group_by(!!as.name(gene_set_col)) %>%
      dplyr::summarise(gene = list(gene))
    gene_sets_list <- temp_list$gene
    names(gene_sets_list) <- temp_list %>%
      base::as.data.frame() %>%
      .[,gene_set_col]
    gene_sets <- gene_sets_list
    if(!is.list(gene_sets)){
      if(file_ext == "gmt") {
        gene_sets <- fgsea::gmtPathways(gene_sets)
      } else {
        if( file_ext == "rds") {
          gene_sets <- readr::read_rds(gene_sets)
        }
      }
    }
  }
  gene_sets <- purrr::map(gene_sets, unique)

  if(isTRUE(enrichment)){
    if(isTRUE(title)){
      title <- top_sets
    }
    fgsea::plotEnrichment(gene_sets[[top_sets]],
                          ranked_genes) +
      ggplot2::labs(title = title) +
      ggpubr::theme_pubr()
    #labs(title="Programmed Cell Death")
  } else {
    topPathwaysUp <- gsea_res[ES > 0, ][utils::head(order(pval), n=top_sets), pathway]
    topPathwaysDown <- gsea_res[ES < 0, ][utils::head(order(pval), n=top_sets), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    fgsea::plotGseaTable(gene_sets[topPathways],
                         ranked_genes,
                         gsea_res,
                         colwidths = c(9, 3, 0.8, 1.2, 1.2)
    )
  }
}

#' Multiple Plot
#'
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#' - cols:   Number of columns in layout
#' - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' adapted from ggplot2 cookbook.
#'
#' @param ... ggplot objects
#'
#' @param plotlist list of ggplot objects
#'
#' @param cols number of columns in layout
#'
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored
#'
#' @export
#'
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  #library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
  }
}

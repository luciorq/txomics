#' @keywords internal
#' @importFrom dplyr %>%
"_PACKAGE"

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
retrieve_metadata_sra <- function() {
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
#' @examples
#' \dontrun{
#' # Example with path to folder input
#' imported_transcripts <- import_tx("data/salmon/res_quants/")
#' }
#'
#' @export
#'
import_tx <- function(dir, source = "salmon", names = "vectorbase" ){
  ## Function used to import transcripts abundance data
  files <- file.path(dir, list.files(dir), "quant.sf")
  if (!all(file.exists(files))) {
    stop("Not all files exists")
    #files <- files[file.exists(files)]
  }
  names(files) <- list.files(dir) %>%
    stringr::str_remove("_quant$")
  tx_vector <- readr::read_delim(
      files[1],"\t", escape_double = FALSE, trim_ws = TRUE
    ) %>%
    dplyr::pull(Name)
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
      dplyr::pull(expected_format)
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
  #de_res$log2FoldChange[is.na(de_res$log2FoldChange)] <- 0
  de_res %>%
    dplyr::mutate(
      log2FoldChange = dplyr::if_else(is.na(log2FoldChange), 0, log2FoldChange )
    )
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
  gsea_res <- fgsea::fgsea(pathways = gene_sets,
               stats = ranked_genes,
               nperm = 100000)
  gsea_res <- gsea_res %>%
    dplyr::as_tibble() %>%
    dplyr::rename(pvalue = pval)
  gsea_res
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
    dplyr::pull(leadingEdge)
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
      dplyr::pull(gene)
  } else {
    deg_genes <- de_res %>%
      dplyr::filter(pvalue <= alpha) %>%
      dplyr::pull(gene)
  }
  deg_genes
}

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
#' @param scale logical, default = "row", normalize data by row, column or none
#'
#' @param show_rownames logical, default = FALSE
#'
#' @param ... additional plot parameters, check pheatmap documentation for
#'            details
#'
#' @export
#'
#tx <- imported_transcripts; sample_table <- samples_metadata
#color_by = c("condition","population"); num = 1000
plot_heatmap <- function(tx, sample_table, color_by = NULL,
                         num = NULL, scale = "row",
                         show_rownames = FALSE, ...) {
  calc_var_by_row <- function(x, ...) {
    base::rowSums((x - base::rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
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
  ## Scaling matrix
  if(scale == "row"){
    mat <- t(scale(t(mat),center = TRUE, scale = TRUE))
  }else{
    if(scale == "column"){
      mat <- scale(mat, center = TRUE, scale = TRUE)
    }
  }
  filtered_mat <- mat %>%
    dplyr::as_data_frame(rownames = "gene") %>%
    dplyr::mutate( row_var = calc_var_by_row(.[2:length(.)])) %>%
    dplyr::arrange(-row_var)
  if( isTRUE(num > 0) ) {
    filtered_mat <- filtered_mat %>%
      utils::head(num)
  }
  if( isTRUE(num < 0) ) {
    num <- abs(num)
    filtered_mat <- filtered_mat %>%
      utils::tail(num)
  }
  if( isTRUE(num == 0) ) {
    stop("Number of genes can't be zero.")
  }
  filtered_mat <- filtered_mat %>%
    dplyr::select( -row_var)
  mat <- filtered_mat %>%
    dplyr::select(-gene) %>%
    as.matrix()
  rownames(mat) <- filtered_mat$gene
  if ( is.null(color_by) ) {
    annotation_df <- NA
    annotation_pallete <- NA
  } else {
    annotation_df <- sample_table %>%
      dplyr::select( sample,!!color_by )
    row_names <- annotation_df$sample
    annotation_df <- annotation_df %>%
      dplyr::select( -sample ) %>%
      as.data.frame()
    rownames(annotation_df) <-  row_names
    ## Generate Annotation discrete color pallette
    df_names <- colnames(annotation_df)
    df_lengths <- df_names %>%
      purrr::map_dbl(~length(unique(annotation_df[[.x]])))
    annotation_pallete <- purrr::map(
      seq_along(df_lengths),
      ~{
        discrete_pallete <- df_lengths[.x] %>%
          base::sum() %>%
          viridis::viridis()
        discrete_pallete
      }
    )
    names(annotation_pallete) <- df_names
    #names(mat_colors$group) <- unique(col_groups)
    for( i in colnames(annotation_df)){
      names(annotation_pallete[[i]]) <- unique(annotation_df[,i])
    }
  }

  ## Sorting dendograms
  sort_hclust <- function(...) {
    stats::as.hclust(dendsort::dendsort(stats::as.dendrogram(...)))
  }
  mat_cluster_cols <- stats::hclust(stats::dist(t(mat)))
  mat_cluster_cols <- sort_hclust(mat_cluster_cols)
  mat_cluster_rows <- sort_hclust(stats::hclust(stats::dist(mat)))

  utils::assignInNamespace(
    x = "draw_colnames",
    value = "draw_colnames_45",
    ns = asNamespace("pheatmap")
  )

  pheatmap::pheatmap(
    mat = mat,
    color = viridis::viridis(100, direction = 1),
    #scale = "row",
    annotation_col = annotation_df,
    annotation_colors = annotation_pallete,
    show_rownames = show_rownames,
    border_color = NA,
    treeheight_col = 30,
    treeheight_row = 30,
    drop_levels = TRUE,
    cluster_cols = mat_cluster_cols,
    cluster_rows = mat_cluster_rows,
    ...
    )
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

  calc_var_by_row <- function(x, ...) {
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
  ## Scaling matrix by row
  mat <- t(scale(t(mat),center = TRUE, scale = TRUE))

  filtered_mat <- mat %>%
    dplyr::as_data_frame(rownames = "gene") %>%
    dplyr::mutate( row_var = calc_var_by_row(.[2:length(.)])) %>%
    dplyr::arrange(-row_var)
  if( isTRUE(num > 0) ) {
    filtered_mat <- filtered_mat %>%
      utils::head(num)
  }
  if( isTRUE(num < 0) ) {
    num <- abs(num)
    filtered_mat <- filtered_mat %>%
      utils::tail(num)
  }
  if( isTRUE(num == 0) ) {
    stop("Number of genes can't be zero.")
  }
  filtered_mat <- filtered_mat %>%
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
  pc_2_df <- dplyr::data_frame(
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
    ggplot2::scale_color_viridis_d( direction = 1 ) +
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
        dplyr::pull(log2FoldChange) %>%
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
  color_pallete <- viridis::viridis(3)
  p <- df %>%
    ggplot2::ggplot( ggplot2::aes(log2FoldChange, -log10(pvalue)) ) +
    ggplot2::geom_point( ggplot2::aes( color = sig ) ) +
    ggplot2::scale_color_manual(
      values=c(`0` = color_pallete[1],`1` = color_pallete[2],`2` = color_pallete[2],`3` = color_pallete[3])
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
#gsea_res <- gsea_res_batch; de_res <- de_res_rio_2013; gene_sets <- aaegdata::kegg_gene_sets;
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
    top_paths_up <- gsea_res %>%
      dplyr::filter(ES > 0) %>%
      dplyr::arrange(pvalue) %>%
      utils::head(n = top_sets) %>%
      dplyr::pull(1)
    top_paths_down <- gsea_res %>%
      dplyr::filter(ES < 0) %>%
      dplyr::arrange(pvalue) %>%
      utils::head(n = top_sets) %>%
      dplyr::pull(1)
    gsea_res <- gsea_res %>%
      dplyr::rename(pval = pvalue)
    top_paths <- c(top_paths_up, rev(top_paths_down))
    fgsea::plotGseaTable(gene_sets[top_paths],
                         ranked_genes,
                         gsea_res,
                         colwidths = c(9, 3, 0.8, 1.2, 1.2)
    )
  }
}
#' Plot Venn Diagram
#'
#' Plot Venn diagrams, from de_analysis or gsea_analysis results
#' **Currentry, this version only plots venn diagram with four groups and just
#' save files to file, don't plot directly**
#'
#' @param ... de_analysis or gsea_analysis data frames
#'
#' @param pvalue p-value cut-off for enrichment or differential consideration
#'
#' @param fdr logical, default = FALSE should adjusted pvalue be used
#'               instead of pvalue
#'
#' @param names vector with names
#'
#' @param filename path to were the image should be saved
#'
#' @param width value to be used as the width of the plot, the venn diagram
#'              will be plotted as a square with side of the smaller
#'              value, between width and height
#'
#' @param height value to be used as the height of the plot, the venn diagram
#'               will be plotted as a square of the smaller value, between
#'               width and height. heigth cam be ommited and the plot will be
#'               a square.
#'
#' @export
#'
#res_list <- list(de_res_rio_2013,de_res_rio_2015,de_res_bot_2014,de_res_neo_2014)
#filename = "plots/deg_venn_diagram.pdf";width = 9; height = NULL
plot_venn_diagram <- function(..., pvalue = 0.05, fdr = FALSE, names = NULL,
                              filename, width = 9, height = NULL){
  file_extension <- filename %>% stringr::str_extract("\\.[a-zA-Z]*$")
  if(!isTRUE(identical(file_extension, ".pdf"))){
    ## remove to introduce another extensions
    stop("filename is not PDF")
  }
  res_list <- list(...)
  if(is.null(names)){
    ##
  } else{
    if( isTRUE(identical(length(res_list),length(names)))) {
      names(res_list) <- names
    } else{
      stop("length of names should be the same of the objects passed to '...'")
    }
  }
  file_path <- filename %>% stringr::str_replace(file_extension,"")
  file_name_svg <- paste0(file_path,".svg")

  pvalue_var <- dplyr::quo(pvalue)
  if(isTRUE(fdr)){
    pvalue_var <- dplyr::quo(padj)
  }
  pvalue_cutoff <- pvalue
  input_list <- res_list %>%
    purrr::map(~{
      .x %>%
        dplyr::filter(!!pvalue_var < pvalue_cutoff) %>%
        dplyr::pull(1)
    })
  names(input_list) <- names
  output_height <- height
  output_width <- width
  if(is.null(height)){
    output_height <- width
  }
  venn_width <- min(output_height,output_width)
  VennDiagram::venn.diagram(
    x = input_list,
    filename = file_name_svg, imagetype = "svg",
    height = venn_width, width = venn_width, col = "transparent",
    fill = c("cornflowerblue","green","yellow","darkorchid1"),
    label.col = c("orange", "white", "darkorchid4", "white", "white",
                  "white", "white", "white", "darkblue", "white", "white", "white", "white",
                  "darkgreen", "white"),
    fontfamily = "serif", fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), cat.cex = 1.5,
    cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "serif", alpha = 0.50, cex = 1.5,
    rotation.degree = 270, margin = 0.2
  )
  if(isTRUE(identical(file_extension, ".pdf"))) {
    file_name_pdf <- paste0(file_path,".pdf")
    rsvg::rsvg_pdf(svg = file_name_svg, file = file_name_pdf)
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

  num_plots <- length(plots)


## #' @param filename file path to save. Filetype is decided by the extension in the path.
## #'                 Currently following formats are supported: only pdf.
##                    filename = NA

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(num_plots/cols)),
                     ncol = cols, nrow = ceiling(num_plots/cols))
  }

  if (num_plots == 1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:num_plots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
  }
}

###############################################
##                                            #
##              Aux Functions                 #
##                                            #
###############################################
#' aux for heatmap
#'
#' aux for heatmap
#' @export
draw_colnames_45 <- function (coln, gaps, ...) {
  ## Correct labels orientation, thanks to https://slowkow.com/notes/heatmap-tutorial/
  find_coordinates_temp <- utils::getFromNamespace("find_coordinates", "pheatmap")
  coord <- find_coordinates_temp(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = grid::unit(1, "npc") - grid::unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}

#' @export
retrieve_deg_table <- function(..., pvalue = 0.05, sample_names){
  pvalue_cutoff <- pvalue
  deg_table_population <- list(...) %>%
    purrr::map_df(
      ~{temp_res <- .x
      up_deg <- temp_res %>%
        filter(pvalue < pvalue_cutoff) %>%
        filter(log2FoldChange > 0) %>%
        nrow()
      down_deg <- temp_res %>%
        filter(pvalue < pvalue_cutoff) %>%
        filter(log2FoldChange < 0) %>%
        nrow()
      data_frame(up = up_deg,
                 down = down_deg,
                 total = up_deg + down_deg)
      })
  deg_table_population %>%
    dplyr::mutate(population = sample_names ) %>%
    dplyr::select(population, up, down, total )
}

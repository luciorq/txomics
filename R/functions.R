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
#' @param sample_table table containing experimental design and sequencing info
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
#' @return A \code{tibble} object
#'
#' @importFrom stats as.formula
#'
#' @export
#'
#txi <- txi_rio2013
#txi <- filtered_txi
#sample_table <- metadata_df
#contrast_var <- "condition"
#numerator <- "resistant"
#denominator <- "susceptible"
DE_analysis <- function(txi, sample_table, contrast_var,
                        numerator, denominator,
                        batch_var = NULL, beta_prior = TRUE){
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
      dplyr::mutate(log2FoldChange = log2(FC)) %>%
      dplyr::mutate(log2FoldChange = dplyr::if_else(is.finite(log2FoldChange), log2FoldChange,NA_real_)) %>%
      dplyr::mutate(pvalue = NA_real_) %>%
      dplyr::mutate(padj = NA_real_)
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
#' @param DE_res DESeqResults object containing DE analysis results
#'
#' @param gene_sets named list of genes named by gene set
#'
#' @param file_ext File extension where the gene sets should be loaded GMT is the default for GSEA
#'   rds to load rds saved files.
#'
#' @return A \code{tibble} containing GSEA results
#'
#' @export
#'
#DE_res <- de_res_rio2013
#gene_sets <- aaegdata::kegg_gene_sets
GSEA_analysis <- function(DE_res, gene_sets, file_ext = "gmt") {
  #library(fgsea)
  DE_res$log2FoldChange[is.na(DE_res$log2FoldChange)] <- 0
  sorted_res <- DE_res %>%
    #base::as.data.frame() %>%
    #dplyr::as_tibble(rownames = "gene") %>%
    dplyr::select(gene, log2FoldChange) %>%
    dplyr::arrange(log2FoldChange)
  ## ordered from strongest down-regulated to strongest upregulated
  ranked_genes <- sorted_res$log2FoldChange
  names(ranked_genes) <- sorted_res$gene
  if(isTRUE("data.frame" %in% class(DE_res))) {
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
#' @param GSEA_res fgsea result data frame
#'
#' @param set gene set to be extracted
#'
#' @return A \code{vector} of leading edge genes for the set
#'
#' @export
#'
LE_analysis <- function(GSEA_res, set){
  GSEA_res %>%
    tidyr::unnest() %>%
    dplyr::filter(pathway == set) %>%
    .$leadingEdge
}

#' Retrieve Differentially Expressed Genes
#'
#' Extract differentially expressed genes from
#' DESeq2 results based on p-value cuttof.
#'
#' @param DE_res DESeq2 results object
#'
#' @return a \code{vector} of DEGs
#'
#' @export
#'
retrieve_DEG <- function(DE_res){
  DE_res %>%
    base::as.data.frame() %>%
    dplyr::as_tibble(rownames = "gene") %>%
    dplyr::filter(pvalue <= 0.05) %>%
    .$gene
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

#' Volcano Plot
#'
#' A volcano plot from DE genes analysis results.
#'
#' @param DE_res DESeq2 results object
#'
#' @param lfc_threshold absolute value of log2FoldChange to limit the significance
#'
#' @param FDR use FDR adjusted p-value, logical
#'
#' @return A \code{ggplot} object
#'
#' @export
#'
#DE_res <- de_res_rio2013
#DE_res <- de_res
#DE_res <- res
plot_volcano <- function(DE_res, lfc_threshold = NULL, FDR = FALSE ) {
  DE_res <- DE_res[order(DE_res$pvalue), ]
  results <- DE_res %>%
    dplyr::rename(Gene = 1) %>%
    dplyr::select(Gene, log2FoldChange, pvalue, padj) %>%
    dplyr::mutate(sig = dplyr::if_else(pvalue < 0.05, "p-value < 0.05", "Not significant")) %>%
    dplyr::mutate(sig = dplyr::if_else( (pvalue < 0.05) & (abs(log2FoldChange) <  2), "DEG",sig)) %>%
    dplyr::mutate(sig = dplyr::if_else(is.na(sig),"Not significant",sig))
  #results %>% dplyr::filter(!is.na(sig))
  results$sig <- as.factor(results$sig)
  levels(results$sig) <- c( "Not significant","p-value < 0.05","DEG")
  #levels(results$sig)
  if(is.null(lfc_threshold)){
    lfc_threshold <- results %>%
      dplyr::filter(pvalue > 0.05 ) %>%
      .$log2FoldChange %>%
      max()
  }
  #library(ggrepel)
  results %>%
    #dplyr::filter( !is.na(sig) ) %>%
    ggplot2::ggplot( ggplot2::aes(log2FoldChange, -log10(pvalue)) ) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(results, (pvalue < 0.05)&(abs(log2FoldChange) > lfc_threshold)),
      ggplot2::aes( label = Gene )) +
    ggplot2::geom_point(ggplot2::aes(col=sig)) +
    ggplot2::scale_color_manual(values=c("red","black","orange")) +
    ggplot2::geom_hline(yintercept = -log10(0.05), alpha = 0.5) +
    ggplot2::geom_vline(xintercept = lfc_threshold, alpha = 0.5 ) +
    ggplot2::geom_vline(xintercept = -lfc_threshold, alpha = 0.5 ) +
    #ggtitle( "Volcano Plot") +
    #coord_cartesian(xlim = c( min(DE_res$log2FoldChange) - 2, max(DE_res$log2FoldChange) + 2 ),
    #               ylim = c(min(-log10(DE_res$pvalue)) , max(-log10(DE_res$pvalue)) + 2)) +
    #xlim( c( min(DE_res$log2FoldChange) - 2, max(DE_res$log2FoldChange) + 2 )) +
    #ylim( 0, max(-log10(DE_res$pvalue)) + 2 ) +
    ggpubr::theme_pubr()
}

#' Multiple Plot
#
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

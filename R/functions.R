#' @keywords internal
#' @importFrom dplyr %>%
"_PACKAGE"

# quiets concerns of R CMD check re: the .'s that appear in pipelines
# reference: https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
if (getRversion() >= "2.15.1") utils::globalVariables(c("."))

## To test if a suggested package is installed
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching txomics package")
}

#' retrieve metadata from NCBI SRA
#'
#' If you are using SRA accessions, it's possible to retrieve
#' sample_table or matadata based on the deposited information
#' retrieve metadata from NCBI SRA
#'
retrieve_metadata_sra <- function() {
  # This function retrieve metadata from sequencing run
  # Using NCBI SRA
  # INCOMPLETE
}

#' Import Expression Data
#'
#' Function used to import transcripts count and abundance data
#'
#' @param path path to files, or vector of paths,
#'             if \code{is_dir} is set to FALSE
#'
#' @param source quantification method used.
#'   Currently supported: "salmon", "kallisto".
#'   If method is unknown use "none".
#'
#' @param accession source of gene names identifier
#'              supported formats:"vectorbase", "gene", "custom"
#'              If using "gene", the trancript names will be used as genes;
#'              If using "custom", a conversion table should be
#'              provided in the gene_table parameter;
#'              not supported yet: "ensembl", "ncbi", "custom"
#'
#' @param is_dir if \code{path} is a directory containing
#'               subdirectories with the files,
#'               if \code{is_dir} is set to FALSE, path should be
#'               a vector of paths;
#'               default = TRUE
#'
#' @param names .vector of names to be used to each sample; default = NULL
#'
#' @param gene_table a transcript to gene conversion table,
#'                   where the first column should be transcript name
#'                  and the second column should be gene name,
#'                  should only be used with accession = "custom"
#'
#' @return A \code{txi} object containing
#'
#' @examples
#' \dontrun{
#' # Example with path to folder input
#' imported_transcripts <- import_tx("data/salmon/quants/")
#' }
#'
#' @export
#'
import_tx <- function(path,
                      source = "salmon",
                      accession = "vectorbase",
                      is_dir = TRUE,
                      names = NULL,
                      gene_table = NULL) {
  ## Function used to import transcripts abundance data

  if ( isTRUE(is_dir)) {
    if (isTRUE(source == "salmon")) {
      path_to_files <- c(
        fs::path(fs::dir_ls(path), "quant.sf"),
        fs::path(path, "quant.sf")
      )
      path_to_files <- path_to_files[fs::file_exists(path_to_files) == TRUE]
    } else {
      path_to_files <- fs::dir_ls(path)
      path_to_files <- path_to_files[fs::is_dir(path_to_files)]
      path_to_files <- fs::dir_ls(path_to_files)
    }
  } else {
    path_to_files <- path
  }
  message("Using files:\n", paste(path_to_files, collapse = "\n"), "\n")

  if (!is.null(names)) {
    if (isTRUE(length(names) == length(path_to_files))) {
      names(path_to_files) <- names
    } else {
      stop("Number of names (sample names) should be the same as files to read.")
    }
  } else {
    names(path_to_files) <- path_to_files %>%
      stringr::str_extract("/(.*)/") %>%
      stringr::str_remove("/$") %>%
      stringr::str_remove(".*/")
  }

  tx_accession <- readr::read_delim(
    path_to_files[1], "\t",
    escape_double = FALSE, trim_ws = TRUE
  ) %>%
    dplyr::pull(1)

  if (isTRUE(accession == "vectorbase")) {
    gene_vectorbase <- tx_accession %>% stringr::str_remove("-R.*")
    tx_to_gene <- dplyr::data_frame(
      "TXNAME" = tx_accession,
      "GENEID" = gene_vectorbase
    )
  }
  if (isTRUE(accession == "gene")) {
    tx_to_gene <- dplyr::data_frame(
      "TXNAME" = tx_accession,
      "GENEID" = tx_accession
    )
  }
  if (isTRUE(accession == "custom")) {
    tx_custom <- gene_table %>%
      dplyr::pull(1)
    gene_custom <- gene_table %>%
      dplyr::pull(2)
    tx_to_gene <- dplyr::data_frame(
      "TXNAME" = tx_custom,
      "GENEID" = gene_custom
    )
  }

  if (isTRUE(source == "none")) {
    #col_names <- c("target_id","eff_length","est_counts","tpm")
    tx <- tximport::tximport(
      files = path_to_files,
      type = source,
      tx2gene = tx_to_gene#,
      #txIdCol = col_names[1],
      #abundanceCol = col_names[4],
      #countsCol = col_names[3],
      #lengthCol = col_names[2]
    )
  }else {
    tx <- tximport::tximport(
      files = path_to_files,
      type = source,
      tx2gene = tx_to_gene
    )
  }
  tx
}

#' Infer Library Construction Type
#'
#' Retrieve optimal library preparation
#' protocol used for Salmon method of
#' mapping.
#'
#' @param path path to Salmon results files
#'
#' @return A \code{tibble} object
#'
#' @export
#'
# path <- "data/salmon/res_quants/"
salmon_libtype <- function(path) {
  files <- list.files(path)
  files %>%
    purrr::map_df(
      ~{
        lib_name <- .x
        lib_type <- paste0(path, "/", lib_name, "/lib_format_counts.json") %>%
          jsonlite::read_json() %>%
          .$`expected_format`
        dplyr::data_frame(lib, lib_type)
      },
      path
    )
}

#' Retrieve Mapping Rate
#'
#' Retrieve total number of fragments,
#' total mapped reads and
#' mapping rate from a Salmon mapped
#' sample.
#'
#' @param path path to Salmon results files
#'
#' @return A \code{tibble} object
#'
#' @export
#'
retrieve_mapping_rate <- function(path) {
  libs <- list.files(path)
  purrr::map_df(libs, ~{
    salmon_log <- readr::read_lines(paste0(path, "/", .x, "/logs/salmon_quant.log"))
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
    dplyr::data_frame(sample = .x, mapped_reads, total_reads, mapping_pct) %>%
      dplyr::distinct()
  })
}

#' Filter Imported Transcripts
#'
#' Filter salmon quantification files
#' in a \code{txi} format
#' based on metadata column and levels
#'
#' @param tx txi object
#'
#' @param sample_table metadata data frame, containing description
#'                     of each sample and experimental design
#'
#' @param var_column variable to filter the input
#'
#' @param var_levels levels of var_column to filter the input
#'
#' @return A \code{txi} object
#'
#' @export
#'
# tx <- imported_transcripts;sample_table <- metadata_df;var_column <- "population";var_levels <- c("Botucatu2014", "Rio2013")
filter_tx <- function(tx, sample_table, var_column = NULL, var_levels = NULL) {
  first_column <- colnames(sample_table)[1]
  if (!is.null(var_levels)) {
    extracted_cols <- sample_table %>%
      dplyr::filter((!!as.name(var_column)) %in% var_levels) %>%
      dplyr::select(1) %>%
      dplyr::filter((!!as.name(first_column)) %in% colnames(tx$counts)) %>%
      unlist()
  } else {
    extracted_cols <- sample_table %>%
      dplyr::select(1) %>%
      unlist()
  }

  tx$counts <- tx$counts[, extracted_cols]
  tx$abundance <- tx$abundance[, extracted_cols]
  tx$length <- tx$length[, extracted_cols]
  tx
}

#' Perform Differential Expression Analysis
#'
#' Run DESeq2 pipeline for differentially expressed genes.
#'
#' @param tx imported transcripts object containing reads counts
#'
#' @param sample_table metadata data frame, containing description
#'                     of each sample and experimental design
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
#' @param parallel Use all available cores to run.
#'   Uses significantly more RAM. Default = FALSE.
#'
#' @return A \code{tibble} object
#'
#' @importFrom stats as.formula
#'
#' @export
#'
# tx <- tx_rio_2013; sample_table = metadata_df = samples_metadata
# contrast_var <- "condition"; numerator <- "resistant"; denominator <- "susceptible"
de_analysis <- function(tx, sample_table, contrast_var,
                        numerator, denominator,
                        batch_var = NULL,
                        beta_prior = TRUE,
                        force_rep = FALSE,
                        parallel = FALSE) {
  lib_names <- colnames(tx$counts)
  var_column <- colnames(sample_table[1])
  sample_table <- sample_table %>%
    dplyr::filter(!!as.name(var_column) %in% lib_names)
  if (is.null(batch_var)) {
    design_formula <- as.formula(paste0("~", contrast_var))
  } else {
    design_formula <- as.formula(paste0("~", batch_var, " + ", contrast_var))
  }
  dds1 <- DESeq2::DESeqDataSetFromTximport(tx,
    sample_table,
    design = design_formula
  )

  number_of_numerator_samples <- sample_table %>%
    dplyr::filter(!!as.name(contrast_var) == numerator) %>%
    base::nrow()
  number_of_denominator_samples <- sample_table %>%
    dplyr::filter(!!as.name(contrast_var) == denominator) %>%
    base::nrow()
  sample_replicates <- TRUE
  if (number_of_numerator_samples < 2) {
    sample_replicates <- FALSE
  }
  if (number_of_denominator_samples < 2) {
    sample_replicates <- FALSE
  }

  ## START Experimental
  if (isTRUE(force_rep)) {
    dds <- DESeq2::DESeq(dds1, betaPrior = TRUE)
    res <- DESeq2::results(dds,
      alpha = 0.05,
      parallel = parallel,
      contrast = c(contrast_var, numerator, denominator)
    )
    if (!isTRUE("data.frame" %in% class(res))) {
      res <- res %>%
        base::as.data.frame() %>%
        dplyr::as_tibble(rownames = "gene")
    }
    return(res)
  }
  ## END Experimental
  if (isTRUE(sample_replicates)) {
    if (isTRUE(beta_prior)) {
      dds <- DESeq2::DESeq(dds1,
                           parallel = parallel,
                           betaPrior = TRUE
                           )
      res <- DESeq2::results(dds,
        alpha = 0.05,
        parallel = parallel,
        contrast = c(contrast_var, numerator, denominator)
      )
    } else {
      dds <- DESeq2::DESeq(dds1,
                           parallel = parallel,
                           betaPrior = FALSE)
      res <- DESeq2::lfcShrink(dds,
        parallel = TRUE,
        coef = paste0(
          contrast_var, "_",
          numerator,
          "_vs_", denominator
        ),
        parallel = parallel,
        type = "apeglm"
      )
    }
  } else {
    numerator_name <- sample_table %>%
      dplyr::filter(!!as.name(contrast_var) == numerator)
    numerator_name <- as.character(numerator_name[, 1])
    denominator_name <- sample_table %>%
      dplyr::filter(!!as.name(contrast_var) == denominator)
    denominator_name <- as.character(denominator_name[, 1])

    res <- tx$abundance %>%
      dplyr::as_tibble(rownames = "gene") %>%
      dplyr::mutate(`FC` = !!as.name(numerator_name) /
        !!as.name(denominator_name)) %>%
      dplyr::arrange(-`FC`) %>%
      dplyr::mutate(`log2FoldChange` = log2(`FC`)) %>%
      dplyr::mutate(
        log2FoldChange = dplyr::if_else(
          is.finite(`log2FoldChange`), `log2FoldChange`, NA_real_
        )
      ) %>%
      dplyr::mutate(pvalue = NA_real_) %>%
      dplyr::mutate(padj = NA_real_) %>%
      dplyr::select(-(!!numerator_name)) %>%
      dplyr::select(-(!!denominator_name)) %>%
      dplyr::select(-`FC`) %>%
      dplyr::arrange(`gene`)
  }
  if (!isTRUE("data.frame" %in% class(res))) {
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
# de_res <- de_res_rio2013
# gene_sets <- aaegdata::kegg_gene_sets
gsea_analysis <- function(de_res, gene_sets, file_ext = "gmt") {
  de_res %>%
    dplyr::mutate(
      log2FoldChange = dplyr::if_else(is.na(log2FoldChange), 0, log2FoldChange)
    )
  sorted_res <- de_res %>%
    dplyr::select(gene, log2FoldChange) %>%
    dplyr::arrange(log2FoldChange)
  # ordered from strongest down-regulated to strongest upregulated
  ranked_genes <- sorted_res$log2FoldChange
  names(ranked_genes) <- sorted_res$gene
  if (isTRUE("data.frame" %in% class(de_res))) {
    gene_set_col <- base::colnames(gene_sets)[1]
    temp_list <- gene_sets %>%
      dplyr::group_by(!!as.name(gene_set_col)) %>%
      dplyr::summarise(gene = list(gene)) %>%
      dplyr::ungroup()
    gene_sets_list <- temp_list$gene
    names(gene_sets_list) <- temp_list %>%
      base::as.data.frame() %>%
      .[, gene_set_col]
    gene_sets <- gene_sets_list
    if (!is.list(gene_sets)) {
      if (file_ext == "gmt") {
        gene_sets <- fgsea::gmtPathways(gene_sets)
      } else {
        if (file_ext == "rds") {
          gene_sets <- readr::read_rds(gene_sets)
        }
      }
    }
  }
  gene_sets <- purrr::map(gene_sets, unique)
  gsea_res <- fgsea::fgsea(
    pathways = gene_sets,
    stats = ranked_genes,
    nperm = 100000
  )
  gsea_res <- gsea_res %>%
    dplyr::as_tibble() %>%
    dplyr::rename(pvalue = pval) %>%
    dplyr::rename(leading_edge = leadingEdge)
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
le_analysis <- function(gsea_res, set) {
  gsea_res %>%
    tidyr::unnest() %>%
    dplyr::filter(pathway == set) %>%
    dplyr::pull(leading_edge)
}

#' Retrieve Differentially Expressed Genes
#'
#' Extract differentially expressed genes from
#' DESeq2 results based on p-value cuttof.
#'
#' @param de_res de_analysis result
#'
#' @param pvalue p-value cut-off for enrichment or differential consideration,
#'               default = 0.05
#
#' @param fdr logical, default = FALSE should False discovery rate adjusted
#'               pvalue be used instead of pvalue
#'
#' @return a \code{vector} of gene names
#'
#' @export
#'
retrieve_deg <- function(de_res, pvalue = 0.05, fdr = FALSE) {
  alpha <- pvalue
  if (!isTRUE("data.frame" %in% class(de_res))) {
    de_res <- de_res %>%
      base::as.data.frame() %>%
      dplyr::as_tibble(rownames = "gene")
  }
  if (isTRUE(fdr)) {
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
#' @param sample_table metadata data frame, containing description
#'                     of each sample and experimental design
#'
#' @param color_by variable from sample table to group samples in the plot
#'
#' @param num  number of genes to plot, filtering by variance between samples,
#'            if num is negative, it is used the absolute num of genes with
#'            smaller, if num is NULL all genes are used.
#'            default = NULL.
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
plot_heatmap <- function(tx, sample_table, color_by = NULL,
                         num = NULL, scale = "row",
                         show_rownames = FALSE, ...) {
  expression_matrix <- txomics::prepare_tx_mat(tx)
  if (is.null(num)) {
    num <- nrow(expression_matrix)
  }
  mat <- expression_matrix[base::rowSums(expression_matrix) > 1, ]
  ## Scaling matrix
  if (scale == "row") {
    mat <- base::t(base::scale(base::t(mat), center = TRUE, scale = TRUE))
  } else {
    if (scale == "column") {
      mat <- base::scale(mat, center = TRUE, scale = TRUE)
    }
  }
  filtered_mat <- mat %>%
    dplyr::as_data_frame(rownames = "gene") %>%
    dplyr::mutate(row_var = txomics::calc_var_by_row(.[2:length(.)])) %>%
    dplyr::arrange(-row_var)
  if (isTRUE(num > 0)) {
    filtered_mat <- filtered_mat %>%
      utils::head(num)
  }
  if (isTRUE(num < 0)) {
    num <- abs(num)
    filtered_mat <- filtered_mat %>%
      utils::tail(num)
  }
  if (isTRUE(num == 0)) {
    stop("Number of genes can't be zero.")
  }
  filtered_mat <- filtered_mat %>%
    dplyr::select(-row_var)
  mat <- filtered_mat %>%
    dplyr::select(-gene) %>%
    as.matrix()
  rownames(mat) <- filtered_mat$gene
  if (is.null(color_by)) {
    annotation_df <- NA
    annotation_pallete <- NA
  } else {
    annotation_df <- sample_table %>%
      dplyr::select(1, !!color_by)
    row_names <- annotation_df %>% dplyr::pull(1)
    annotation_df <- annotation_df %>%
      dplyr::select(-1) %>%
      base::as.data.frame()
    base::rownames(annotation_df) <- row_names
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
    for (i in colnames(annotation_df)) {
      names(annotation_pallete[[i]]) <- unique(annotation_df[, i])
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
  if (!isTRUE(scale == "none")) {
    color_function <- function() {
      grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(
        n = 7, name =
          "PRGn"
      )))(100)
    }
  } else {
    color_function <- function() {
      viridis::viridis(100, direction = 1)
    }
  }
  pheatmap::pheatmap(
    mat = mat,
    color = color_function(),
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
#' Plot Flat Three Dimensional PCA
#'
#' Plot the results of a pseudo/flat tri-dimensional
#' principal component analysis, using 2 two dimensional plots.
#'
#' @param tx txi object, matrix or data frame with trancript abundance
#'
#' @param sample_table metadata data frame, containing description
#'                     of each sample and experimental design
#'
#' @param color_by variable from sample table to group samples in the plot
#'
#' @param num  number of genes to plot, filtering by variance between samples,
#'            if num is negative, it is used the absolute num of genes with
#'            smaller, if num is NULL all genes are used.
#'            default = NULL.
#'
#' @importFrom stats prcomp
#'
#' @importFrom utils head
#'
#' @export
#'
plot_pca <- function(tx, sample_table,
                     color_by = NULL, num = NULL) {
  expression_matrix <- txomics::prepare_tx_mat(tx)
  if (is.null(num)) {
    num <- nrow(expression_matrix)
  }
  mat <- expression_matrix[rowSums(expression_matrix) > 1, ]
  ## Scaling matrix by row
  mat <- base::t(base::scale(base::t(mat), center = TRUE, scale = TRUE))

  filtered_mat <- mat %>%
    dplyr::as_data_frame(rownames = "gene") %>%
    dplyr::mutate(row_var = txomics::calc_var_by_row(.[2:length(.)])) %>%
    dplyr::arrange(-row_var)
  if (isTRUE(num > 0)) {
    filtered_mat <- filtered_mat %>%
      utils::head(num)
  }
  if (isTRUE(num < 0)) {
    num <- abs(num)
    filtered_mat <- filtered_mat %>%
      utils::tail(num)
  }
  if (isTRUE(num == 0)) {
    stop("Number of genes can't be zero.")
  }
  filtered_mat <- filtered_mat %>%
    dplyr::select(-row_var)
  mat <- filtered_mat %>%
    dplyr::select(-gene) %>%
    base::as.matrix()
  base::rownames(mat) <- filtered_mat$gene

  pca <- mat %>%
    base::t() %>%
    stats::prcomp()

  var_percent <- pca$sdev^2 / base::sum(pca$sdev^2)

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
  pc_3 <- pca$x[, 3]
  pc_3_df <- dplyr::data_frame(
    sample = names(pc_3),
    PC3 = pc_3
  )
  if (is.null(color_by)) {
    ## plot without colors
    df <- pc_1_df %>%
      dplyr::left_join(pc_2_df, by = "sample") %>%
      dplyr::left_join(pc_3_df, by = "sample") %>%
      tidyr::gather(pcomp, value, -c(sample, PC1))
    breaks_x <- base::pretty(df$value)
    breaks_y <- base::pretty(df$PC1)
    plot_obj <- df %>%
      ggplot2::ggplot(ggplot2::aes(x = value, y = PC1, color = sample)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::xlab(paste0(
        "PC2: ", round(var_percent[2] * 100), "% variance",
        "\t  PC3: ", round(var_percent[3] * 100),
        "% variance"
      )) +
      ggplot2::ylab(paste0(
        "PC1: ", round(var_percent[1] * 100),
        "% variance"
      )) +
      ggplot2::facet_grid(cols = ggplot2::vars(pcomp), scales = "free") +
      ggplot2::scale_color_viridis_d(direction = 1) +
      ggplot2::annotate("segment",
        x = -Inf, xend = -Inf,
        y = -Inf, yend = Inf, color = "black"
      ) +
      ggplot2::annotate("segment",
        x = -Inf, xend = -Inf,
        y = -Inf, yend = Inf, color = "black"
      ) +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text.x = ggplot2::element_blank()
      ) +
      ggplot2::scale_x_continuous(
        breaks = breaks_x,
        limits = base::range(breaks_x)
      ) +
      ggplot2::scale_y_continuous(
        breaks = breaks_y,
        limits = base::range(breaks_y)
      )
    return(plot_obj)
  }
  if (!isTRUE(all(color_by %in% names(sample_table)))) {
    stop("the values of 'color_by' should specify columns of sample_table metadata")
  }
  group_df <- sample_table %>%
    dplyr::select(sample, !!color_by)

  df <- pc_1_df %>%
    dplyr::left_join(pc_2_df, by = "sample") %>%
    dplyr::left_join(pc_3_df, by = "sample") %>%
    tidyr::gather(pcomp, value, -c(sample, PC1)) %>%
    dplyr::left_join(group_df, by = "sample") %>%
    tidyr::unite(group, !!color_by)

  breaks_x <- base::pretty(df$value)
  breaks_y <- base::pretty(df$PC1)
  df %>%
    ggplot2::ggplot(ggplot2::aes(x = value, y = PC1, color = group)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::xlab(paste0(
      "PC2: ", round(var_percent[2] * 100), "% variance",
      "\t\t\t\t  PC3: ", round(var_percent[3] * 100), "% variance"
    )) +
    ggplot2::ylab(paste0(
      "PC1: ", round(var_percent[1] * 100),
      "% variance"
    )) +
    ggplot2::facet_grid(cols = ggplot2::vars(pcomp), scales = "free") +
    ggplot2::scale_color_viridis_d(direction = 1) +
    ggplot2::annotate("segment",
      x = -Inf, xend = -Inf,
      y = -Inf, yend = Inf, color = "black"
    ) +
    ggplot2::annotate("segment",
      x = -Inf, xend = -Inf,
      y = -Inf, yend = Inf, color = "black"
    ) +
    ggpubr::theme_pubr() +
    # ggplot2::coord_fixed() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_continuous(breaks = breaks_x, limits = base::range(breaks_x)) +
    ggplot2::scale_y_continuous(breaks = breaks_y, limits = base::range(breaks_y))
  # ggplot2::xlab(paste0("PC2: ", round(var_percent[2] * 100), "% variance")) +
  #  ggplot2::ylab(paste0("PC1: ", round(var_percent[1] * 100), "% variance")) +
  #  ggplot2::scale_color_viridis_d( direction = 1 ) +
  #  ggpubr::theme_pubr()
}
#' Volcano Plot
#'
#' A volcano plot from DE genes analysis results.
#'
#' @param de_res de_analysis results
#'
#' @param pvalue p-value cut-off for enrichment or differential consideration,
#'               default = 0.05
#'
#' @param lfc_threshold absolute value of log2FoldChange to limit the significance
#'
#' @param fdr logical, default = FALSE should False discovery rate adjusted
#'               pvalue be used instead of pvalue
#'
#' @return A \code{ggplot} object
#'
#' @export
#'
plot_volcano <- function(de_res, pvalue = 0.05,
                         lfc_threshold = NULL, fdr = FALSE) {
  # scatter_plot <- function(data, x, y) {
  #   x <- enquo(x)
  #   y <- enquo(y)
  #
  #   ggplot(data) + geom_point(aes(!!x, !!y))
  # }
  # scatter_plot(mtcars, disp, drat)
  alpha <- pvalue
  pvalue_var <- dplyr::quo(pvalue)
  if (isTRUE(fdr)) {
    pvalue_var <- dplyr::quo(padj)
  }
  df <- de_res %>%
    dplyr::select(gene, log2FoldChange, pvalue, padj)
  if (isTRUE(fdr)) {
    df <- df %>%
      dplyr::filter(!is.na(padj)) %>%
      dplyr::mutate(sig_alpha = dplyr::if_else(padj < alpha, 1, 0))
  } else {
    df <- df %>%
      dplyr::filter(!is.na(pvalue)) %>%
      dplyr::mutate(sig_alpha = dplyr::if_else(pvalue < alpha, 1, 0))
  }
  if (is.null(lfc_threshold)) {
    if (isTRUE(base::nrow(dplyr::filter(df, sig_alpha == 1)) == 0)) {
      lfc_cutoff <- 2
    } else {
      lfc_cutoff <- df %>%
        dplyr::filter(sig_alpha == 1) %>%
        dplyr::pull(log2FoldChange) %>%
        abs() %>%
        max()
    }
  } else {
    lfc_cutoff <- lfc_threshold
  }
  df <- df %>%
    dplyr::mutate(sig_lfc = dplyr::if_else(
      abs(log2FoldChange) >= lfc_cutoff, 2, 0
    ))
  df <- df %>%
    dplyr::mutate(sig = as.character(sig_alpha + sig_lfc)) %>%
    dplyr::select(-c(sig_alpha, sig_lfc)) %>%
    dplyr::filter(!is.na(log2FoldChange))
  pvalue_line <- alpha
  if (isTRUE(fdr)) {
    pvalue_line <- df %>%
      dplyr::filter(padj < alpha)
    if (isTRUE(nrow(pvalue_line) == 0)) {
      pvalue_line <- df %>%
        dplyr::pull(pvalue) %>%
        base::min()
    } else {
      pvalue_line <- pvalue_line %>%
        dplyr::pull(pvalue) %>%
        base::max()
    }
  }
  if (isTRUE(max(abs(df$log2FoldChange)) <= abs(lfc_cutoff))) {
    breaks_x <- base::pretty(base::range(-lfc_cutoff, lfc_cutoff))
  } else {
    breaks_x <- base::pretty(df$log2FoldChange)
  }
  if (isTRUE(max(-log10(df$pvalue)) <= 4)) {
    range_y <- c(0, 4)
  } else {
    range_y <- -log10(df$pvalue)
  }
  breaks_y <- base::pretty(range_y)
  pvalue_cutoff <- pvalue
  sig_text <- c(
    paste0(
      as.character(pvalue_var)[2],
      " > ", as.character(pvalue_cutoff)
    ),
    paste0(
      as.character(pvalue_var)[2],
      " <= ", as.character(pvalue_cutoff)
    ),
    paste0(
      "abs(log2FC) > ",
      as.character(lfc_cutoff)
    ),
    paste0(
      as.character(pvalue_var)[2], " <= ",
      as.character(pvalue_cutoff), " & ",
      "abs(log2FC) >= ", as.character(lfc_cutoff)
    )
  )
  df <- df %>%
    dplyr::mutate(sig = dplyr::if_else(
      sig == "0", sig_text[1], sig
    )) %>%
    dplyr::mutate(sig = dplyr::if_else(
      sig == "1", sig_text[2], sig
    )) %>%
    dplyr::mutate(sig = dplyr::if_else(
      sig == "2", sig_text[3], sig
    )) %>%
    dplyr::mutate(sig = dplyr::if_else(
      sig == "3", sig_text[4], sig
    ))
  df$sig <- factor(df$sig, levels = sig_text)
  color_pallete <- viridis::viridis(4)
  names(color_pallete) <- sig_text
  p <- df %>%
    ggplot2::ggplot(ggplot2::aes(log2FoldChange, -log10(pvalue))) +
    ggplot2::geom_hline(yintercept = -log10(pvalue_line), alpha = 0.5) +
    ggplot2::geom_vline(xintercept = lfc_cutoff, alpha = 0.5) +
    ggplot2::geom_vline(xintercept = -lfc_cutoff, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = sig)) +
    ggplot2::scale_color_manual(
      values = color_pallete
    ) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(df, sig %in% sig_text[4]),
      ggplot2::aes(label = gene)
    ) +
    ggplot2::scale_x_continuous(
      breaks = breaks_x, limits = base::range(breaks_x * 1.1)
    ) +
    ggplot2::scale_y_continuous(
      breaks = breaks_y, limits = base::range(breaks_y)
    ) +
    ggpubr::theme_pubr() +
    ggplot2::xlab(~Log[2] ~ "(FC)") +
    ggplot2::ylab(~-Log[10] ~ "(p-value)") +
    ggplot2::labs(color = " ")
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
# gsea_res <- gsea_res_batch; de_res <- de_res_rio_2013; gene_sets <- aaegdata::kegg_gene_sets;
# top_sets <- 10; enrichment <- FALSE; title = TRUE
plot_gsea <- function(gsea_res, de_res, gene_sets,
                      top_sets = 10, enrichment = FALSE,
                      title = NULL) {
  de_res$log2FoldChange[is.na(de_res$log2FoldChange)] <- 0
  sorted_res <- de_res %>%
    dplyr::select(gene, log2FoldChange) %>%
    dplyr::arrange(log2FoldChange)
  # ordered from strongest down-regulated to strongest upregulated
  ranked_genes <- sorted_res$log2FoldChange
  names(ranked_genes) <- sorted_res$gene
  if (isTRUE("data.frame" %in% class(de_res))) {
    gene_set_col <- base::colnames(gene_sets)[1]
    temp_list <- gene_sets %>%
      dplyr::group_by(!!as.name(gene_set_col)) %>%
      dplyr::summarise(gene = list(gene)) %>%
      dplyr::ungroup()
    gene_sets_list <- temp_list$gene
    names(gene_sets_list) <- temp_list %>%
      base::as.data.frame() %>%
      .[, gene_set_col]
    gene_sets <- gene_sets_list
    if (!is.list(gene_sets)) {
      if (file_ext == "gmt") {
        gene_sets <- fgsea::gmtPathways(gene_sets)
      } else {
        if (file_ext == "rds") {
          gene_sets <- readr::read_rds(gene_sets)
        }
      }
    }
  }
  gene_sets <- purrr::map(gene_sets, unique)

  if (isTRUE(enrichment)) {
    if (isTRUE(title)) {
      title <- top_sets
    }
    fgsea::plotEnrichment(
      gene_sets[[top_sets]],
      ranked_genes
    ) +
      ggplot2::labs(title = title) +
      ggpubr::theme_pubr()
    # labs(title="Programmed Cell Death")
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
#' @param pvalue p-value cut-off for enrichment or differential consideration,
#'               default = 0.05
#'
#' @param fdr logical, default = FALSE should False discovery rate adjusted
#'               pvalue be used instead of pvalue
#'
#' @param names vector with names to be used for each group
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
plot_venn_diagram <- function(..., pvalue = 0.05, fdr = FALSE, names = NULL,
                              filename, width = 9, height = NULL) {
  file_extension <- filename %>%
    stringr::str_extract("\\.[a-zA-Z]*$")
  if (!isTRUE(identical(file_extension, ".pdf"))) {
    # remove to introduce another extensions
    stop("filename is not PDF")
  }
  res_list <- list(...)
  if (isTRUE(length(res_list) == 1)) {
    if (!is.list(res_list)) {
      stop("can't plot a venn diagram for one group")
    }
  }
  if (!is.null(names(res_list))) {
    names <- names(res_list)
  }
  if (is.null(names)) {

  } else {
    if (isTRUE(identical(length(res_list), length(names)))) {
      names(res_list) <- names
    } else {
      stop("length of names should be the same of the objects passed to '...'")
    }
  }
  file_path <- filename %>% stringr::str_replace(file_extension, "")
  file_name_svg <- paste0(file_path, ".svg")

  pvalue_var <- dplyr::quo(pvalue)
  if (isTRUE(fdr)) {
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
  if (is.null(height)) {
    output_height <- width
  }
  venn_width <- min(output_height, output_width)
  VennDiagram::venn.diagram(
    x = input_list,
    filename = file_name_svg, imagetype = "svg",
    height = venn_width, width = venn_width, col = "transparent",
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
    label.col = c(
      "orange", "white", "darkorchid4", "white", "white",
      "white", "white", "white", "darkblue", "white", "white", "white", "white",
      "darkgreen", "white"
    ),
    fontfamily = "serif", fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"), cat.cex = 1.5,
    cat.pos = 0, cat.dist = 0.07, cat.fontfamily = "serif", alpha = 0.50, cex = 1.5,
    rotation.degree = 270, margin = 0.2
  )
  if (isTRUE(identical(file_extension, ".pdf"))) {
    file_name_pdf <- paste0(file_path, ".pdf")
    rsvg::rsvg_pdf(svg = file_name_svg, file = file_name_pdf)
  }
}
#' Plot a Heatmap From gsea_analysis
#'
#' Plot a heatmap of enriched gene sets from gsea_analysis
#'
#' @param ... list of gsea_analisys results
#'
#' @param pvalue p-value cut-off for enrichment or differential consideration,
#'               default = 0.05
#'
#' @param fdr logical, default = FALSE should False discovery rate adjusted
#'               pvalue be used instead of pvalue
#'
#' @param names vector with names to be used at the x axis label
#'
#' @param occurrence Number of samples that the gene set
#'                    needs to be enriched
#'
#' @export
#'
# gsea_res_list <- list(gsea_go_bot_2014,gsea_go_neo_2014,gsea_go_rio_2013,gsea_go_rio_2015)
# pvalue = 0.2; fdr = TRUE; names = c("Botucatu 2014","Neopolis 2014","Rio 2013","Rio 2015")
#  occurrence = NULL
plot_gsea_heatmap <- function(..., pvalue = 0.05, fdr = FALSE, names = NULL,
                              occurrence = NULL) {
  gsea_res_list <- list(...)
  pvalue_cutoff <- pvalue
  sample_names <- names
  pvalue_var <- dplyr::quo(pvalue)
  if (isTRUE(fdr)) {
    pvalue_var <- dplyr::quo(padj)
  }
  if (is.null(occurrence)) {
    occurrence <- length(gsea_res_list)
  }
  enriched_gene_sets <- gsea_res_list %>%
    purrr::map(~{
      .x %>%
        dplyr::filter(!!pvalue_var < pvalue_cutoff) %>%
        dplyr::pull(1)
    })
  enriched_vector <- enriched_gene_sets %>%
    unlist() %>%
    table()
  top_sets <- dplyr::data_frame(
    gene_sets = names(enriched_vector),
    occurrency = enriched_vector
  ) %>%
    dplyr::arrange(-occurrency)

  top_sets <- top_sets %>%
    dplyr::filter(occurrency >= occurrence) %>%
    dplyr::pull(1)

  gene_set_var <- gsea_res_list[[1]] %>% colnames()
  gene_set_var <- gene_set_var[1]
  gene_set_var <- dplyr::sym(gene_set_var)

  pvalue_df <- gsea_res_list %>%
    purrr::map2_df(sample_names, ~{
      .x %>%
        dplyr::filter(!!gene_set_var %in% top_sets) %>%
        dplyr::select(1, !!pvalue_var) %>%
        dplyr::mutate(sample_name = .y)
    })

  mat_df <- pvalue_df %>%
    tidyr::spread(sample_name, !!pvalue_var, -1)
  mat_rownames <- mat_df %>%
    dplyr::select(1) %>%
    dplyr::pull()
  mat <- mat_df %>%
    dplyr::select(-1) %>%
    base::as.matrix()
  rownames(mat) <- mat_rownames
  mat <- mat %>%
    base::log10() %>%
    base::abs()
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

  # Because -log10(0.001) = 3 and
  # if p-value is smaller than 0.01, it's already much significant
  p_value_max <- 3
  # p_value_max <- max(mat)

  pheatmap::pheatmap(
    mat = mat,
    color = viridis::viridis(100, direction = 1),
    show_rownames = TRUE,
    border_color = NA,
    treeheight_col = 30,
    treeheight_row = 30,
    cluster_cols = mat_cluster_cols,
    cluster_rows = mat_cluster_rows,
    breaks = seq(0, p_value_max, (p_value_max / 100))
  )
}

#' Plot Gene Expression
#'
#' Plot gene expression points
#'
#' @param gene gene or vector of genes
#'
#' @param tx \code{txi} object, object imported with
#'            `import_tx()` or data frame wih gene abundance
#'
#' @param sample_table metadata data frame, containing description
#'                     of each sample and experimental design
#'
#' @param color_by variable from sample table to group samples in the plot
#'
#' @param filter samples to filter, dafault = NULL
#'
#' @export
#'
# tx=imported_expression;gene="AAEL000080";color_by="day"
# tx <- imported_transcripts;gene <- "AAEL006469";sample_table <- samples_metadata;color_by_var <- "condition";
plot_gene_expression <- function(gene, tx, sample_table = NULL,
                                 color_by = NULL, filter = NULL) {
  if (isTRUE(is.list(tx))) {
    expr_df <- tx$abundance
  }
  metadata_df <- sample_table
  if (is.null(sample_table)) {
    metadata_df <- dplyr::data_frame()
  } else {
    sample_var <- metadata_df %>%
      dplyr::select(1) %>%
      names()
  }
  if (!is.null(color_by)) {
    color_by_var <- color_by
  }
  gene_query <- gene
  expr_df %>%
    dplyr::as_tibble(rownames = "gene") %>%
    dplyr::filter(gene %in% gene_query) %>%
    tidyr::gather(sample, tpm, -gene) %>%
    dplyr::left_join(metadata_df, by = c("sample" = sample_var)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(
      ggplot2::aes(sample, tpm, color = !!as.name(color_by_var)),
      size = 3
    ) +
    ggpubr::theme_pubr() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    )) +
    ggplot2::scale_color_viridis_d()
}

#' Plot GSEA Results and Leading Edge
#'
#' plot gsea results with leading edge genes for enriched gene sets
#'
#' @param gsea_res result from gsea_analysis
#'
#' @param pvalue p-value cut-off for enrichment or differential consideration,
#'               default = 0.05
#'
#' @param fdr logical, default = FALSE should False discovery rate adjusted
#'               pvalue be used instead of pvalue
#'
#' @export
#'
# gsea_res <- gsea_kegg_res; pvalue <- 0.2; fdr = TRUE;
# gsea_res <- gsea_kegg_res; pvalue <- 0.05; fdr = FALSE;
plot_gsea_res <- function(gsea_res, pvalue = 0.05, fdr = FALSE) {
  pvalue_cutoff <- pvalue
  pvalue_var <- dplyr::quo(pvalue)
  pvalue_label <- dplyr::quo(`P-value`)
  if (isTRUE(fdr)) {
    pvalue_var <- dplyr::quo(padj)
    pvalue_label <- dplyr::quo(`Adjusted p-value`)
  }
  gsea_res_sig <- gsea_res %>%
    dplyr::filter(!!pvalue_var < pvalue_cutoff)
  leading_edge_number <- c()
  for (paths in gsea_res_sig$leading_edge) {
    leading_edge_number <- c(leading_edge_number, length(paths))
  }
  gsea_res_df <- dplyr::data_frame(
    gene_set = gsea_res_sig %>% dplyr::pull(1),
    pvalue_temp = gsea_res_sig %>% dplyr::pull(!!pvalue_var),
    NES = gsea_res_sig$NES,
    `Number of LE genes` = leading_edge_number
  )  %>%
    dplyr::arrange(NES)

  if (isTRUE(fdr)){
    gsea_res_df <- gsea_res_df %>%
      dplyr::rename(`Adjusted p-value` = pvalue_temp)
  } else {
    gsea_res_df <- gsea_res_df %>%
      dplyr::rename(`P-value` = pvalue_temp)
  }

  # Plotting
  plot <- gsea_res_df %>%
    ggplot2::ggplot(ggplot2::aes(NES, gene_set)) +
      ggplot2::geom_point(
      ggplot2::aes(colour = !!pvalue_label, size = `Number of LE genes`)) +
      ggplot2::scale_color_gradientn(
        colours = viridis::viridis(4, direction = -1),
        limits = c(0, pvalue_cutoff)
      ) +
      ggplot2::geom_vline(xintercept = 0, size = 0.5, colour = "gray50") +
      ggplot2::theme(
        panel.background = ggplot2::element_rect(
          fill = "gray95", colour = "gray95"
        ),
        panel.grid.major = ggplot2::element_line(
          size = 0.25, linetype = "solid", colour = "gray90"
        ),
        panel.grid.minor = ggplot2::element_line(
          size = 0.25, linetype = "solid", colour = "gray90"
        ),
        axis.title.y = ggplot2::element_blank()
      ) +
      ggplot2::expand_limits(x = c(-3, 3)) +
      ggplot2::scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
      ggplot2::scale_y_discrete(limits = rev(gsea_res_df$gene_set)) +
      ggplot2::ylab( "Gene sets") +
      ggpubr::theme_pubr()
  plot
}

#' Plot Leading Edge Heatmap
#'
#' plot Ledading Edge analysis results enriched gene sets
#'
#' @param gsea_res result from gsea_analysis
#'
#' @param pvalue p-value cut-off for enrichment or differential consideration,
#'               default = 0.05
#'
#' @param fdr logical, default = FALSE should False discovery rate adjusted
#'               pvalue be used instead of pvalue
#'
#' @param num  number of genes to plot, based on the leading edge occurence,
#'            if num is negative, it is used the absolute num of genes with
#'            lesser occurence, if num is NULL all genes are used.
#'            default = NULL.
#' @export
#'
# gsea_res <- gsea_kegg_res; pvalue <- 0.2; fdr = TRUE; num <- 30;
# gsea_res <- gsea_reactome_res; pvalue <- 0.05; fdr = FALSE; num <- -20;
plot_le_heatmap <- function(gsea_res, pvalue = 0.05, fdr = FALSE,
                            num = NULL) {
  pvalue_cutoff <- pvalue
  pvalue_var <- dplyr::quo(pvalue)
  if (isTRUE(fdr)) {
    pvalue_var <- dplyr::quo(padj)
  }
  gsea_res_sig <- gsea_res %>%
    dplyr::filter(!!pvalue_var < pvalue_cutoff)
  le_df <- gsea_res_sig %>%
    dplyr::select(pathway, leading_edge) %>%
    tidyr::unnest()

  le_occurence_df <- le_df %>%
    dplyr::group_by(leading_edge) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(-n)

  if (is.null(num)) {
    le_genes <- unique(le_occurence_df$leading_edge)
  } else if (isTRUE(num < 0)) {
    le_genes <- le_occurence_df %>%
      dplyr::slice((dplyr::n() - num):dplyr::n()) %>%
      dplyr::pull(leading_edge)
  } else if (isTRUE(num > 0)) {
    le_genes <- le_occurence_df %>%
      dplyr::slice(1:num) %>%
      dplyr::pull(leading_edge)
  } else {
    stop("`num` should be an integer different from 0 or NULL")
  }

  le_wide_df <- le_df %>%
    dplyr::filter(leading_edge %in% le_genes) %>%
    tidyr::spread(leading_edge, -pathway) %>%
    dplyr::distinct()
  le_mat_row_names <- le_wide_df$pathway
  le_wide_df <- le_wide_df %>%
    dplyr::select(-pathway) %>%
    as.data.frame()
  for (column in colnames(le_wide_df)) {
    le_wide_df[[column]] <- dplyr::if_else(!is.na(le_wide_df[[column]]), 1, 0)
  }
  rownames(le_wide_df) <- le_mat_row_names

  utils::assignInNamespace(
    x = "draw_colnames",
    value = "draw_colnames_45",
    ns = asNamespace("pheatmap")
  )

  pheatmap::pheatmap(
    mat = le_wide_df,
    color = viridis::viridis(2, direction = 1),
    show_rownames = TRUE,
    #border_color = NA,
    treeheight_col = 0,
    treeheight_row = 0,
    drop_levels = TRUE,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    legend = TRUE,
    legend_breaks = c(0.25, 0.75),
    legend_labels = c("Not present", "Leading Edge")
  )
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
multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
  # library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(base::list(...), plotlist)

  num_plots <- base::length(plots)

  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- base::matrix(seq(1, cols * base::ceiling(num_plots / cols)),
      ncol = cols, nrow = base::ceiling(num_plots / cols)
    )
  }

  if (num_plots == 1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), base::ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:num_plots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- base::as.data.frame(base::which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(
        layout.pos.row = matchidx$row,
        layout.pos.col = matchidx$col
      ))
    }
  }
}

###############################################
##                                            #
##              Aux Functions                 #
##                                            #
###############################################
#' Pipe, as in Pipeline
#'
#' Exported from dplyr, use \code{\%>\%} to turn
#' function composition into a series of imperative statements.
#'
#' @importFrom dplyr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs,rhs An object and a function to apply to it
#' @examples
#' # Instead of
#' abs(nrow(head(mtcars)))
#' # you can write
#' mtcars %>% head() %>% nrow() %>% abs()
NULL

#' plot aux function
#'
#' Auxiliar function for rotating plot_heatmap x axis
#'
#' @param coln number of columns
#' @param gaps gaps
#' @param ... additional parameters
#' @export
draw_colnames_45 <- function(coln, gaps, ...) {
  # Correct labels orientation, thanks to https://slowkow.com/notes/heatmap-tutorial/
  find_coordinates_temp <- utils::getFromNamespace("find_coordinates", "pheatmap")
  coord <- find_coordinates_temp(length(coln), gaps)
  x <- coord$coord - 0.5 * coord$size
  res <- grid::textGrob(
    coln,
    x = x, y = grid::unit(1, "npc") - grid::unit(3, "bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}

#' Calculate Variance by Row
#'
#' Calculate variance by row in a matrix or data frame
#'
#' @param x matrix or data frame
#' @param ... additional paramaters addressed to rowSums or rowMeans
#' @export
calc_var_by_row <- function(x, ...) {
  base::rowSums((x - base::rowMeans(x, ...))^2, ...) / (base::dim(x)[2] - 1)
}

#' plot aux function 2
#'
#' Prepare data for plotting PCA or Expression heatmap
#' @param tx expression data or result from import_tx
#' @export
prepare_tx_mat <- function(tx) {
  if (is.list(tx)) {
    if (!dplyr::is.tbl(tx) & !is.data.frame(tx)) {
      expression_matrix <- tx$abundance
    } else {
      expression_matrix <- tx
    }
  } else {
    expression_matrix <- tx
  }
  if (!is.matrix(expression_matrix)) {
    if (!isTRUE(all(!is.na(suppressWarnings(as.numeric(as.matrix(expression_matrix))))))) {
      temp_mat <- base::as.matrix(expression_matrix[, 2:length(expression_matrix)])
      base::rownames(temp_mat) <- expression_matrix %>%
        dplyr::pull(1)
      expression_matrix <- temp_mat
    }
  }
  expression_matrix
}

#' Differential Expression results
#'
#' Filter upregulated and downregulated genes from
#' Differential expression analysis in a table
#'
#' @param ... de_analysis or gsea_analysis results data frames
#'
#' @param pvalue p-value cut-off for enrichment or differential consideration,
#'               default = 0.05
#'
#' @param fdr logical, default = FALSE should False discovery rate adjusted
#'               pvalue be used instead of pvalue
#'
#' @param names vector with names to be used for each group
#'
#' @export
retrieve_deg_table <- function(..., pvalue = 0.05, fdr = FALSE, names = NULL) {
  pvalue_var <- dplyr::quo(pvalue)
  if (isTRUE(fdr)) {
    pvalue_var <- dplyr::quo(padj)
  }
  pvalue_cutoff <- pvalue
  temp_var <- list(...)
  if (is.list(temp_var)) {
    deg_table_population <- temp_var
  }
  else {
    deg_table_population <- list(temp_var)
  }
  if (!is.null(base::names(deg_table_population))) {
    names <- base::names(deg_table_population)
  }
  deg_table_population <- deg_table_population %>%
    purrr::map_df(~{
      temp_res <- .x
      up_deg <- temp_res %>%
        dplyr::filter(!!pvalue_var < pvalue_cutoff) %>%
        dplyr::filter(`log2FoldChange` > 0) %>%
        nrow()
      down_deg <- temp_res %>%
        dplyr::filter(!!pvalue_var < pvalue_cutoff) %>%
        dplyr::filter(`log2FoldChange` < 0) %>%
        nrow()
      dplyr::data_frame(
        up = up_deg,
        down = down_deg,
        total = up_deg + down_deg
      )
    })
  deg_table_population %>%
    dplyr::mutate(population = names) %>%
    dplyr::select(population, up, down, total)
}
#' Leading Edge Genes By Gene Set
#'
#' Retrieve table with leading edge genes by gene set, for multiple
#' experiments the column occurrence will have the number of hits
#'
#' @param ... gsea_analysis results data frames
#'
#' @param names vector with names to be used for each group
#'
#' @export
# gsea_res_list <- list(gsea_complex_bot_2014,gsea_complex_neo_2014,gsea_complex_rio_2013,gsea_complex_rio_2015)
# names <- names <- c("Botucatu 2014", "Neopolis 2014", "Rio 2013", "Rio 2015")
retrieve_le_table <- function(..., names = NULL) {
  gsea_res_list <- list(...)
  if (is.null(names)) {
    names <- ""
  }
  le_res <- gsea_res_list %>%
    purrr::map2_df(names, ~{
      gene_set_var <- .x %>%
        dplyr::select(1) %>%
        names()
      sample_temp <- .y
      .x %>%
        tidyr::unnest() %>%
        dplyr::select(!!gene_set_var, `leading_edge`) %>%
        dplyr::mutate(sample = sample_temp)
    })
  le_res %>%
    dplyr::group_by(!!as.name(gene_set_var), `leading_edge`) %>%
    dplyr::summarise(occurrence = n()) %>%
    dplyr::arrange(!!as.name(gene_set_var), -occurrence) %>%
    dplyr::ungroup()
}

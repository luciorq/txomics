# txomics 0.1.9

* added kallisto to `import_tx()`.
* `gsea_analysis()` now has an option to use fold-change and pvalue to rank elements.
* added `plot_le_heatmap()` function.
* added `plot_deg_heatmap()` function.
* added `label` to `plot_pca()`. **BROKEN**
* added `max` to `plot_volcano()`.

# txomics 0.1.8

* handlers for expression tables in plots.
* added `plot_gene_expression()`.
* added PC3 to `plot_pca()`, for a pseudo 3d view.
* using `fs` for consistent cross-platform file operations.
* `import_tx()` now accepts custom gene to transcripts conversion table.

# txomics 0.1.7

* Added `plot_gene_abundance()`.
* Added `plot_gsea_res()`.
* Added `retrieve_le_table()`.

# txomics 0.1.6

* `txomics` now exports `dplyr` pipe operator `%>%`.
* Some consistency changes regarding parameter names.
* `plot_heatmap()` supporting results.

# txomics 0.1.5

* Added a `NEWS.md` file to track changes to the package.

# txomics 0.1.0

* Initial development

# bulkseqviz: Streamlined Visualization for Bulk RNA-Seq Data

# [![GitHub release](https://img.shields.io/github/v/release/MobinKhoramjoo/bulkseqviz?color=brightgreen)](https://github.com/MobinKhoramjoo/bulkseqviz/releases) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`bulkseqviz` is an R package designed to simplify the downstream analysis and visualization of bulk RNA-sequencing data. Built on an S3 object-oriented framework, it wraps complex workflows, including normalization, differential expression (DESeq2), and dimensionality reduction (PCA, UMAP, t-SNE) into intuitive, one-line functions.

## 📦 Installation

You can install the latest version of `bulkseqviz` from GitHub with:

```         
# install.packages("devtools")
devtools::install_github("MobinKhoramjoo/bulkseqviz")
library(bulkseqviz)
```

## 💡 Real Dataset Example

See `bulkseqviz` in action with a complete case study using public GEO data.

[**View the Tutorial**](https://github.com/MobinKhoramjoo/bulkseqviz/blob/master/vignettes/geo_analysis.md)

## 📄 Cheatsheet

A high level overview of available functions.

![](man/figures/bulkseqviz.png)

## 🛠️ Functions

The first function to run

```         
bs_obj <- create_bulkseqvis_object(counts = counts, # A matrix or data frame of raw counts
                                   metadata = metadata # A data frame of sample information
                                   )

# Print object summary
print(bs_obj)
```

### EDA

#### Overall Count Distribution

```         
overall_count_boxplot(bs_obj,
                      color_by = "condition",
                      group_order = NULL,
                      group_colors = NULL
                      )
```

#### Sample-to-Sample Distance

```         
sample_distance_heatmap(bs_obj, 
                        color_by = "condition", 
                        fontsize = 9,
                        min_gene_counts = 10,
                        colors = NULL,
                        title = "Sample-to-sample distance (Sorted by Group)",
                        ...)
```

### Dimensionality Reduction

#### PCA (2D)

```         
# 2D PCA colored by condition, shaped by batch
plot_pca_2d(bs_obj, 
            color_by = "condition",
             ellipse_by = NULL,
             subset_samples = NULL,
             min_gene_counts = 10,
             ntop = 500,
             pcs = c(1, 2),
             ellipse_level = 0.99,
             group_colors = NULL
             )
```

#### PCA (3D)

```         
# Generates an interactive plotly object
plot_pca_3d(bs_obj, 
            color_by = "condition",
            ntop = 500,
            colors = NULL,
            size = 5
            )
```

#### UMAP

```         
plot_umap_2d(bs_obj, 
             color_by = "condition",
             shape_by = NULL,
             subset_samples = NULL,
             min_gene_counts = 100,
             ntop = 1000,
             umap_seed = 123,
             n_neighbors = 15,
             group_colors = NULL
             )
```

#### t-SNE

```         
plot_tsne_2d(bs_obj, 
             color_by = "condition",
             shape_by = NULL,
             subset_samples = NULL,
             min_gene_counts = 100,
             ntop = 1000,
             tsne_perplexity = 20,
             tsne_seed = 123,
             group_colors = NULL
             )
```

### Differential Expression Analysis (DEA)

We perform Differential Expression (DE) analysis using `DESeq2` wrapped inside the `DEG()` function.

```         
# Contrast 1: Treatment A vs Control
bs_obj <- DEG(bs_obj, 
              design_col = "condition", 
              compare_levels = c("TreatA", "Control"),
              biomart_dataset = NULL, # use one of BioMart databases to map IDs to symbols for visualization.
              ...)

# Contrast 2: Treatment B vs Control
bs_obj <- DEG(bs_obj, 
              design_col = "condition", 
              compare_levels = c("TreatB", "Control"),
              biomart_dataset = NULL, # use one of BioMart databases to map IDs to symbols for visualization.
              ...)
```

### Visualizing DEA Results

#### Volcano Plots

```         
plot_volcano(bs_obj, 
             comparison_id = "TreatA_vs_Control", 
             fc_thresh = 2,
             padj_thresh = 0.05,
             repel_force = 1.5,
             repel_max_overlaps = 15,
             top_labels = 0,
             plot_title = NULL
             )
```

#### Summary of Differential Expression (Barplot)

```         
plot_deg_bar(bs_obj,
             method = "padj",
             cutoff = 0.05,
             fc_cutoff = 1,
             custom_colors = NULL
             )
```

#### Comparing Contrasts (Log2FC vs Log2FC)

```         
plot_fcvsfc(bs_obj, 
            name1 = "TreatA_vs_Control", 
            name2 = "TreatB_vs_Control",
            fc_cutoff = 1, 
            padj_cutoff = 0.05
            )
```

#### Log2FC Dotplot

```         
# Select genes to visualize: draw the gene symbol/id into a vector
library(dplyr)
target_genes <- bs_obj$DE_results$TreatA_vs_Control %>% 
  dplyr::arrange(desc(log2FoldChange)) %>% 
  dplyr::slice_head(n=5) %>% 
  dplyr::pull(gene_id)

plot_fc_dotplot(bs_obj, 
                genes = target_genes,
                padj_cutoff = 0.05, 
                id_col = NULL
                )
```

### Single Gene Expression (Boxplot)

```         
plot_gene_boxplot(bs_obj, 
                  genes = target_genes, 
                  x_var = "condition", 
                  region_var = NULL,
                  x_order = NULL,
                  region_order = NULL,
                  x_colors = NULL,
                  jitter_points = TRUE,
                  normalize_to_baseline = FALSE,
                  log2_transform = FALSE,
                  show_ns = TRUE,
                  biomart_dataset = "hsapiens_gene_ensembl"
                  )
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request or open an Issue for bug reports and feature requests.

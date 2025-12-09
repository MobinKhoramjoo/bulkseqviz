# bulkseqvis: Streamlined Visualization for Bulk RNA-Seq Data

<!-- badges: start -->

<!-- badges: end -->

**bulkseqvis** is an R package designed to simplify the downstream analysis and visualization of bulk RNA-sequencing data. Built on an S3 object-oriented framework, it wraps complex workflows, including normalization, differential expression (DESeq2), and dimensionality reduction (PCA, UMAP, t-SNE) into intuitive, one-line functions.

## 📦 Installation

You can install the development version of bulkseqvis from GitHub with:

```         
# install.packages("devtools") devtools::install_github("MobinKhoramjoo/bulkseqvis") 
```

## 🚀 Key Features

-   **Unified Data Object:** Bundles counts and metadata into a single validated `bulkseq` object.

-   **Differential Expression:** Wrappers for `DESeq2` that handle design formulas and results extraction automatically.

-   **Publication-Ready Plots:** adjustable with `ggplot2`.

## 🏁 Quick Start

Here is a minimal example of how to go from raw counts to a PCA plot.

```         
library(bulkseqvis)  
# 1. Load your data (Example using dummy data) 
counts <- matrix(sample(1:1000, 1000), ncol = 10) 
colnames(counts) <- paste0("Sample", 1:10) 
rownames(counts) <- paste0("Gene", 1:100)  

metadata <- data.frame(   
    condition = rep(c("Control", "Treatment"), each = 5),   
    batch = rep(c("A", "B"), 5) 
) 
    
rownames(metadata) <- colnames(counts)  

# 2. Create the bulkseq object 

bs_obj <- create_bulkseqvis_object(counts, metadata)  

# 3. 2D PCA plot 
plot_pca_2d(bs_obj, color_by = "condition", shape_by = "batch") 
```

## 📚 Documentation

For a full walkthrough of all functions, please refer to the Get Started vignette or run:

```         
vignette("tutorial", package = "bulkseqvis") 
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request or open an Issue for bug reports and feature requests.

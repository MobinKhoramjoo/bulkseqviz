#' Calculate Differential Expression Results
#'
#' Performs differential expression analysis using DESeq2 and stores the results
#' in the bulkseq object. Optionally fetches gene annotations via BioMart.
#'
#' @param bs_obj A \code{bulkseq} object created by \code{create_bulkseqvis_object}.
#' @param design_col Character string. The metadata column to use for the design formula (e.g., "condition").
#' @param compare_levels Character vector of length 2. The levels to compare in the format c("Treatment", "Control").
#' @param covariates Character vector. Optional. Additional metadata columns to include in the design formula.
#' @param min_count Integer. Minimum number of counts required for a gene to be considered expressed. Default 10.
#' @param min_sample Integer. Minimum number of samples that must have at least \code{min_count} reads. Default 3.
#' @param padj_cutoff Numeric. Filter output for padj < this value (for console summary). Default 1.
#' @param log2fc_cutoff Numeric. Filter output for |log2FC| > this value (for console summary). Default 0.
#' @param biomart_dataset Character string. The BioMart dataset to use (e.g., "hsapiens_gene_ensembl"). If NULL, annotation is skipped.
#' @param id_type Character string. The type of gene ID in the count matrix rownames (e.g., "ensembl_gene_id").
#' @param output_dir Character string. Optional. Path to save CSV/rank files. If NULL, files are not saved.
#'
#' @return An updated \code{bulkseq} object with results stored in \code{obj$DE_results}.
#' @export
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
#' @importFrom SummarizedExperiment assay
#' @importFrom edgeR cpm
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom dplyr filter left_join rename distinct arrange select mutate desc
#' @importFrom tibble tibble as_tibble enframe
#' @importFrom stringr str_to_upper
#' @importFrom stats as.formula na.omit
#' @importFrom utils write.csv write.table
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#'   # Basic DE analysis
#'   bs_obj <- DEG(bs_obj, design_col = "condition",
#'                 compare_levels = c("Treated", "Control"))
#'
#'   # With BioMart annotation and file saving
#'   bs_obj <- DEG(bs_obj, design_col = "condition",
#'                 compare_levels = c("Treated", "Control"),
#'                 biomart_dataset = "hsapiens_gene_ensembl",
#'                 output_dir = "./results")
#' }
DEG <- function(bs_obj,
                design_col,
                compare_levels,
                covariates = NULL,
                min_count = 10,
                min_sample = 3,
                padj_cutoff = 1,
                log2fc_cutoff = 0,
                biomart_dataset = "hsapiens_gene_ensembl",
                id_type = "ensembl_gene_id",
                output_dir = NULL) {

  # 1. Validation
  if (!inherits(bs_obj, "bulkseq")) stop("Input must be a 'bulkseq' object.")
  if (!design_col %in% names(bs_obj$metadata)) stop(paste("Column", design_col, "not found in metadata."))

  # Normalize inputs
  metadata <- as.data.frame(bs_obj$metadata)
  counts <- bs_obj$counts

  # Filter metadata to only include the comparison groups
  # This prevents DESeq2 from trying to model groups we aren't interested in
  keep_samples <- metadata[[design_col]] %in% compare_levels
  metadata_sub <- metadata[keep_samples, , drop = FALSE]
  counts_sub <- counts[, keep_samples, drop = FALSE]

  # Set factor levels (Control should be the reference, i.e., second element)
  # DESeq2 uses the first level as reference, so we set compare_levels[2] (Control) as reference
  metadata_sub[[design_col]] <- factor(metadata_sub[[design_col]], levels = rev(compare_levels))

  # 2. Pre-filtering counts
  # Drop genes with NA or low counts
  counts_sub <- counts_sub[rowSums(is.na(counts_sub)) == 0, , drop = FALSE]

  # 3. Build DESeq2 Object
  terms <- c(design_col, covariates)
  design_formula <- stats::as.formula(paste("~", paste(terms, collapse = " + ")))

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_sub,
                                        colData = metadata_sub,
                                        design = design_formula)

  # 4. Low count filter (per group strategy)
  # Keep genes that have >= min_count in at least min_sample samples
  # This logic ensures that a gene is expressed in at least one group (assuming min_sample <= group size)
  # but allows genes to be OFF in the other group (essential for DE analysis).
  keep <- rowSums(SummarizedExperiment::assay(dds) >= min_count) >= min_sample
  dds <- dds[keep, ]

  # 5. Run DESeq
  message("Running DESeq2...")
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  # 6. Results & CPM
  message("Calculating CPM and results...")
  raw_counts <- SummarizedExperiment::assay(dds)
  # Calculate avg CPM using edgeR logic as requested
  avg_cpm <- rowMeans(edgeR::cpm(raw_counts))
  cpms_tbl <- tibble::enframe(avg_cpm, name = "gene_id", value = "avg_cpm")

  # Extract results
  # contrast: c("condition", "Treated", "Control")
  res <- DESeq2::results(dds,
                         contrast = c(design_col, compare_levels[1], compare_levels[2]),
                         cooksCutoff = TRUE,
                         independentFiltering = TRUE)

  res_tbl <- tibble::as_tibble(as.data.frame(res), rownames = "gene_id")

  # 7. BioMart Annotation (Optional)
  if (!is.null(biomart_dataset)) {
    message(paste("Fetching annotation from BioMart:", biomart_dataset))

    ann <- tryCatch({
      mart <- biomaRt::useEnsembl("ensembl", dataset = biomart_dataset)
      # attributes to retrieve
      attrs <- c(id_type, "external_gene_name", "gene_biotype", "description")

      biomaRt::getBM(attributes = attrs,
                     filters = id_type,
                     values = res_tbl$gene_id,
                     mart = mart)
    }, error = function(e) {
      warning(paste("BioMart connection failed or invalid dataset/ID:", e$message))
      return(NULL)
    })

    if (!is.null(ann)) {
      # Rename columns to standard names for joining
      # We assume the first column retrieved is the ID we passed
      colnames(ann)[1] <- "gene_id"
      if("external_gene_name" %in% colnames(ann)) colnames(ann)[which(colnames(ann)=="external_gene_name")] <- "gene_name"

      # Join annotation to results
      # Distinct is needed because biomart can return duplicates
      # Fix: Use .data$gene_id to avoid global variable note
      ann <- dplyr::distinct(ann, .data$gene_id, .keep_all = TRUE)
      res_tbl <- dplyr::left_join(res_tbl, ann, by = "gene_id")
    }
  }

  # 8. Final Join & Sort
  res_tbl <- dplyr::left_join(res_tbl, cpms_tbl, by = "gene_id")
  res_tbl <- dplyr::arrange(res_tbl, .data$padj)

  # 9. Save to Object
  comparison_name <- paste0(compare_levels[1], "_vs_", compare_levels[2])
  bs_obj$DE_results[[comparison_name]] <- res_tbl

  # 10. Write to Disk (Optional)
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    out_prefix <- file.path(output_dir, comparison_name)

    # Save Full CSV
    utils::write.csv(res_tbl, paste0(out_prefix, "_DEGs.csv"), row.names = FALSE)

    # Save Rank File (for GSEA) - requires gene_name and log2FoldChange
    if ("gene_name" %in% colnames(res_tbl) && "log2FoldChange" %in% colnames(res_tbl)) {
      rank_df <- res_tbl %>%
        dplyr::filter(!is.na(.data$gene_name), !is.na(.data$log2FoldChange)) %>%
        dplyr::select("gene_name", "log2FoldChange") %>%
        dplyr::mutate(gene_name = stringr::str_to_upper(.data$gene_name)) %>%
        dplyr::arrange(dplyr::desc(.data$log2FoldChange))

      utils::write.table(rank_df, paste0(out_prefix, "_rank.rnk"),
                         sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
    }
  }

  # 11. Console Summary
  n_sig <- nrow(dplyr::filter(res_tbl, .data$padj < padj_cutoff))
  n_sig_fc <- nrow(dplyr::filter(res_tbl, .data$padj < padj_cutoff & abs(.data$log2FoldChange) > log2fc_cutoff))

  cat("--------------------------------------\n")
  cat(paste("Comparison:", comparison_name, "\n"))
  cat(paste("Total genes:", nrow(res_tbl), "\n"))
  cat(paste("Signif (padj <", padj_cutoff, "):", n_sig, "\n"))
  cat(paste("Signif + FC (|log2FC| >", log2fc_cutoff, "):", n_sig_fc, "\n"))
  cat("--------------------------------------\n")

  return(bs_obj)
}


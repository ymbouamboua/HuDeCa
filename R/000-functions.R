
# Customized R functions
# Author= Yvon Mbouamboua (yvon.mbouamboua@inserm.fr)

# Utility: %||%
`%||%` <- function(a, b) if(!is.null(a)) a else b

#' Remove duplicated cell barcodes within and/or across Seurat objects
#'
#' This function cleans duplicated cell barcodes in a list of Seurat objects.
#' It can remove duplicates occurring **within each object**, **across objects**,
#' or **both**. A detailed log is returned reporting how many cells were removed
#' per object at each step.
#'
#' @param objs Named list of \code{Seurat} objects.
#' @param mode Character string specifying which duplicates to remove:
#'   \itemize{
#'     \item \code{"within"} – remove duplicated barcodes within each object
#'     \item \code{"across"} – remove duplicated barcodes shared across objects
#'     \item \code{"both"} – remove both within- and across-object duplicates (default)
#'   }
#' @param save Logical; whether to save cleaned Seurat objects to disk (default: \code{TRUE}).
#' @param outdir Output directory where cleaned objects will be written.
#' @param suffix Filename suffix for saved objects (default: \code{"_clean.rds"}).
#' @param verbose Logical; whether to print a summary table to the console (default: \code{TRUE}).
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item \code{objs} – cleaned list of Seurat objects
#'     \item \code{log} – data.frame summarizing barcode removal per object
#'   }
#'
#' @details
#' Duplicate removal strategy:
#' \enumerate{
#'   \item Remove duplicated cell barcodes *within* each Seurat object using
#'   \code{duplicated(colnames(object))}.
#'   \item Identify barcodes shared across multiple objects and remove them
#'   from all objects.
#' }
#'
#' The function does not rename barcodes; duplicates are removed entirely.
#' This behavior is appropriate for demultiplexed or merged datasets where
#' barcode reuse indicates technical duplication.
#'
#' @examples
#' \dontrun{
#' res <- clean_seurat_barcodes(
#'   objs,
#'   mode = "both",
#'   outdir = "cleaned_objects"
#' )
#'
#' objs_clean <- res$objs
#' removal_log <- res$log
#' }
#'
#' @seealso \code{\link[Seurat]{Seurat}}, \code{\link[base]{duplicated}}
#'
#' @export
#' 
#' Clean duplicated cell barcodes in Seurat objects
#'
#' Removes duplicated cell barcodes within individual Seurat objects,
#' across multiple Seurat objects, or both.
#'
#' @param objs A Seurat object or a named list of Seurat objects
#' @param mode One of "within", "across", or "both"
#' @param outdir Optional output directory to save cleaned objects
#' @param verbose Logical; print summary messages
#' @param save Logical; save data
#' @return A list with cleaned Seurat objects and a log table
#' @export
#' 
clean_seurat_barcodes <- function(
    objs,
    mode = c("within", "across", "both"),
    outdir = NULL,
    save = FALSE,
    verbose = TRUE
) {
  mode <- match.arg(mode)
  
  # CRITICAL FIX — MUST BE FIRST
  if (inherits(objs, "Seurat")) {
    objs <- list(object = objs)
  }
  
  if (!is.list(objs)) {
    stop("`objs` must be a Seurat object or a list of Seurat objects.")
  }
  
  # Initialize log
  log <- data.frame(
    object = names(objs),
    removed_within = 0L,
    removed_across = 0L,
    kept_cells = 0L,
    stringsAsFactors = FALSE
  )
  
  # Remove duplicates WITHIN each object
  if (mode %in% c("within", "both")) {
    for (i in seq_along(objs)) {
      obj <- objs[[i]]
      cn <- colnames(obj)
      
      dup_within <- duplicated(cn)
      n_removed <- sum(dup_within)
      
      if (n_removed > 0) {
        obj <- obj[, !dup_within]
      }
      
      objs[[i]] <- obj
      log$removed_within[i] <- n_removed
    }
  }
  
  # Remove duplicates ACROSS objects
  if (mode %in% c("across", "both") && length(objs) > 1) {
    
    all_cells <- unlist(lapply(objs, colnames))
    dup_cells <- all_cells[duplicated(all_cells)]
    
    if (length(dup_cells) > 0) {
      dup_cells <- unique(dup_cells)
      
      for (i in seq_along(objs)) {
        obj <- objs[[i]]
        cn <- colnames(obj)
        
        to_remove <- cn %in% dup_cells
        n_removed <- sum(to_remove)
        
        if (n_removed > 0) {
          obj <- obj[, !to_remove]
        }
        
        objs[[i]] <- obj
        log$removed_across[i] <- n_removed
      }
    }
  }
  
  # Final counts
  for (i in seq_along(objs)) {
    log$kept_cells[i] <- ncol(objs[[i]])
  }
  
  # Save cleaned objects
  if (save) {
    if (!is.null(outdir)) {
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      for (nm in names(objs)) {
        saveRDS(objs[[nm]], file.path(outdir, paste0(nm, ".rds")))
      }
    }
  }
  
  if (verbose) {
    message("✔ Barcode cleaning complete")
    print(log)
  }
  
  return(list(
    objects = objs,
    log = log
  ))
}


#' Run Quality Control on Seurat Objects
#'
#' Performs automated cell filtering based on gene/UMI counts, mitochondrial
#' and ribosomal content, dropout rate, and optional doublet detection.
#' Works on single or multiple Seurat objects.
#'
#' @param x Seurat object or list of Seurat objects.
#' @param min_feat Minimum number of detected genes.
#' @param min_umi Minimum number of UMIs.
#' @param mad_n Number of MADs for upper cutoffs (if method = "MAD").
#' @param max_mito Maximum mitochondrial percentage.
#' @param calc_ribo Logical; compute ribosomal gene percentage.
#' @param max_ribo Maximum ribosomal percentage.
#' @param rm_dbl Logical; perform doublet removal using scDblFinder.
#' @param calc_drop Logical; compute dropout fraction.
#' @param max_drop Maximum dropout fraction.
#' @param method Thresholding method: "MAD", "fixed", or "none".
#' @param fixed_thr Named list of fixed upper cutoffs.
#' @param mito_pat Regex for mitochondrial genes.
#' @param ribo_pat Regex for ribosomal genes.
#' @param species "human" or "mouse" (sets default gene patterns).
#' @param log_g2u Logical; filter by log10(genes/UMI).
#' @param min_g2u Minimum genes/UMI ratio.
#' @param sample_col Metadata column for sample ID.
#' @param outdir Output directory.
#' @param merge Logical; merge outputs into a single Seurat object.
#'
#' @return A filtered Seurat object, list, or merged object.
#' @export
#'
#' @examples
#' qc_data <- run_qc(seurat_list, calc_ribo = TRUE, rm_dbl = TRUE)
#'
run_qc <- function(
    x,
    min_feat = 200,
    min_umi = 500,
    mad_n = 5,
    max_mito = 5,
    calc_ribo = FALSE,
    max_ribo = 3,
    calc_drop = FALSE,
    max_drop = 0.95,
    rm_dbl = FALSE,
    method = c("MAD", "fixed", "none"),
    fixed_thr = list(max_feat = 6000, max_umi = 20000),
    mito_pat = NULL,
    ribo_pat = NULL,
    species = c("human", "mouse"),
    sample_col = "orig.ident",
    outdir = "QC",
    merge_output = TRUE,
    verbose = TRUE
) {
  library(Seurat)
  library(dplyr)
  library(Matrix)
  
  method <- match.arg(method)
  species <- match.arg(species)
  
  # Default gene patterns
  if (is.null(mito_pat)) mito_pat <- if (species == "human") "^MT-" else "^mt-"
  if (is.null(ribo_pat)) ribo_pat <- if (species == "human") "^RP[LS]" else "^Rp[ls]"
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Input validation
  if (inherits(x, "Seurat")) {
    obj_list <- list(x)
    single <- TRUE
  } else if (is.list(x) && all(sapply(x, inherits, "Seurat"))) {
    obj_list <- x
    single <- FALSE
  } else stop("Input must be a Seurat object or list of Seurat objects.")
  
  qc_sum <- data.frame()
  dbl_sum <- data.frame()
  
  if (verbose) cat("\n>>> Starting quality control for", length(obj_list), "sample(s)...\n")
  
  for (i in seq_along(obj_list)) {
    obj <- obj_list[[i]]
    sample <- if (sample_col %in% colnames(obj@meta.data))
      unique(obj@meta.data[[sample_col]]) else paste0("Sample_", i)
    
    pre_cells <- ncol(obj)
    if (verbose) cat("\n────────────────────────────\n Sample:", sample, "| Cells:", pre_cells, "\n")
    
    # Compute QC metrics
    obj$percent_mito <- PercentageFeatureSet(obj, pattern = mito_pat)
    if (calc_ribo) obj$percent_ribo <- PercentageFeatureSet(obj, pattern = ribo_pat)
    if (calc_drop) {
      counts <- GetAssayData(obj, slot = "counts")
      obj$dropout <- Matrix::colSums(counts == 0) / nrow(counts)
    }
    
    # Compute thresholds
    if (method == "MAD") {
      max_feat <- median(obj$nFeature_RNA) + mad_n * mad(obj$nFeature_RNA)
      max_umi  <- median(obj$nCount_RNA) + mad_n * mad(obj$nCount_RNA)
    } else if (method == "fixed") {
      max_feat <- fixed_thr$max_feat
      max_umi  <- fixed_thr$max_umi
    } else {
      max_feat <- Inf; max_umi <- Inf
    }
    
    # Apply permissive filters
    filt <- subset(
      obj,
      nFeature_RNA > min_feat &
        nFeature_RNA < max_feat &
        nCount_RNA > min_umi &
        percent_mito < max_mito &
        (if (calc_ribo) percent_ribo < max_ribo else TRUE) &
        (if (calc_drop) dropout < max_drop else TRUE)
      # No max UMI cutoff — preserves high-depth nuclei
    )
    
    post_cells <- ncol(filt)
    pct_rm <- round((1 - post_cells / pre_cells) * 100, 2)
    if (verbose) cat(">>> Retained:", post_cells, "cells (", pct_rm, "% removed)\n")
    
    # Optional doublet filtering
    if (rm_dbl) {
      if (verbose) cat(">>> Running scDblFinder...\n")
      sce <- scDblFinder::scDblFinder(as.SingleCellExperiment(filt))
      filt$scDblFinder <- sce$scDblFinder.class
      filt <- subset(filt, scDblFinder == "singlet")
      dbl_sum <- rbind(dbl_sum, data.frame(Sample = sample, Doublets = sum(sce$scDblFinder.class == "doublet")))
    }
    
    # Summarize
    qc_sum <- rbind(qc_sum, data.frame(
      Sample = sample,
      Pre_Cells = pre_cells,
      Post_Cells = post_cells,
      Removed_Pct = pct_rm,
      min_feat, max_feat, min_umi,
      max_mito,
      Method = method, MADs = mad_n
    ))
    
    obj_list[[i]] <- filt
  }
  
  # Save summary
  write.csv(qc_sum, file.path(outdir, "QC_summary.csv"), row.names = FALSE)
  if (rm_dbl && nrow(dbl_sum) > 0)
    write.csv(dbl_sum, file.path(outdir, "Doublet_summary.csv"), row.names = FALSE)
  
  # Return object
  if (merge_output && length(obj_list) > 1) {
    merged <- merge(x = obj_list[[1]], y = obj_list[-1], merge.data = TRUE)
    if (verbose) cat("\n>>> Merged object returned.\n")
    if (verbose) cat("\n>>> Join leyers.\n")
    #merged <- JoinLayers(merged)
    #merged <- fix_seurat_matrix_names(merged)
    return(merged)
  } else if (single) {
    if (verbose) cat("\n>>> Returning single filtered object.\n")
    #obj_list[[1]] <- fix_seurat_matrix_names(obj_list[[1]])
    return(obj_list[[1]])
  } else {
    if (verbose) cat("\n>>> Returning list of filtered objects.\n")
    return(obj_list)
  }
}



#' Sex and Erythroid Contamination QC for Seurat Objects
#'
#' This function computes per-cell sex scores based on XIST and Y-linked gene expression,
#' flags sex-discordant cells, calculates erythroid contamination using canonical hemoglobin genes,
#' and filters cells that are both sex-discordant and above a specified erythroid threshold. Optional QC plots are produced.
#'
#' @param obj A Seurat object containing single-cell RNA-seq data.
#' @param sample_col Character. Column name in `obj@meta.data` corresponding to sample IDs. Default is "sample".
#' @param y_genes Character vector. List of Y-chromosome genes to use for sex score. Default is c("UTY","RPS4Y1","ZFY","DDX3Y","KDM5D").
#' @param xist_gene Character. Gene name for XIST. Default is "XIST".
#' @param eps Numeric. Small pseudocount added to avoid division by zero. Default is 1e-6.
#' @param female_thresh Numeric. Sex score threshold above which a cell is classified as female-like. Default is 1.
#' @param male_thresh Numeric. Sex score threshold below which a cell is classified as male-like. Default is -1.
#' @param eryth_genes Character vector. Canonical hemoglobin genes for erythroid contamination assessment. Default is c("HBB","HBA1","HBA2","HBE1","HBG1","HBG2","HBM").
#' @param eryth_percentile Numeric. Percentile to define high erythroid signal for filtering (0-1). Default is 0.95.
#' @param make_plots Logical. Whether to produce QC plots. Default is TRUE.
#' @param plot_prefix Character. Prefix for saved QC plot filenames. Default is "qc_".
#'
#' @return A list with the following components:
#' \item{object}{The input Seurat object with added metadata columns: sex_score, sex_call, eryth_sum, inferred_sex, sex_discordant.}
#' \item{filtered}{A Seurat object filtered to remove cells that are both sex-discordant and above the erythroid threshold.}
#' \item{summary}{A per-sample summary table including cell counts, fractions female/male/ambiguous, median sex score, fraction of cells with erythroid expression, and median erythroid sum.}
#' \item{eryth_threshold}{Numeric value corresponding to the threshold for high erythroid signal (computed from `eryth_percentile`).}
#' \item{n_removed}{Number of cells removed by filtering.}
#'
#' @examples
#' \dontrun{
#' qc <- sex_contamination_qc(obj, sample_col = "sample")
#' qc$summary
#' filtered_obj <- qc$filtered
#' }
#'
#' @import Seurat
#' @import dplyr
#' @import ggplot2
#' @export
sex_contamination_qc <- function(
    obj,
    sample_col = "sample",
    y_genes = c("UTY","RPS4Y1","ZFY","DDX3Y","KDM5D"),
    xist_gene = "XIST",
    eps = 1e-6,
    female_thresh = 1,
    male_thresh = -1,
    eryth_genes = c("HBB","HBA1","HBA2","HBE1","HBG1","HBG2","HBM"),
    majority_threshold = 0.5,
    eryth_percentile = 0.95,
    remove_ambiguous = TRUE,
    min_sex_score_diff = NULL,
    out_file = NULL
){
  cat("Extracting expression matrix (slot = data)...")
  Seurat::DefaultAssay(obj) <- "RNA"
  expr <- Seurat::GetAssayData(obj, slot = "data")
  meta_df <- obj@meta.data
  
  # Sex score
  cat("Computing sex scores...")
  
  meta_df$XIST <- if (xist_gene %in% rownames(expr)) expr[xist_gene, ] else 0
  
  y_present <- intersect(y_genes, rownames(expr))
  meta_df$Ysum <- if (length(y_present) > 0)
    Matrix::colSums(expr[y_present, , drop = FALSE])
  else 0
  
  meta_df$Ysum[is.na(meta_df$Ysum)] <- 0
  meta_df$XIST[is.na(meta_df$XIST)] <- 0
  
  meta_df$sex_score <- log2((meta_df$XIST + eps)/(meta_df$Ysum + eps))
  
  meta_df$sex_call <- ifelse(
    meta_df$sex_score > female_thresh, "Female-like",
    ifelse(meta_df$sex_score < male_thresh, "Male-like", "Ambiguous")
  )
  
  # Erythroid score
  cat("Computing erythroid scores...")
  
  eryth_present <- intersect(eryth_genes, rownames(expr))
  meta_df$eryth_sum <- if (length(eryth_present) > 0)
    Matrix::colSums(expr[eryth_present, , drop = FALSE])
  else 0
  
  meta_df$eryth_sum[is.na(meta_df$eryth_sum)] <- 0
  
  eryth_threshold <- as.numeric(
    stats::quantile(meta_df$eryth_sum, eryth_percentile, na.rm = TRUE)
  )
  
  # Summary per sample
  cat("Summarizing per-sample sex composition...")
  
  tbl_summary <- meta_df %>%
    dplyr::group_by(.data[[sample_col]]) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      frac_female_like = mean(sex_call == "Female-like"),
      frac_male_like   = mean(sex_call == "Male-like"),
      frac_ambiguous   = mean(sex_call == "Ambiguous"),
      median_sex_score = stats::median(sex_score),
      median_eryth_sum = stats::median(eryth_sum),
      .groups = "drop"
    )
  
  # Infer sample sex
  inferred <- tbl_summary %>%
    dplyr::mutate(inferred_sex =
                    ifelse(frac_female_like > majority_threshold, "Female", "Male")) %>%
    dplyr::select(.data[[sample_col]], inferred_sex) %>%
    tibble::deframe()
  
  meta_df$inferred_sex <- unname(inferred[meta_df[[sample_col]]])
  
  # Discordance
  meta_df$sex_discordant <-
    (meta_df$inferred_sex == "Female" & meta_df$sex_call == "Male-like") |
    (meta_df$inferred_sex == "Male"   & meta_df$sex_call == "Female-like")
  
  if (!is.null(min_sex_score_diff)) {
    global_median <- stats::median(meta_df$sex_score, na.rm = TRUE)
    meta_df$strong_discordance <- meta_df$sex_discordant &
      abs(meta_df$sex_score - global_median) >= min_sex_score_diff
  } else {
    meta_df$strong_discordance <- meta_df$sex_discordant
  }
  
  # Flag cells for removal
  cat("Flagging cells for removal...")
  
  remove_cells <- meta_df$strong_discordance & (meta_df$eryth_sum > eryth_threshold)
  
  if (remove_ambiguous)
    remove_cells <- remove_cells | meta_df$sex_call == "Ambiguous"
  
  cat(sum(remove_cells), " cells flagged.")
  
  # Write metadat back
  obj <- Seurat::AddMetaData(obj, metadata = meta_df)
  
  # Subset filtered object
  cat("Creating filtered object...")
  
  keep_cells <- rownames(meta_df)[!remove_cells]
  filt_obj <- subset(obj, cells = keep_cells)
  
  # Summary for output
  removal_summary <- data.frame(
    sample = meta_df[[sample_col]],
    removed = remove_cells
  ) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      n_removed = sum(removed),
      frac_removed = mean(removed),
      .groups = "drop"
    )
  
  tbl_summary <- tbl_summary %>%
    dplyr::left_join(removal_summary, by = "sample") %>%
    dplyr::mutate(eryth_threshold = eryth_threshold)
  
  # Output file
  if (!is.null(out_file)) {
    utils::write.table(tbl_summary, file = out_file, sep = "\t",
                       quote = FALSE, row.names = FALSE)
  }
  
  # Fix seurat matrix names
  obj <- fix_seurat_matrix_names(obj)
  filt_obj <- fix_seurat_matrix_names(filt_obj)
  
  # Return output
  return(list(
    object = obj,
    filtered = filt_obj,
    tbl_summary = tbl_summary,
    eryth_threshold = eryth_threshold,
    n_removed = sum(remove_cells)
  ))
}


genefilter <- function(
    obj,
    mode = c("B","A","C"),
    keep.genes = c("MEG3"),
    excl.genes = NULL,
    assay = "RNA",
    normalize = TRUE,
    norm.method = "LogNormalize",
    scale.fac = 10000,
    find.hvg = TRUE,
    n.hvg = 2000,
    scale.data = TRUE,
    verbose = TRUE
){
  stopifnot(inherits(obj, "Seurat"))
  mode <- match.arg(mode)
  DefaultAssay(obj) <- assay
  
  genes <- rownames(obj[[assay]])
  
  # robust meta.features accessor (works for Seurat v5) 
  get_gene_biotype <- function(obj, assay) {
    assay_obj <- obj[[assay]]
    
    # Assay5: meta.features only exists if user added it
    if (!"meta.features" %in% slotNames(assay_obj)) {
      return(NULL)
    }
    mf <- assay_obj@meta.features
    if (is.null(mf) || nrow(mf) == 0) return(NULL)
    
    mf <- as.data.frame(mf)
    
    if (!"gene_biotype" %in% colnames(mf)) return(NULL)
    
    x <- mf$gene_biotype
    names(x) <- rownames(mf)
    return(x)
  }
  
  gb <- get_gene_biotype(obj, assay)
  
  # filtering rules 
  filters <- list(
    biotypes_rm = c(
      "lincRNA","lncRNA","antisense_RNA","processed_transcript",
      "pseudogene","processed_pseudogene","unprocessed_pseudogene",
      "transcribed_unprocessed_pseudogene","transcribed_processed_pseudogene",
      "unitary_pseudogene","Mt_rRNA","Mt_tRNA","rRNA"
    ),
    regex_rm = c(
      "^MT-",
      "^(RPL|RPS)[0-9]+",
      "^HB[ABDGEMQZ]",
      "^LINC[0-9]+",
      "^(AC|AL|AP|RP)[0-9]+",
      "^(CTD|CTC)-",
      "^LOC[0-9]+",
      "^C[0-9]+orf", "-AS[0-9]+$",
      "^FAM[0-9]+","^DAZ[0-9]+$",
      "^XXBAC", "^RPS4Y[12]$",
      "\\.[0-9]+$", "^TTTY[0-9]+$",
      "^ENSG[0-9]+$","^USP9Y$",
      "^XIST$", "^TSIX$",
      "^SRY$", "^ZFY$", "^UTY$"
    )
  )
  
  genes_up <- toupper(genes)
  names(genes_up) <- genes
  
  rm_biotype <- c()
  if (!is.null(gb)) {
    for (bt in filters$biotypes_rm)
      rm_biotype <- c(rm_biotype, names(gb)[gb == bt])
  }
  
  rm_regex <- c()
  for (rx in filters$regex_rm)
    rm_regex <- c(rm_regex, genes[grepl(rx, genes_up)])
  
  must_remove <- unique(c(rm_biotype, rm_regex))
  
  # keep sets 
  if (!is.null(gb)) {
    prot_ok <- names(gb)[gb == "protein_coding"]
  } else {
    prot_ok <- setdiff(genes, must_remove)
  }
  
  if (mode == "A") {
    keep <- prot_ok
  } else if (mode == "C") {
    keep <- setdiff(genes, must_remove)
  } else {
    keep <- union(prot_ok, keep.genes)
  }
  
  if (!is.null(excl.genes))
    keep <- setdiff(keep, excl.genes)
  
  keep <- setdiff(keep, must_remove)
  removed <- setdiff(genes, keep)
  
  if (verbose) {
    message("=== filter_genes (Seurat v5-safe) ===")
    message("Genes:", length(genes))
    message("Kept :", length(keep))
    message("Removed:", length(removed))
  }
  
  # subset + normalization 
  obj_f <- subset(obj, features = keep)
  
  if (normalize)
    obj_f <- NormalizeData(obj_f, normalization.method = norm.method, scale.factor = scale.fac)
  
  if (find.hvg)
    obj_f <- FindVariableFeatures(obj_f, nfeatures = n.hvg)
  
  if (scale.data)
    obj_f <- ScaleData(obj_f, features = VariableFeatures(obj_f))
  
  return(obj_f)
}



#' Filter Outliers in Seurat Clusters Using Local Density
#'
#' @description
#' Removes low-density (outlier) cells within specified clusters of a Seurat object,
#' based on k-nearest neighbor (KNN) distances in a given dimensional reduction (e.g., UMAP or PCA).
#' Optionally visualizes before/after filtering for each cluster and summarizes retained/removed cell counts.
#'
#' @param obj A Seurat object.
#' @param group.by Character string. Metadata column used to define clusters (default: `"seurat_clusters"`).
#' @param reduction Character string. Dimensional reduction to use (e.g., `"umap"`, `"pca"`).
#' @param density.threshold Numeric in (0,1). Quantile of local density defining cutoff (default: `0.95`).
#' @param k Integer. Number of neighbors for KNN-based local density (default: `30`).
#' @param relax Numeric multiplier (>1) to loosen the density threshold (default: `2`).
#' @param clean.all Logical. If `TRUE`, filter all clusters; otherwise specify `target.clusters`.
#' @param target.clusters Vector of cluster identities to clean. Ignored if `clean.all = TRUE`.
#' @param verbose Logical. Print progress and summary (default: `TRUE`).
#' @param color Color for highlighted cells in plots (default: `"red"`).
#' @param pt.size Point size for plotting (default: `1`).
#' @param base.size Base text size for plot titles (default: `10`).
#' @param ncol Integer. Number of columns for the combined plot layout (default: `2`).
#' @param return.plot Logical. Return combined patchwork plot (default: `FALSE`).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{obj}: The filtered Seurat object.
#'   \item \code{summary}: A data.frame summarizing before/after cell counts per cluster.
#'   \item \code{plot}: (optional) Combined before/after UMAP visualization if `return.plot = TRUE`.
#' }
#'
#' @details
#' The function estimates local cell density for each cluster using mean KNN distances.
#' Cells below the density cutoff (low-density regions) are considered outliers and removed.
#' Clusters not listed in `target.clusters` remain untouched.
#'
#' @examples
#' \dontrun{
#' res <- filter_outliers(
#'   obj = obj,
#'   group.by = "seurat_clusters",
#'   target.clusters = c(0,1,2,3),
#'   reduction = "umap",
#'   density.threshold = 0.95,
#'   k = 30,
#'   relax = 2,
#'   verbose = TRUE,
#'   return.plot = TRUE
#' )
#' res$summary
#' print(res$plot)
#' }
#'
#' @import Seurat
#' @import FNN
#' @import patchwork
#' @export
filter_outliers <- function(obj,
                            group.by = "seurat_clusters",
                            reduction = "umap",
                            density.threshold = 0.95,
                            k = 30,
                            relax = 2,
                            clean.all = FALSE,
                            target.clusters = NULL,
                            verbose = TRUE,
                            color = "red",
                            pt.size = 1,
                            base.size = 10,
                            ncol = 2,
                            return.plot = FALSE) {
  
  require(Seurat)
  require(FNN)
  require(patchwork)
  
  #Checks
  if (!group.by %in% colnames(obj@meta.data))
    stop("Column ", group.by, " not found in metadata.")
  if (!reduction %in% Reductions(obj))
    stop("Reduction ", reduction, " not found in object.")
  
  obj <- SetIdent(obj, value = group.by)
  uniq.clust <- unique(na.omit(as.character(obj@meta.data[[group.by]])))
  
  if (clean.all) {
    clusts <- uniq.clust
  } else if (!is.null(target.clusters)) {
    bad <- setdiff(target.clusters, uniq.clust)
    if (length(bad) > 0) stop("Invalid clusters: ", paste(bad, collapse = ", "))
    clusts <- target.clusters
  } else {
    stop("Specify either clean.all = TRUE or target.clusters.")
  }
  
  if (verbose) cat("Processing clusters:", paste(clusts, collapse = ", "), "\n")
  
  #Loop
  all.keep <- c()
  plots <- list()
  summ <- data.frame()
  emb <- Embeddings(obj, reduction)
  
  for (cl in clusts) {
    cells <- rownames(obj@meta.data)[obj@meta.data[[group.by]] == cl]
    if (length(cells) < 2) next
    
    coord <- emb[cells, , drop = FALSE]
    kk <- min(k, length(cells) - 1)
    dens <- rowMeans(knn.dist(coord, k = kk), na.rm = TRUE)
    obj@meta.data[cells, "dens"] <- dens
    
    thr <- quantile(dens, density.threshold, na.rm = TRUE)
    keep <- cells[dens <= thr * relax]
    all.keep <- c(all.keep, keep)
    
    summ <- rbind(summ, data.frame(
      cluster = cl,
      before = length(cells),
      after = length(keep),
      removed = length(cells) - length(keep),
      pct.rm = round((1 - length(keep)/length(cells)) * 100, 2)
    ))
    
    p1 <- DimPlot(obj, reduction = reduction, group.by = group.by,
                  cells.highlight = cells, cols.highlight = color,
                  sizes.highlight = pt.size) +
      NoAxes() + NoLegend() +
      ggtitle(paste("Cluster", cl, "(n=", length(cells), ")")) +
      theme(plot.title = element_text(hjust = 0.5, size = base.size))
    
    p2 <- DimPlot(obj, reduction = reduction, group.by = group.by,
                  cells.highlight = keep, cols.highlight = color,
                  sizes.highlight = pt.size) +
      NoAxes() + NoLegend() +
      ggtitle(paste("Filtered", cl, "(n=", length(keep),
                    ")\nCutoff=", round(thr, 2), " K=", kk)) +
      theme(plot.title = element_text(hjust = 0.5, size = base.size))
    
    plots[[as.character(cl)]] <- p1 + p2 + patchwork::plot_layout(ncol = 2)
  }
  
  #Keep filtered + untouched clusters
  all.target <- rownames(obj@meta.data)[obj@meta.data[[group.by]] %in% clusts]
  other.cells <- setdiff(colnames(obj), all.target)
  keep.cells <- unique(c(all.keep, other.cells))
  filt <- subset(obj, cells = keep.cells)
  
  if (verbose) {
    cat("\n=== FILTER SUMMARY ===\n")
    print(summ, row.names = FALSE)
    cat("Total before:", ncol(obj), "\n")
    cat("Total after:", ncol(filt), "\n")
    cat("Removed:", ncol(obj) - ncol(filt), "\n")
  }
  
  if (return.plot) {
    combo <- wrap_plots(plots, ncol = ncol)
    return(list(obj = filt, summary = summ, plot = combo))
  } else {
    return(list(obj = filt, summary = summ))
  }
}

#' Compute cluster statistics from Seurat object (multiple group.by)
#'
#' @param object Seurat object
#' @param ident Metadata column for cluster identities (default: "seurat_clusters")
#' @param group.by One or more metadata columns to group by (default: "orig.ident")
#' @param rm.sum Logical; remove total row (default: FALSE)
#' @param rm.pct Logical; hide percentage columns (default: FALSE)
#' @return Data frame of cluster counts and percentages
#' @export
cluster_stats <- function(object, ident = "seurat_clusters", 
                          group.by = "orig.ident", rm.sum = FALSE, rm.pct = FALSE) {
  
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.")
  if (!missing(ident)) object <- SetIdent(object, value = ident)
  
  group.by <- as.character(group.by)
  missing.cols <- setdiff(group.by, colnames(object@meta.data))
  if(length(missing.cols)) stop(sprintf("'%s' not in meta.data.", paste(missing.cols, collapse = ", ")))
  
  # Base cluster totals
  tbl <- table(object@active.ident)
  stats <- data.frame(cluster = names(tbl), count = as.numeric(tbl), pct = as.numeric(prop.table(tbl)*100))
  
  # Function to compute counts and percent per group.by
  compute_group <- function(g) {
    tab <- table(object@active.ident, object@meta.data[[g]])
    df <- as.data.frame.matrix(tab) %>% tibble::rownames_to_column("cluster")
    pct <- prop.table(tab, margin = 2) * 100
    pct.df <- as.data.frame.matrix(pct) %>% tibble::rownames_to_column("cluster")
    colnames(pct.df)[-1] <- paste0(colnames(pct.df)[-1], ".", g, ".pct")
    cbind(df[,-1], pct.df[,-1])
  }
  
  if(length(group.by) > 0) {
    stats <- cbind(stats, purrr::map_dfc(group.by, compute_group))
  }
  
  stats <- janitor::adorn_totals(stats, "row")
  if(rm.sum) stats <- stats[1:(nrow(stats)-1), ]
  if(rm.pct) stats <- stats[, !grepl("\\.pct$", colnames(stats))]
  
  stats
}


#' Automatically Select the Optimal Number of PCs
#'
#' This function computes an optimal number of principal components (PCs)
#' to use for downstream analysis based on variance explained. It supports
#' any dimensional reduction stored in a Seurat object (e.g., PCA, Harmony).
#'
#' @param object A Seurat object.
#' @param reduction.name Name of the reduction slot to inspect
#'   (default: \code{"pca"}). Must exist in \code{object[[reduction.name]]}.
#'
#' @details
#' The selection heuristic uses two criteria:
#'
#' \itemize{
#'   \item \strong{co1}: First PC where cumulative variance exceeds 90\% AND
#'         the PC-specific variance drops below 5\%.
#'
#'   \item \strong{co2}: Last PC where the drop in variance from the previous PC
#'         is greater than 0.1\%.
#' }
#'
#' The function returns the larger of the two values, ensuring sufficient
#' dimensionality for clustering and UMAP.
#'
#' @return An integer giving the recommended number of dimensions.
#'
#' @examples
#' \dontrun{
#' pcs <- get_pcs(object, reduction.name = "pca")
#' }
#'
#' @export
#' 
get_pcs <- function(object, reduction.name="pca") {
  # Check seurat object
  if(class(object) != "Seurat") {
    message("WARNING: this rds file does not contain a Seurat object! STOP RUNNING THIS SCRIPT")
    message("Check the data type by running:")
    message("class(obj)")
    stop()
  }
  # Determine percent of variation associated with each PC
  pct <- object[[reduction.name]]@stdev / sum(object[[reduction.name]]@stdev)*100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % 
  # variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  co1
  # Determine the difference between variation of PC and subsequent PC and
  # selecting last point where change of % of variation is more than 0.1%.
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # Minimum of the two calculation
  #pcs <- min(co1, co2)
  c(co1, co2)
  pcs <- max(co1, co2)
  #pcs <- max(co1, co2)
  return(pcs)
}


#' Compute Unique Gene Expression Statistics
#'
#' This function evaluates a set of genes in a Seurat object and returns:
#' \enumerate{
#'   \item The number of selected genes expressed by each cell.
#'   \item The number of cells expressing each gene.
#' }
#'
#' @param object A Seurat object containing an RNA assay.
#' @param gene.list Character vector of gene names to evaluate.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{per.cell}}{A data frame with columns \code{Cell} and \code{Unique}.}
#'   \item{\code{per.gene}}{A data frame with columns \code{gene} and \code{cells.expressing}.}
#' }
#'
#' @export
get_unique_gene_table <- function(object, gene.list, assay = "RNA", layer = "counts") {
  if (!inherits(object, "Seurat"))
    stop("`object` must be a Seurat object.")
  
  if (!assay %in% names(object@assays))
    stop("Assay not found in object.")
  
  assay_obj <- object[[assay]]
  
  # Access the matrix via @layers (Seurat v5)
  if (!layer %in% names(assay_obj@layers))
    stop(paste0("Layer '", layer, "' not found in the assay."))
  
  mat <- assay_obj@layers[[layer]]
  
  # Ensure dimnames are set
  if (is.null(rownames(mat))) rownames(mat) <- rownames(assay_obj)
  if (is.null(colnames(mat))) colnames(mat) <- colnames(object)
  
  genes.present <- intersect(gene.list, rownames(mat))
  if (length(genes.present) == 0)
    stop("None of the provided genes exist in the object.")
  
  mat <- mat[genes.present, , drop = FALSE]
  mat <- mat[Matrix::rowSums(mat > 0) > 0, , drop = FALSE]
  
  per.cell <- data.frame(
    Cell = colnames(mat),
    Unique = Matrix::colSums(mat > 0)
  )
  
  per.gene <- data.frame(
    gene = rownames(mat),
    cells.expressing = Matrix::rowSums(mat > 0)
  )
  
  per.gene <- per.gene[order(per.gene$cells.expressing), ]
  
  list(per.cell = per.cell, per.gene = per.gene)
}


run_decontx <- function(
    input,
    assay_name = "RNA",
    verbose = TRUE,
    ...
) {
  library(Seurat)
  library(SingleCellExperiment)
  library(scater)
  library(celda)
  library(Matrix)
  
  # Load input
  if (inherits(input, "Seurat")) {
    if (verbose) cat("Input is a Seurat object\n")
    
    sce <- as.SingleCellExperiment(input)
    
    # Extract original UMAP from Seurat
    if ("umap" %in% names(input@reductions)) {
      reducedDims(sce)$UMAP_before <- Embeddings(input, "umap")
      if (verbose) cat("UMAP_before extracted from Seurat\n")
    } else {
      reducedDims(sce)$UMAP_before <- NULL
      if (verbose) cat("No UMAP found in Seurat\n")
    }
    
  } else if (inherits(input, "SingleCellExperiment")) {
    
    sce <- input
    if (verbose) cat("Input is SingleCellExperiment\n")
    
    # Preserve existing UMAP
    if ("UMAP" %in% names(reducedDims(sce))) {
      reducedDims(sce)$UMAP_before <- reducedDims(sce)$UMAP
      if (verbose) cat("UMAP_before extracted from SCE\n")
    } else {
      reducedDims(sce)$UMAP_before <- NULL
      if (verbose) cat("No UMAP present in SCE\n")
    }
    
  } else if (is.matrix(input) || inherits(input, "dgCMatrix")) {
    
    if (verbose) cat("Input is counts matrix\n")
    sce <- SingleCellExperiment(assays=list(counts=input))
    reducedDims(sce)$UMAP_before <- NULL
    
  } else {
    stop("Input must be Seurat, SCE, or matrix.")
  }
  
  # Run decontx
  if (verbose) cat("Running DecontX...\n")
  sce <- decontX(sce, ...)
  
  cleaned <- decontXcounts(sce)
  cleaned <- round(cleaned)
  
  # Compute umap after
  if (verbose) cat("Computing UMAP_after from DecontX-cleaned counts...\n")
  
  sce <- runUMAP(sce, exprs_values = "decontXcounts")
  reducedDims(sce)$UMAP_after <- reducedDims(sce)$UMAP
  reducedDims(sce)$UMAP <- NULL    # avoid confusion
  
  # Return both
  return(list(
    sce = sce,
    UMAP_before = reducedDims(sce)$UMAP_before,
    UMAP_after  = reducedDims(sce)$UMAP_after,
    cleaned_counts = cleaned
  ))
}


save_heatmap <- function(
    ht,
    filename,
    width = 8,
    height = 6,
    formats = c("pdf"),
    dpi = 600
) {
  
  stopifnot(!missing(ht))
  stopifnot(is.character(filename))
  stopifnot(is.numeric(width), is.numeric(height))
  
  for (fmt in formats) {
    
    fmt <- tolower(fmt)
    file_out <- paste0(filename, ".", fmt)
    
    if (fmt == "pdf") {
      
      pdf(
        file_out,
        width = width,
        height = height,
        useDingbats = FALSE
      )
      
    } else if (fmt == "png") {
      
      png(
        file_out,
        width = width,
        height = height,
        units = "in",
        res = dpi
      )
      
    } else if (fmt == "tiff") {
      
      tiff(
        file_out,
        width = width,
        height = height,
        units = "in",
        res = dpi,
        compression = "lzw"
      )
      
    } else if (fmt %in% c("jpg", "jpeg")) {
      
      jpeg(
        file_out,
        width = width,
        height = height,
        units = "in",
        res = dpi,
        quality = 100
      )
      
    } else {
      stop("Unsupported format: ", fmt)
    }
    
    # Draw AFTER opening device
    ComplexHeatmap::draw(ht)
    
    dev.off()
  }
  
  invisible(TRUE)
}



#' Fix matrix names in Seurat v5 objects (robust)
#'
#' Ensures that all matrices in a Seurat RNA assay (`layers`, `data`, `scale.data`)
#' have proper row and column names. Skips slots that do not exist or are empty.
#'
#' @param sobj A Seurat object
#' @return A Seurat object with all matrix names fixed
fix_seurat_matrix_names <- function(sobj) {
  assay <- sobj[["RNA"]]
  
  for (slot.name in c("layers", "data", "scale.data")) {
    # Handle layers specially
    if (slot.name == "layers" && !is.null(assay@layers$counts)) {
      mat <- assay@layers$counts
      if (!is.null(mat) && length(dim(mat)) == 2) {
        if (is.null(rownames(mat))) rownames(mat) <- rownames(assay)
        if (is.null(colnames(mat))) colnames(mat) <- colnames(sobj)
        assay@layers$counts <- mat
      }
    } 
    # Other slots
    else if (slot.name %in% slotNames(assay) && !is.null(slot(assay, slot.name))) {
      mat <- slot(assay, slot.name)
      if (!is.null(mat) && length(dim(mat)) == 2) {
        if (is.null(rownames(mat))) rownames(mat) <- rownames(assay)
        if (is.null(colnames(mat))) colnames(mat) <- colnames(sobj)
        slot(assay, slot.name) <- mat
      }
    }
  }
  
  sobj[["RNA"]] <- assay
  return(sobj)
}


savefig <- function(
    filename,
    fig = NULL,
    layout = c("single", "double"),
    type = c("ggplot", "heatmap", "baseplot"),
    format = c("pdf", "png", "tiff"),
    width = NULL,
    height = NULL,
    dpi = 600
){
  layout <- match.arg(layout)
  type <- match.arg(type)
  formats <- format   # allow vector of formats
  
  # Nature Communications sizing defaults
  default_width <- ifelse(layout == "single", 3.5, 7.2)
  default_height <- ifelse(type %in% c("heatmap", "baseplot"), 8, 5)
  
  width <- width %||% default_width
  height <- height %||% default_height
  
  for(fmt in formats) {
    
    out_file <- filename
    if(!grepl(paste0("\\.", fmt, "$"), out_file))
      out_file <- paste0(out_file, ".", fmt)
    
    cat("Saving ", out_file,
        "\n layout=", layout,
        "\n type=", type,
        "\n size=", width, "×", height, " inches @ DPI=", dpi,
        "\n format=", fmt, "\n\n")
    
    # PDF
    if(fmt == "pdf") {
      pdf(out_file, width = width, height = height, useDingbats = FALSE)
      if(type == "ggplot") print(fig)
      if(type %in% c("heatmap", "baseplot")) eval(fig)
      dev.off()
      next
    }
    
    # PNG / TIFF
    if(fmt %in% c("png", "tiff")) {
      if(type == "ggplot") {
        ggsave(out_file, plot = fig,
               width = width, height = height,
               units = "in", dpi = dpi)
        
      } else if(type %in% c("heatmap", "baseplot")) {
        
        if(fmt == "png")
          png(out_file, width = width, height = height, units = "in", res = dpi)
        
        if(fmt == "tiff")
          tiff(out_file, width = width, height = height, units = "in",
               res = dpi, compression = "lzw")
        
        eval(fig)
        dev.off()
      }
    }
  }
  
  invisible(filename)
}



#' Colorize Values Using Various Color Palettes
#'
#' This function assigns colors to unique values in a vector using a variety of distinct and continuous color palettes.  
#' It supports both discrete and continuous color schemes, including RColorBrewer, Viridis, and custom palettes.
#'
#' @param x A vector of values to be colorized.
#' @param theme A string specifying the color theme. Available themes:
#'   * **Distinct palettes**: `"alphabet"`, `"tube"`, `"tol"`, `"glasbey"`, `"tableau"`, `"bright"`, `"pastel"`, `"cList"``
#'   * **RColorBrewer palettes**:
#'     - *Sequential*: `"Blues"`, `"BuGn"`, `"BuPu"`, `"GnBu"`, `"Greens"`, `"Greys"`, `"Oranges"`, `"OrRd"`, `"PuBu"`,
#'                     `"PuBuGn"`, `"PuRd"`, `"Purples"`, `"RdPu"`, `"Reds"`, `"YlGn"`, `"YlGnBu"`, `"YlOrBr"`, `"YlOrRd"`
#'     - *Qualitative*: `"Accent"`, `"Dark2"`, `"Paired"`, `"Pastel1"`, `"Pastel2"`, `"Set1"`, `"Set2"`, `"Set3"`
#'     - *Divergent*: `"BrBG"`, `"PiYG"`, `"PRGn"`, `"PuOr"`, `"RdBu"`, `"RdGy"`, `"RdYlBu"`, `"RdYlGn"`, `"Spectral"`
#'   * **Viridis palettes**: `"viridis"`, `"magma"`, `"plasma"`, `"cividis"`, `"inferno"`
#'   * **Random colors**: `"random"`
#' @param bgval A value in `x` that should be assigned a background color (`"#CCCCCC"` by default). Default is `NULL` (no special background value).
#' @param n_colors Number of colors to generate for continuous palettes. Default is `100`.
#'
#' @return A named vector of hex color codes, where names correspond to unique values in `x`.
#'
#' @examples
#' # Apply Viridis color palette
#' colorize(1:10, theme = "viridis")
#'
#' # Apply RColorBrewer "Blues" palette
#' colorize(letters[1:10], theme = "Blues")
#'
#' # Use a distinct color palette (Glasbey)
#' colorize(1:15, theme = "glasbey")
#'
#' @import viridis
#' @import RColorBrewer
#' @export
colorize <- function(x, theme = "alphabet", bgval = NULL, n_colors = 100, reverse = FALSE) {
  
  require(viridis)
  require(RColorBrewer)
  
  # Define distinct color palettes
  color_alphabet <- matrix(c(
    240, 163, 255, 0, 117, 220, 153, 63, 0, 76, 0, 92, 0, 92, 49, 43, 206, 72, 255, 204, 153,
    128, 128, 128, 148, 255, 181, 143, 124, 0, 157, 204, 0, 194, 0, 136, 0, 51, 128, 255, 164, 5,
    255, 168, 187, 66, 102, 0, 255, 0, 16, 94, 241, 242, 0, 153, 143, 224, 255, 102, 116, 10, 255,
    153, 0, 0, 255, 255, 128, 255, 255, 0, 255, 80, 5
  ), ncol = 3, byrow = TRUE) / 256
  
  tube_colors <- c(
    "#B36305", "#E32017", "#FFD300", "#00782A", "#F3A9BB",
    "#A0A5A9", "#9B0056", "#000000", "#003688", "#0098D4",
    "#95CDBA", "#00A4A7", "#EE7C0E", "#84B817", "#E21836",
    "#7156A5"
  )
  
  tol_colors <- c(
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
    "#DDCC77", "#CC6677", "#882255", "#AA4499", "#DDDDDD",
    "#E69F00", "#56B4E9"
  )
  
  glasbey_colors <- c(
    "#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF",
    "#00FFFF", "#800000", "#808000", "#008000", "#800080",
    "#808080", "#C0C0C0", "#008080", "#000080", "#FFA500"
  )
  
  tableau_colors <- c(
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#BCBD22", "#17BECF",
    "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
    "#C49C94", "#F7B6D2", "#DBDB8D", "#9EDAE5"
  )
  
  bright_colors <- c(
    "#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",
    "#911EB4", "#42D4F4", "#F032E6", "#BFEF45", "#FABEBE"
  )
  
  pastel_colors <- c(
    "#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF",
    "#E6C0E9", "#D9C3A1", "#C4E5F5", "#F6D7A7", "#D1E8E2"
  )
  
  cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84",
                 "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"),
               c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF",
                 "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)],
               c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C",
                 "#2C728E","#3B528B","#472D7B","#440154"))
  
  # cList = c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000", "#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
  #                "#FEE090","#FDAE61","#F46D43","#D73027", "#FDE725","#AADC32","#5DC863","#27AD81","#21908C", "#2C728E","#3B528B","#472D7B","#440154")
  
  viridis_colors <- viridis(n_colors, option = "viridis")
  magma_colors <- viridis(n_colors, option = "magma")
  plasma_colors <- viridis(n_colors, option = "plasma")
  inferno_colors <- viridis(n_colors, option = "inferno")
  cividis_colors <- viridis(n_colors, option = "cividis")
  
  brewer_palettes <- c('Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd', 'PuBu', 
                       'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 
                       "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3", 
                       "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
  
  set.seed(42)
  random_colors <- rgb(runif(n_colors), runif(n_colors), runif(n_colors))
  
  # Choose color theme
  if (theme == "alphabet") {
    colors <- apply(rbind(color_alphabet, 1 - (1 - color_alphabet) / 2, color_alphabet / 2), 1, 
                    function(row) rgb(row[1], row[2], row[3]))
  } else if (theme == "tube") {
    colors <- tube_colors
  } else if (theme == "cList") {
    colors <- cList
  } else if (theme == "tol") {
    colors <- tol_colors
  } else if (theme == "glasbey") {
    colors <- glasbey_colors
  } else if (theme == "tableau") {
    colors <- tableau_colors
  } else if (theme == "bright") {
    colors <- bright_colors
  } else if (theme == "pastel") {
    colors <- pastel_colors
  } else if (theme == "random") {
    colors <- random_colors
  } else if (theme %in% brewer_palettes) {
    colors <- colorRampPalette(brewer.pal(8, theme))(n_colors) 
  } else if (theme == "viridis") {
    colors <- viridis_colors
  } else if (theme == "magma") {
    colors <- magma_colors
  } else if (theme == "plasma") {
    colors <- plasma_colors
  } else if (theme == "cividis") {
    colors <- cividis_colors
  } else if (theme == "inferno") {
    colors <- inferno_colors
  } else {
    stop("Invalid theme.")
  }
  
  # Reverse colors if needed
  if (reverse) {
    colors <- rev(colors)
  }
  
  unique_vals <- unique(x)
  encoded_vals <- match(x, unique_vals) - 1
  colors <- colors[(encoded_vals %% length(colors)) + 1]
  
  if (!is.null(bgval)) {
    colors[x == bgval] <- "#CCCCCC"
  }
  
  names(colors) <- unique_vals
  return(colors)
}


#' Safe subset of Seurat object
#'
#' Wrapper around Seurat::subset that fixes matrix names automatically.
#'
#' @param sobj Seurat object
#' @param ... Arguments passed to Seurat::subset
#' @return Subsetted Seurat object
get_subset <- function(sobj, group.by = NULL, invert = FALSE,...) {
  if (!is.null(group.by)) sobj <- SetIdent(sobj, value = group.by)
  sobj2 <- subset(sobj, invert = invert,...)
  sobj2 <- fix_seurat_matrix_names(sobj2)
  return(sobj2)
}

#' Safe merge of Seurat objects
#'
#' Wrapper around Seurat::merge that fixes matrix names automatically.
#'
#' @param x Seurat object
#' @param y Seurat object
#' @param ... Additional arguments to Seurat::merge
#' @return Merged Seurat object
get_merge <- function(obj.list) {
  if (length(obj.list) == 0) stop("No Seurat objects provided")
  merged <- Reduce(merge, obj.list)
  merged <- JoinLayers(merged)
  assay <- merged[["RNA"]]
  assay@layers$counts <- GetAssayData(merged, slot = "counts")
  merged[["RNA"]] <- assay
  merged <- fix_seurat_matrix_names(merged)
  return(merged)
}

#' Load multiple 10X/Seurat samples with optional merge and save
#'
#' @param sample.dirs Character vector of sample folder names
#' @param parent.dir Parent directory containing sample folders
#' @param use.filtered Logical, whether to use filtered matrices (default TRUE)
#' @param save.dir Directory to save individual Seurat objects (default NULL)
#' @param merge Logical, whether to merge all Seurat objects (default FALSE)
#' @return List of Seurat objects or a single merged Seurat object
#' @examples
#' objs <- load_seurat(c("FN_S1256","FN_S3478"), parent.dir = "data/processed")
load_seurat <- function(sample.dirs,
                        parent.dir,
                        use.filtered = TRUE,
                        save.dir = NULL,
                        merge = FALSE,
                        verbose = TRUE) {
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
  })
  
  load_one <- function(sample.name) {
    
    log <- function(...) if (verbose) cat(..., "\n")
    
    log("\n[INFO] Loading sample:", sample.name)
    
    samp.path <- file.path(parent.dir, sample.name)
    if (!dir.exists(samp.path)) {
      stop("[ERROR] Sample folder not found: ", samp.path)
    }
    
    # Patterns
    if (use.filtered) {
      h5.pattern  <- "filtered_feature_bc_matrix\\.h5$"
      dir.pattern <- "filtered_feature_bc_matrix$"
    } else {
      h5.pattern  <- "raw_feature_bc_matrix\\.h5$"
      dir.pattern <- "raw_feature_bc_matrix$"
    }
    
    # Find H5 files 
    h5.files <- list.files(
      samp.path,
      pattern = h5.pattern,
      recursive = TRUE,
      full.names = TRUE
    )
    
    # Find matrix directories
    dirs <- list.dirs(
      samp.path,
      recursive = TRUE,
      full.names = TRUE
    )
    dirs <- dirs[grepl(dir.pattern, dirs)]
    
    log("  Found", length(h5.files), "H5 file(s)")
    log("  Found", length(dirs), "matrix folder(s)")
    
    if (length(h5.files) == 0 && length(dirs) == 0) {
      stop("[ERROR] No valid 10X input found for sample: ", sample.name)
    }
    
    # Read counts 
    counts <- if (length(h5.files) > 0) {
      log("  → Reading H5:", basename(h5.files[1]))
      Read10X_h5(h5.files[1])
    } else {
      log("  → Reading MTX folder:", dirs[1])
      Read10X(dirs[1])
    }
    
    # Create Seurat object 
    log("  Creating Seurat object...")
    sobj <- CreateSeuratObject(
      counts = counts,
      project = sample.name
    )
    
    # Optional matrix-name fix
    if (exists("fix_seurat_matrix_names")) {
      log("  Fixing matrix names...")
      sobj <- fix_seurat_matrix_names(sobj)
    }
    
    log("  Done:",
        ncol(sobj), "cells x",
        nrow(sobj), "features")
    
    return(sobj)
  }
  
  cat("[INFO] Loading samples...\n")
  objs <- setNames(lapply(sample.dirs, load_one), sample.dirs)
  cat("[INFO] All samples loaded\n")
  
  # Save individual objects 
  if (!is.null(save.dir)) {
    dir.create(save.dir, recursive = TRUE, showWarnings = FALSE)
    cat("[INFO] Saving Seurat objects to:", save.dir, "\n")
    for (n in names(objs)) {
      saveRDS(objs[[n]], file.path(save.dir, paste0(n, ".rds")))
    }
  }
  
  # Merge if requested 
  if (merge && length(objs) > 1) {
    cat("[INFO] Merging objects...\n")
    merged <- Reduce(function(x, y) merge(x, y), objs)
    
    if (exists("fix_seurat_matrix_names")) {
      merged <- fix_seurat_matrix_names(merged)
    }
    
    cat("[INFO] Merged object:",
        ncol(merged), "cells x",
        nrow(merged), "features\n")
    
    return(merged)
  }
  
  return(objs)
}


#' Split a Seurat object by metadata and fix matrix names
#'
#' @param object A Seurat object
#' @param split.by Metadata column name to split by
#' @return A named list of Seurat objects with @layers$counts fixed
#' @examples
#' seurat.list <- split_objectects(merged, split.by = "BulkSample")
split_objectects <- function(object, split.by) {
  require(Seurat)
  
  objs <- Seurat::SplitObject(object, split.by = split.by)
  
  objs <- lapply(objs, function(sobj) {
    # Ensure @layers$counts exist
    assay <- sobj[["RNA"]]
    assay@layers$counts <- Seurat::GetAssayData(sobj, slot = "counts")
    sobj[["RNA"]] <- assay
    
    # Fix row/colnames
    sobj <- fix_seurat_matrix_names(sobj)
    return(sobj)
  })
  
  return(objs)
}


subset_seurat <- function(object,
                          group.by = NULL,
                          use.harmony = TRUE,
                          harmony.group = "orig.ident",
                          ndims = 50,
                          resolution = 0.5,
                          min.dist = 0.3,
                          spread = 1,
                          preprocess = TRUE,
                          ...) {
  
  stopifnot(inherits(object, "Seurat"))
  
  # Set identities
  if (!is.null(group.by)) {
    object <- Seurat::SetIdent(object, value = object[[group.by]][,1])
  }
  
  # Subset
  result <- tryCatch(subset(object, ...), error = function(e) {
    stop("Subsetting failed: ", e$message)
  })
  
  if (preprocess) {
    # Normalize + variable features
    result <- Seurat::NormalizeData(result)
    result <- Seurat::FindVariableFeatures(result)
    
    if (length(unique(Seurat::Idents(result))) > 1) {
      result <- Seurat::ScaleData(result)
    }
    
    # PCA / Harmony
    if (use.harmony && length(unique(result@meta.data[[harmony.group]])) > 1) {
      result <- harmony::RunHarmony(result, group.by.vars = harmony.group)
      reduction <- "harmony"
      n_reduc_dims <- ncol(result@reductions[[reduction]]@cell.embeddings)
    } else {
      result <- Seurat::RunPCA(result, npcs = min(ndims, 50))
      reduction <- "pca"
      n_reduc_dims <- ncol(result@reductions[[reduction]]@cell.embeddings)
    }
    
    dims_use <- seq_len(min(ndims, n_reduc_dims))
    
    # Neighbors + Clusters + UMAP
    graph_name <- "RNA_nn"
    result <- Seurat::FindNeighbors(result, dims = dims_use, reduction = reduction, graph.name = graph_name)
    result <- Seurat::FindClusters(result, resolution = resolution, graph.name = graph_name)
    result <- Seurat::RunUMAP(result, dims = dims_use, reduction = reduction, min.dist = min.dist, spread = spread)
  }
  
  return(result)
}


#' Select cells from a Seurat object based on gene expression, clusters, and metadata
#'
#' @description
#' `cellpick()` selects cells from a Seurat object using positive and/or negative
#' gene expression criteria, optional cluster restrictions, and flexible metadata
#' filtering (categorical or numeric ranges).
#'
#' @param obj A Seurat object.
#' @param pos.genes Character vector of genes that should be expressed.
#' @param neg.genes Character vector of genes that must NOT be expressed (default NULL).
#' @param slot Assay slot to use for expression values (default "data").
#' @param expr.thresh Numeric expr.thresh above which a gene is considered expressed (default 0).
#' @param min.genes Minimum number of `pos.genes` that must be expressed
#'   (default 1; use length(pos.genes) to require all).
#' @param clusters Optional vector of cluster IDs to retain (default NULL).
#' @param cluster.col Metadata column containing cluster identities
#'   (default NULL; uses `Idents(obj)`).
#' @param meta.col Metadata column to filter on (default NULL).
#' @param meta.vals Values or numeric range (length 2) used to filter `meta.col`.
#' @param return.obj Logical; if TRUE, return a Seurat object subset
#'   instead of cell names (default FALSE).
#' @param verbose Logical; print filtering diagnostics (default FALSE).
#'
#' @return
#' Character vector of selected cell names, or a Seurat object if `return.obj = TRUE`.
#'
#' @export
#'
#' @examples
#' # Olfactory HBC selection
#' cells <- cellpick(
#'   obj,
#'   pos.genes = c("TP63", "KRT5", "KRT14"),
#'   neg.genes = c("KRT13"),
#'   min.genes = 2,
#'   cluster.col = "seurat_clusters",
#'   clusters = 16,
#'   meta.col = "group",
#'   meta.vals = c("PCW10", "PCW12"),
#'   verbose = TRUE
#' )
cellpick <- function(
    obj,
    pos.genes,
    neg.genes = NULL,
    slot = "data",
    expr.thresh = 0,
    min.genes = 1,
    clusters = NULL,
    cluster.col = NULL,
    meta.col = NULL,
    meta.vals = NULL,
    return.obj = FALSE,
    verbose = FALSE
) {
  
  stopifnot(inherits(obj, "Seurat"))
  stopifnot(length(pos.genes) >= 1)
  
  # Expression matrix
  mat <- Seurat::GetAssayData(obj, slot = slot)
  
  missing <- setdiff(pos.genes, rownames(mat))
  if (length(missing) > 0) {
    stop("Positive gene(s) not found: ", paste(missing, collapse = ", "))
  }
  
  # 2. Positive gene filtering
  expr.pos <- mat[pos.genes, , drop = FALSE] > expr.thresh
  keep.pos <- Matrix::colSums(expr.pos) >= min.genes
  cells <- colnames(mat)[keep.pos]
  
  if (verbose) message("After positive genes: ", length(cells))
  
  # Negative gene filtering
  if (!is.null(neg.genes)) {
    missing.neg <- setdiff(neg.genes, rownames(mat))
    if (length(missing.neg) > 0) {
      stop("Negative gene(s) not found: ", paste(missing.neg, collapse = ", "))
    }
    
    expr.neg <- mat[neg.genes, , drop = FALSE] > expr.thresh
    keep.neg <- Matrix::colSums(expr.neg) == 0
    cells <- intersect(cells, colnames(mat)[keep.neg])
    
    if (verbose) message("After negative genes: ", length(cells))
  }
  
  # Cluster restriction
  if (!is.null(clusters)) {
    clust <- if (is.null(cluster.col)) {
      as.character(Seurat::Idents(obj))
    } else {
      if (!cluster.col %in% colnames(obj@meta.data)) {
        stop("cluster.col not found in obj@meta.data")
      }
      as.character(obj@meta.data[[cluster.col]])
    }
    
    names(clust) <- colnames(obj)
    cells <- intersect(cells, names(clust)[clust %in% clusters])
    
    if (verbose) message("After cluster filter: ", length(cells))
  }
  
  # Metadata filtering (generalized)
  if (!is.null(meta.col)) {
    if (!meta.col %in% colnames(obj@meta.data)) {
      stop("meta.col not found in obj@meta.data")
    }
    if (is.null(meta.vals)) {
      stop("meta.vals must be provided when meta.col is used")
    }
    
    meta <- obj@meta.data[[meta.col]]
    names(meta) <- colnames(obj)
    
    # Attempt numeric coercion for character metadata
    meta.num <- suppressWarnings(as.numeric(gsub("[^0-9.-]", "", meta)))
    
    if (all(!is.na(meta.num)) && length(meta.vals) == 2) {
      keep.meta <- meta.num >= meta.vals[1] & meta.num <= meta.vals[2]
    } else {
      keep.meta <- meta %in% meta.vals
    }
    
    cells <- intersect(cells, names(meta)[keep.meta])
    
    if (verbose) message("After metadata filter: ", length(cells))
  }
  
  # Return
  if (return.obj) {
    return(subset(obj, cells = cells))
  }
  
  return(cells)
}


#'#' Compute OR Gene Dominance Scores in Single-Cell RNA-seq Data
#'
#' This function computes a dominance score for olfactory receptor (OR) genes in each cell
#' of a Seurat object. The dominance score quantifies how strongly the top-expressed OR
#' dominates over the second-most expressed OR, adjusted for sparse single-cell expression.
#'
#' @param object A Seurat object containing single-cell RNA-seq data.
#' @param or_genes A character vector of OR gene names to consider.
#' @param assay Character; the assay name in Seurat object (default: "RNA").
#' @param slot Character; the assay slot to use (default: "data").
#' @param epsilon Numeric; a small number to prevent division by zero (default: 1e-9).
#' @param scale_factor Numeric; multiplier for top OR expression to scale dominance (default: 10).
#' @param exponent Numeric; exponent to raise the scaled log expression (default: 1).
#'
#' @return A data frame with one row per cell containing:
#'   \describe{
#'     \item{top_expr}{Expression of the top-expressed OR gene.}
#'     \item{second_expr}{Expression of the second-highest OR gene.}
#'     \item{dominance_score}{Computed dominance score.}
#'     \item{cell}{Cell barcode or identifier.}
#'     \item{n_OR_expressed}{Number of OR genes with nonzero expression in the cell.}
#'     \item{PCW}{PCW metadata from the Seurat object.}
#'     \item{cluster}{Cell identity from Seurat object.}
#'     \item{cell_group}{Coarse cell grouping (mOSN, iOSN, Others).}
#'     \item{dominance_bin}{Binned dominance score (0–1, 1–2, …, >6).}
#'   }
#'
#' @examples
#' \dontrun{
#' Idents(sub) <- "ann_level_2"
#' or_genes <- intersect(ORs$Gene_name, rownames(sub))
#' df <- compute_OR_dominance(sub, or_genes)
#' head(df)
#' table(df$dominance_bin)
#' }
#'
#' @export
#' 
compute_OR_dominance <- function(
    object,
    or_genes,
    assay = "RNA",
    slot = "data",
    epsilon = 1e-9,
    scale_factor = 10,
    exponent = 1
){
  # --- Load required packages ---
  require(dplyr)
  require(Seurat)
  require(rlang)
  
  # --- Filter OR genes present in Seurat object ---
  or_genes <- intersect(or_genes, rownames(object))
  if(length(or_genes) == 0) stop("No OR genes found in the Seurat object.")
  
  # --- Extract expression matrix ---
  or_mat <- as.matrix(GetAssayData(object, assay = assay, slot = slot)[or_genes, , drop = FALSE])
  
  # --- Function to compute scaled dominance per cell ---
  compute_scaled_dominance <- function(x) {
    x_sorted <- sort(x, decreasing = TRUE)
    x1 <- x_sorted[1]
    x2 <- ifelse(length(x_sorted) >= 2, x_sorted[2], 0)
    
    if (x1 == 0) return(c(top_expr = 0, second_expr = 0, dominance_score = 0))
    
    # --- Scaled for sparse single-cell data ---
    score <- ((x1 - x2) / (x1 + x2 + epsilon)) * (log1p(scale_factor * x1) ^ exponent)
    return(c(top_expr = x1, second_expr = x2, dominance_score = score))
  }
  
  # --- Apply dominance calculation to all cells ---
  dom_vals <- t(apply(or_mat, 2, compute_scaled_dominance))
  dom_vals <- as.data.frame(dom_vals)
  
  # --- Ensure correct column names ---
  colnames(dom_vals) <- c("top_expr", "second_expr", "dominance_score")
  
  # --- Construct output dataframe ---
  df <- dom_vals %>%
    mutate(
      cell = colnames(or_mat),
      n_OR_expressed = colSums(or_mat > 0),
      PCW = object$PCW,
      cluster = Idents(object)
    ) %>%
    mutate(
      cell_group = case_when(
        cluster %in% c("GBC") ~ "GBC",
        cluster %in% c("INP") ~ "INP",
        cluster %in% c("iOSN") ~ "iOSN",
        TRUE ~ "Others"
      ),
      dominance_bin = cut(
        dominance_score,
        breaks = c(0, 1, 2, 3, 4, 5, 6, Inf),
        labels = c("0–1", "1–2", "2–3", "3–4", "4–5", "5–6", ">6"),
        right = FALSE,
        include.lowest = TRUE
      )
    )
  
  return(df)
}


run_slingshot <- function(
    obj,
    cluster_col = "ann_level_2",
    reduction = "pca",
    n_pcs = NULL,
    start_clust = "OHBC",
    end_clusts = NULL,
    approx_points = 100,
    verbose = TRUE
) {
  
  suppressPackageStartupMessages({
    library(slingshot)
    library(Seurat)
    library(dplyr)
    library(mgcv) 
    library(scales)
  })
  
  if (verbose) cat("Preparing clustering info...\n")
  Idents(obj) <- obj[[cluster_col, drop = TRUE]]
  clusters <- factor(as.character(Idents(obj)))
  
  if (verbose) cat(paste0("Using ", reduction, " reduction...\n"))
  emb <- Embeddings(obj, reduction)
  
  if (is.null(n_pcs)) {
    if (verbose) cat("Auto-detecting PCs (all available)...\n")
    pca_mat <- emb
  } else {
    pca_mat <- emb[, seq_len(n_pcs), drop = FALSE]
  }
  
  if (verbose) cat("Running Slingshot...\n")
  sds <- slingshot(
    pca_mat,
    clusterLabels = clusters,
    start.clus = start_clust,
    end.clus = end_clusts,
    approx_points = approx_points
  )
  
  if (verbose) cat("Extracting pseudotime...\n")
  pt <- slingPseudotime(sds)
  colnames(pt) <- paste0("Lineage", seq_len(ncol(pt)))
  
  # add pseudotime back into metadata
  pt_df <- as.data.frame(pt)
  obj <- AddMetaData(obj, metadata = pt_df)
  
  if (verbose) {
    cat("Done!\n")
    cat("Lineages detected: ", paste(colnames(pt), collapse = ", "),"\n")
    cat("Cells with pseudotime per lineage:\n")
    print(colSums(!is.na(pt)))
  }
  
  return(list(
    seurat = obj,
    sds = sds,
    pseudotime = pt
  ))
}


#' Split a Seurat object by metadata and fix matrix names
#'
#' @param object A Seurat object
#' @param split.by Metadata column name to split by
#' @return A named list of Seurat objects with @layers$counts fixed
#' @examples
#' seurat.list <- split_seurat_objects(merged, split.by = "BulkSample")
split_seurat <- function(object, split.by) {
  require(Seurat)
  
  objs <- Seurat::SplitObject(object, split.by = split.by)
  
  objs <- lapply(objs, function(sobj) {
    # Ensure @layers$counts exist
    assay <- sobj[["RNA"]]
    assay@layers$counts <- Seurat::GetAssayData(sobj, slot = "counts")
    sobj[["RNA"]] <- assay
    
    # Fix row/colnames
    sobj <- fix_seurat_matrix_names(sobj)
    return(sobj)
  })
  
  return(objs)
}


#' Add combined bulk/demux mapping to a Seurat object
#'
#' This function merges a combined cell-level mapping table (e.g., from souporcell, vireo, or scds)
#' into a Seurat object's metadata. Optionally, it can remove cells labeled as "doublet" or "unassigned".
#'
#' @param object A Seurat object.
#' @param combined_tsv Path to the combined mapping TSV file or a data.frame.
#' @param barcode_col Name of the barcode column in the TSV (default: "Barcode").
#' @param remove_doublets Logical, if TRUE removes cells labeled as "doublet" or "unassigned" in `BulkSample` (default: TRUE).
#' @return A Seurat object with updated metadata.
#' @examples
#' FN_S1256 <- add_bulk_mapping_to_seurat(FN_S1256, combined_tsv = "combined_results_with_bulk_mapping.tsv")
#' @export
add_bulk_mapping_to_seurat <- function(
    object,
    combined_tsv,
    barcode_col = "Barcode",
    remove_doublets = TRUE
) {
  library(Seurat)
  library(dplyr)
  
  # Load combined mapping
  if (is.character(combined_tsv)) {
    df <- read.delim(combined_tsv, stringsAsFactors = FALSE)
  } else if (is.data.frame(combined_tsv)) {
    df <- combined_tsv
  } else {
    stop("combined_tsv must be a file path or a data.frame")
  }
  
  # Ensure barcode column exists
  if (!barcode_col %in% colnames(df)) stop(paste("Column", barcode_col, "not found in combined mapping"))
  
  # Clean barcodes to match Seurat
  df$Barcode_clean <- gsub("_1$", "", df[[barcode_col]])
  rownames(df) <- df$Barcode_clean
  
  # Keep only barcodes present in Seurat
  common_cells <- intersect(colnames(object), rownames(df))
  if (length(common_cells) == 0) stop("No overlapping barcodes between Seurat object and mapping")
  df_sub <- df[common_cells, , drop = FALSE]
  
  # Add metadata
  object <- AddMetaData(object, df_sub)
  
  # Remove doublets/unassigned if requested
  if (remove_doublets && "BulkSample" %in% colnames(object@meta.data)) {
    bad_labels <- c("doublet", "doublets", "Doublet", "Doublets",
                    "unassigned", "Unassigned")
    
    current_idents <- Idents(object)
    to_remove <- intersect(bad_labels, levels(current_idents))
    
    if (length(to_remove) > 0) {
      object <- subset(object, idents = to_remove, invert = TRUE)
    }
  }
  return(object)
}


#' Compute gene-set dominance score per cell
#'
#' Computes a dominance score for a gene set in each cell based on the
#' relative expression of the top-expressed gene versus the second-highest
#' expressed gene. Optionally scales dominance by expression magnitude,
#' bins scores, assigns QC flags, and appends selected metadata for
#' downstream analysis.
#'
#' @param obj A \code{Seurat} object.
#' @param geneset Character vector of gene names defining the gene set.
#' @param assay Assay to use (default: \code{"RNA"}).
#' @param slot Expression slot to use (default: \code{"data"}).
#' @param epsilon Small constant to avoid division by zero.
#' @param scale.log Logical; whether to scale dominance by log-transformed
#'   top gene expression.
#' @param exponent Numeric exponent applied to log-scaled expression.
#' @param bins Numeric vector defining bins for adjusted dominance scores.
#' @param bin.labels Optional character labels for dominance bins.
#' @param strong.dom.thresh Threshold above which a cell is flagged as
#'   strongly dominant.
#' @param group.by Metadata column used for cell grouping
#'   (e.g. \code{"ann2"}).
#' @param idents Character vector of identity values to retain;
#'   all others are labeled \code{"Others"}.
#' @param add.meta Character vector of additional metadata columns
#'   (e.g. \code{c("sample","stage","sex")}) to append to output.
#' @param verbose Logical; whether to print informative messages.
#'
#' @details
#' For each cell, the dominance score is calculated as:
#'
#' \deqn{(Top1 - Top2) / (Top1 + Top2 + epsilon)}
#'
#' If \code{scale.log = TRUE}, the score is multiplied by
#' \code{log1p(Top1)^exponent}.
#'
#' Only the metadata columns specified in \code{group.by} and
#' \code{add.meta} are joined, ensuring minimal memory overhead.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{df}: Data frame with per-cell dominance metrics,
#'         bin assignments, QC flags, grouping identity, and selected metadata.
#'   \item \code{genes_detected}: Character vector of detected genes from
#'         \code{geneset}.
#' }
#'
#' @seealso \code{\link[Seurat]{GetAssayData}}
#'
#' @examples
#' \dontrun{
#' res <- domscore(
#'   obj = object,
#'   geneset = OR_genes,
#'   group.by = "ann2",
#'   idents = c("INP", "iOSN"),
#'   add.meta = c("sample", "stage", "sex"),
#'   bins = c(-Inf, 0.5, 1, 1.5, 2, Inf)
#' )
#'
#' head(res$df)
#' }
#'
#' @export
#' 
domscore <- function(
    obj,
    geneset = NULL,
    assay = "RNA",
    slot = "data",
    epsilon = 1e-9,
    scale.log = TRUE,
    exponent = 1,
    bins = c(-Inf, 0.1, 0.5, 1, 2, Inf),
    bin.labels = NULL,
    strong.dom.thresh = 1.5,
    group.by = NULL,        
    idents = NULL,      
    add.meta = NULL,    
    verbose = TRUE
) {
  
  require(Seurat)
  require(Matrix)
  require(dplyr)
  
  stopifnot(inherits(obj, "Seurat"))
  
  ## Gene detection
  all_genes <- rownames(obj[[assay]])
  detected_genes <- intersect(unique(geneset), all_genes)
  
  if (length(detected_genes) == 0)
    stop("No genes from geneset found in Seurat object.")
  
  expr_mat <- as.matrix(
    GetAssayData(obj, assay = assay, leyer = slot)[
      detected_genes, , drop = FALSE
    ]
  )
  
  ## Top1 / Top2 per cell
  
  compute_top2_cell <- function(v) {
    nz <- which(v != 0)
    if (length(nz) == 0) return(c(top1 = 0, top2 = 0))
    if (length(nz) == 1) return(c(top1 = v[nz], top2 = 0))
    ord <- order(v, decreasing = TRUE)
    c(top1 = v[ord[1]], top2 = v[ord[2]])
  }
  
  top2_mat <- apply(expr_mat, 2, compute_top2_cell)
  if (is.null(dim(top2_mat)))
    top2_mat <- matrix(top2_mat, nrow = 2)
  
  top1_vec <- top2_mat[1, ]
  top2_vec <- top2_mat[2, ]
  
  dominance <- (top1_vec - top2_vec) /
    (top1_vec + top2_vec + epsilon)
  
  adjust_score <- if (scale.log) {
    dominance * (log1p(top1_vec) ^ exponent)
  } else {
    dominance
  }
  
  ## Top gene & counts
  get_top_gene <- function(v) {
    if (all(v == 0)) return(NA_character_)
    rownames(expr_mat)[which.max(v)]
  }
  
  top_gene_vec <- apply(expr_mat, 2, get_top_gene)
  n_genes_expr <- Matrix::colSums(expr_mat > 0)
  
  ## Base dataframe
  df <- data.frame(
    cell = colnames(expr_mat),
    top_gene = top_gene_vec,
    top_expr = as.numeric(top1_vec),
    second_expr = as.numeric(top2_vec),
    dominance = as.numeric(dominance),
    adjust_score = as.numeric(adjust_score),
    n_genes_expr = as.integer(n_genes_expr),
    stringsAsFactors = FALSE
  )
  
  ## Selective metadata join
  meta <- obj@meta.data
  meta$cell <- rownames(meta)
  
  cols_to_add <- unique(c(group.by, add.meta))
  cols_to_add <- intersect(cols_to_add, colnames(meta))
  
  if (length(cols_to_add) > 0) {
    df <- left_join(
      df,
      meta[, c("cell", cols_to_add), drop = FALSE],
      by = "cell"
    )
  } else if (verbose) {
    message("No metadata columns added.")
  }
  
  ## Cell grouping
  ## Cell grouping (SAFE)
  if (!is.null(group.by) && group.by %in% colnames(df)) {
    
    group_vals <- as.character(df[[group.by]])
    
    if (!is.null(idents)) {
      df$ident <- ifelse(
        group_vals %in% idents,
        group_vals,
        "Others"
      )
    } else {
      df$ident <- group_vals
    }
    
  } else {
    df$ident <- "Others"
    if (verbose)
      warning("group.by not found; all cells set to 'Others'.")
  }
  
  ## Binning
  if (is.null(bin.labels)) {
    bin.labels <- sapply(seq_len(length(bins) - 1), function(i) {
      lb <- bins[i]; ub <- bins[i + 1]
      if (is.finite(lb) && is.finite(ub)) paste0(lb, "-", ub)
      else if (!is.finite(lb)) paste0("<", ub)
      else paste0(">=", lb)
    })
  }
  
  df$adjusted_bin <- cut(
    df$adjust_score,
    breaks = bins,
    labels = bin.labels,
    include.lowest = TRUE,
    right = FALSE
  )
  
  ## QC flags
  df$qc_flag <- "ok"
  df$qc_flag[df$n_genes_expressed == 0] <- "no_gene"
  df$qc_flag[df$top_expr == 0] <- "no_top_expr"
  df$qc_flag[df$adjust_score >= strong.dom.thresh] <- "strong_dominant"
  
  ## Output
  list(df = df, genes_detected = detected_genes)
}




#' Find Marker Genes in a Seurat Object
#' This function identifies marker genes for specified clusters in a Seurat object
#' using differential expression testing within each cluster.
#' Run Seurat::FindMarkers across clusters with multiple tests + consensus support + visualization
#'
#' @param object Seurat object
#' @param group.by Metadata column for clustering/identity
#' @param test.use Differential test(s) ("wilcox", "bimod", "roc", "t",
#'   "negbinom", "poisson", "LR", "MAST" or "all")
#' @param only.pos Keep only positive markers (default TRUE)
#' @param min.pct Minimum fraction of cells expressing a gene (default 0.25)
#' @param min.diff.pct Minimum difference in pct (default -Inf)
#' @param man.logfc.threshold LogFC threshold (default 0.25)
#' @param clusters.to.exclude Clusters to skip (default none)
#' @param max.cells.per.ident Downsample max cells per cluster (default Inf)
#' @param consensus Logical, if TRUE build consensus markers across tests
#' @param consensus_min_tests Integer, minimum number of tests a gene must be significant in
#' @param alpha Adjusted p-value threshold for significance (default 0.05)
#' @param return_both Logical, return both raw + consensus results (default FALSE)
#' @param plot_type "none", "bar", or "upset" (default "none")
#' @param plot_cluster Cluster name for upset plot (default first cluster)
#' @param ... Passed to Seurat::FindMarkers
#'
#' @return Data frame(s) of markers, optionally with ggplot object
#' @export
cellmarker <- function(object,
                       group.by,
                       assay = "RNA",
                       features = NULL,
                       test.use = 'wilcox',
                       only.pos = TRUE,
                       min.pct = 0.01,
                       min.diff.pct = -Inf,
                       logfc.threshold = 0.1,
                       clusters.to.exclude = c(),
                       max.cells.per.ident = Inf, 
                       latent.vars = NULL,
                       consensus = FALSE,
                       consensus_min_tests = 2,
                       alpha = 0.05,
                       return_both = FALSE,
                       plot_type = c("none","bar","upset"),
                       plot_cluster = NULL,
                       ...) {
  
  requireNamespace("ggplot2")
  if ("upset" %in% plot_type) requireNamespace("UpSetR")
  
  if (!inherits(object, "Seurat")) stop("The provided object is not a Seurat object.")
  
  if (!missing(group.by)) object <- Seurat::SetIdent(object = object, value = group.by)
  
  valid_tests <- c("wilcox", "bimod", "roc", "t", 
                   "negbinom", "poisson", "LR", "MAST")
  
  if (identical(test.use, "all")) {
    tests_to_run <- valid_tests
  } else {
    if (!all(test.use %in% valid_tests)) stop(paste("Invalid test.use. Choose from:", paste(valid_tests, collapse = ", ")))
    tests_to_run <- test.use
  }
  
  clusters.to.test <- sort(unique(object@active.ident))
  clusters.to.test <- setdiff(clusters.to.test, clusters.to.exclude)
  
  results_all <- list()
  
  for (test in tests_to_run) {
    message("=== Running test: ", test, " ===")
    joined <- data.frame()
    
    for (i in seq_along(clusters.to.test)) {
      cluster <- clusters.to.test[i]
      message(sprintf("[ %d / %d ] %s - cluster: %s ...", 
                      i, length(clusters.to.test), test, cluster))
      
      markers <- Seurat::FindMarkers(
        object,
        ident.1 = cluster, 
        ident.2 = NULL, 
        features = features,
        only.pos = only.pos, 
        assay = assay, 
        slot = "data",
        test.use = test,
        min.pct = min.pct,
        min.cells.group = 3,
        min.diff.pct = min.diff.pct,
        logfc.threshold = logfc.threshold,
        max.cells.per.ident = max.cells.per.ident, 
        latent.vars = latent.vars,
        ...
      )
      
      if (!is.null(markers) && nrow(markers) > 0) {
        markers$cluster <- cluster
        markers$gene <- rownames(markers)
        markers$test.use <- test
        rownames(markers) <- NULL
        joined <- rbind(joined, markers)
      } else {
        message(sprintf("Skipping cluster '%s': no markers found.", cluster))
      }
    }
    
    results_all[[test]] <- joined
  }
  
  all_results <- dplyr::bind_rows(results_all)
  
  #- Single test shortcut-
  if (length(tests_to_run) == 1 && !consensus && !return_both) {
    return(all_results)
  }
  
  #- Consensus / visualization-
  consensus_df <- NULL
  if (consensus || return_both) {
    consensus_df <- all_results %>%
      dplyr::filter(p_val_adj <= alpha) %>%
      dplyr::group_by(cluster, gene) %>%
      dplyr::summarise(
        n_tests = dplyr::n_distinct(test.use),
        mean_logFC = mean(avg_log2FC, na.rm = TRUE),
        min_p_val_adj = min(p_val_adj, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_tests >= consensus_min_tests)
  }
  
  plot_type <- match.arg(plot_type)
  plt <- NULL
  
  if (plot_type == "bar") {
    marker_counts <- all_results %>%
      dplyr::filter(p_val_adj <= alpha) %>%
      dplyr::group_by(cluster, test.use) %>%
      dplyr::summarise(n_genes = dplyr::n(), .groups = "drop")
    
    plt <- ggplot2::ggplot(marker_counts, 
                           ggplot2::aes(x = cluster, y = n_genes, fill = test.use)) +
      ggplot2::geom_bar(stat="identity", position="dodge") +
      plot_theme(theme.type = "classic", x.angle = 45) +
      ggplot2::labs(title="Significant markers per test per cluster", 
                    y="Number of markers", x="Cluster")
  }
  
  if (plot_type == "upset") {
    if (is.null(plot_cluster)) plot_cluster <- clusters.to.test[1]
    df_upset <- all_results %>%
      dplyr::filter(cluster == plot_cluster, p_val_adj <= alpha) %>%
      dplyr::select(gene, test.use) %>%
      dplyr::distinct()
    
    mat <- table(df_upset$gene, df_upset$test.use) > 0
    plt <- UpSetR::upset(UpSetR::fromMatrix(mat), 
                         mainbar.y.label = paste("Overlap of markers -", plot_cluster))
  }
  
  if (return_both) return(list(raw = all_results, consensus = consensus_df, plot = plt))
  if (consensus) return(list(consensus = consensus_df, plot = plt))
  
  return(list(raw = all_results, plot = plt))
}


#' Identify Temporally Regulated Genes Along Pseudotime
#'
#' This function identifies genes whose expression varies significantly
#' along pseudotime within selected cell types. It supports Slingshot
#' pseudotime objects (`sds_obj`) or a user-provided pseudotime matrix (`pt_matrix`).
#' Genes are filtered by variability, minimum expression, pseudotime range,
#' and minimal dynamic range (logFC). Candidate genes are tested via GAM
#' (`mgcv::gam`) for significant smooth changes along pseudotime.
#'
#' @param object Seurat object containing expression data.
#' @param sds_obj Optional \code{SlingshotDataSet} object with pseudotime estimates.
#' @param pt_matrix Optional numeric matrix of pseudotime values
#'   (cells × lineages). Used if \code{sds_obj} is not provided.
#' @param annotation_level Metadata column name specifying cell type labels
#'   (default: `"ann_level_2"`).
#' @param selected_celltypes Vector of cell types to include (default:
#'   `c("OHBC","GBC","INP","iOSN")`).
#' @param assay Assay name in Seurat object (default: `"RNA"`).
#' @param expr_slot Expression slot to use (default: `"data"`).
#' @param n_var_genes Number of highly variable genes to preselect for testing.
#' @param min_cells_expressed Minimum number of cells a gene must be expressed in.
#' @param min_pct_cells Minimum fraction of cells expressing a gene.
#' @param min_pt_range Minimum pseudotime range over which a gene is expressed.
#' @param min_logFC Minimum dynamic range (max − min expression) required.
#' @param gam_k Number of basis functions used in GAM smoothing (default: 6).
#' @param top_n Optional integer. If provided, extracts the top-n significant genes for heatmap use.
#' @param qval_cutoff Adjusted p-value (FDR) threshold for significance.
#' @param verbose Logical; print progress messages.
#'
#' @return A list containing:
#' \describe{
#'   \item{lineage_results}{A list of data frames per lineage,
#'     containing p-values, q-values, and significant genes.}
#'   \item{heatmap_matrices}{A list of matrices per lineage:
#'     raw expression, scaled expression (all significant genes),
#'     scaled expression for top N genes, and ordered cell order.}
#'   \item{pt_used}{The pseudotime matrix used for analysis.}
#' }
#'
#' @details
#' Filtering steps:
#' \enumerate{
#'   \item Highest-variable genes preselected.
#'   \item Genes must be expressed in sufficient cells.
#'   \item Genes must span a minimum pseudotime window.
#'   \item Genes must show sufficient dynamic range (logFC).
#' }
#'
#' GAM testing is performed using: \code{y ~ s(pt, k = gam_k)}.
#'
#' @import mgcv Seurat dplyr matrixStats tibble
#'
#' @examples
#' \dontrun{
#' result <- identify_temporal_genes(
#'   object = object,
#'   sds_obj = slingshot_obj,
#'   selected_celltypes = c("iOSN")
#' )
#' }
#'
#' @export
#' 
identify_temporal_genes <- function(object,
                                    sds_obj = NULL,
                                    pt_matrix = NULL,
                                    annotation_level = "ann_level_2",
                                    selected_celltypes = c("OHBC","GBC","INP","iOSN"),
                                    assay = "RNA",
                                    expr_slot = "data",
                                    n_var_genes = 500,
                                    min_cells_expressed = 10,
                                    min_pct_cells = 0.05,
                                    min_pt_range = 0.15,
                                    min_logFC = 0.5,
                                    gam_k = 6,
                                    top_n = NULL,      # optional
                                    qval_cutoff = 0.05,
                                    verbose = TRUE) {
  
  require(mgcv); require(Seurat); require(dplyr); require(matrixStats); require(tibble)
  
  if (is.null(sds_obj) && is.null(pt_matrix)) stop("Provide sds_obj or pt_matrix")
  
  cells_sel <- WhichCells(object, idents = selected_celltypes)
  pt_full <- if(is.null(pt_matrix)) slingPseudotime(sds_obj) else pt_matrix
  common_cells <- intersect(rownames(pt_full), colnames(object))
  cells_use <- intersect(common_cells, cells_sel)
  pt_full <- pt_full[cells_use,, drop=FALSE]
  
  expr_all <- as.matrix(GetAssayData(object, assay=assay, slot=expr_slot))[,cells_use, drop = FALSE]
  var_feats <- rownames(expr_all)[order(rowVars(expr_all), decreasing=TRUE)][1:min(n_var_genes, nrow(expr_all))]
  
  lineage_results <- list()
  heatmap_matrices <- list()
  
  for(li in seq_len(ncol(pt_full))){
    
    lineage_name <- colnames(pt_full)[li]
    if(verbose) message("Processing ", lineage_name, " ...")
    
    pts <- pt_full[, li]
    valid_cells <- names(pts)[!is.na(pts)]
    pts <- pts[valid_cells]
    expr <- expr_all[var_feats, valid_cells, drop = FALSE]
    
    # filtering
    expressed_counts <- rowSums(expr > 0)
    cells_min <- max(min_cells_expressed, ceiling(min_pct_cells * length(valid_cells)))
    keep1 <- names(expressed_counts)[expressed_counts >= cells_min]
    if(length(keep1) == 0) {
      if(verbose) message("  No genes after expression filter; skipping")
      lineage_results[[lineage_name]] <- tibble()
      heatmap_matrices[[lineage_name]] <- list(expr_raw = NULL, expr_scaled_top = NULL, ordered_cells = names(sort(pts)), top_genes = NULL)
      next
    }
    
    expr_pt_ranges <- t(apply(expr[keep1, , drop=FALSE] > 0, 1, function(x) {
      r <- range(pts[x], na.rm=TRUE)
      if (any(is.na(r))) return(c(NA, NA))
      r
    }))
    expr_pt_range <- expr_pt_ranges[,2] - expr_pt_ranges[,1]
    keep2 <- names(expr_pt_range)[!is.na(expr_pt_range) & expr_pt_range >= min_pt_range]
    if(length(keep2) == 0) {
      if(verbose) message("  No genes after pseudotime-range filter; skipping")
      lineage_results[[lineage_name]] <- tibble()
      heatmap_matrices[[lineage_name]] <- list(expr_raw = NULL, expr_scaled_top = NULL, ordered_cells = names(sort(pts)), top_genes = NULL)
      next
    }
    
    logfc <- matrixStats::rowMaxs(expr[keep2, , drop=FALSE]) - matrixStats::rowMins(expr[keep2, , drop=FALSE])
    keep3 <- names(logfc)[!is.na(logfc) & logfc >= min_logFC]
    if(length(keep3) == 0) {
      if(verbose) message("  No genes after logFC filter; skipping")
      lineage_results[[lineage_name]] <- tibble()
      heatmap_matrices[[lineage_name]] <- list(expr_raw = NULL, expr_scaled_top = NULL, ordered_cells = names(sort(pts)), top_genes = NULL)
      next
    }
    
    if(verbose) message("  Candidate genes: ", length(keep3))
    
    # GAM testing
    pvals <- setNames(rep(NA_real_, length(keep3)), keep3)
    for(g in keep3){
      d <- data.frame(y = as.numeric(expr[g, ]), pt = as.numeric(pts))
      if(sd(d$y, na.rm = TRUE) == 0) next
      m <- tryCatch(mgcv::gam(y ~ s(pt, k = gam_k), data = d), error = function(e) NULL)
      if(!is.null(m)){
        sst <- summary(m)$s.table
        if(!is.null(sst)) {
          pvals[g] <- if("p-value" %in% colnames(sst)) sst[1, "p-value"] else sst[1, ncol(sst)]
        }
      }
    }
    
    res_df <- tibble(id = names(pvals), raw_p = as.numeric(pvals)) %>%
      filter(!is.na(raw_p)) %>%
      mutate(qval = p.adjust(raw_p, method = "fdr")) %>%
      arrange(qval)
    
    # all significant genes (full table)
    sig_df <- res_df %>% filter(qval <= qval_cutoff)
    if(verbose) message("  Significant genes: ", nrow(sig_df))
    
    ord_cells <- names(sort(pts))
    
    # Always store expr_raw for the significant genes (could be zero rows)
    if(nrow(sig_df) > 0) {
      expr_raw_sig <- as.matrix(expr[sig_df$id, ord_cells, drop = FALSE])
      # scale for inspection (not required for plotting; plot function may re-scale or smooth)
      expr_scaled_all_sig <- t(scale(t(expr_raw_sig)))
      expr_scaled_all_sig[is.na(expr_scaled_all_sig)] <- 0
    } else {
      expr_raw_sig <- NULL
      expr_scaled_all_sig <- NULL
    }
    
    # If user requested top_n, compute expr_scaled_top and top_genes
    if(!is.null(top_n) && nrow(sig_df) > 0) {
      top_genes <- head(sig_df$id, top_n)
      expr_sub <- as.matrix(expr[top_genes, ord_cells, drop = FALSE])
      expr_scaled_top <- t(scale(t(expr_sub))); expr_scaled_top[is.na(expr_scaled_top)] <- 0
    } else {
      top_genes <- NULL
      expr_scaled_top <- NULL
    }
    
    lineage_results[[lineage_name]] <- sig_df
    heatmap_matrices[[lineage_name]] <- list(
      expr_raw = expr_raw_sig,            # raw expr for *all* significant genes
      expr_scaled_top = expr_scaled_top,  # scaled expr matrix for top_n (NULL if not requested)
      expr_scaled_all = expr_scaled_all_sig, # scaled expr matrix for all sig genes (may be NULL)
      ordered_cells = ord_cells,
      top_genes = top_genes
    )
  }
  
  return(list(
    lineage_results = lineage_results,
    heatmap_matrices = heatmap_matrices,
    pt_used = pt_full
  ))
}


plot_theme <- function(
    theme.style = c("minimal","classic","bw","test","void","dirty","gray"),
    font.size = 8,
    xy.val = TRUE, 
    x.angle = 0, 
    hjust = NULL,
    vjust = NULL,
    xlab = TRUE, 
    ylab = TRUE,
    xy.lab = TRUE,
    facet.face = "bold", 
    ttl.face = "bold",
    txt.face = c("plain","italic","bold"),
    ttl.pos = c("center","left","right"),
    x.ttl = TRUE, 
    y.ttl = TRUE,
    ticks = NULL,      
    line = NULL,       
    border = NULL,      
    grid.major = NULL, 
    grid.minor = NULL,  
    panel.fill = "white", 
    facet.bg = TRUE,
    mode = c("light","dark"),
    leg.pos = "right",
    leg.dir = "vertical",
    leg.size = 8,
    leg.ttl = 8, 
    leg.ttl.size = 8,
    leg.just = "center",
    leg.ttl.text = NULL,
    ...
) {
  require(ggplot2)
  
  theme.style <- match.arg(theme.style)
  ttl.pos <- match.arg(ttl.pos)
  txt.face <- match.arg(txt.face)
  mode <- match.arg(mode)
  
  # Canonical line width (THIS FIXES YOUR PROBLEM)
  lw <- 0.3
  
  if (is.null(line)) {
    line <- theme.style == "classic"
  }
  
  # Colors
  if (mode == "light") {
    col.txt   <- "#1A1A1A"
    col.grid  <- "#D9D9D9"
    col.panel <- panel.fill
    col.strip <- "#EFEFEF"
  } else {
    col.txt   <- "#DDDDDD"
    col.grid  <- "#444444"
    col.panel <- "#1E1E1E"
    col.strip <- "#383838"
  }
  
  # Automatic x-label alignment
  if (is.null(hjust) || is.null(vjust)) {
    if (x.angle == 0)   { hjust <- .5; vjust <- 1 }
    else if (x.angle == 45) { hjust <- 1; vjust <- 1 }
    else if (x.angle == 90) { hjust <- 1; vjust <- .5 }
    else if (x.angle == 270){ hjust <- 0; vjust <- .5 }
    else { hjust <- 1; vjust <- 1 }
  }
  
  ttl.pos <- switch(ttl.pos, left = 0, center = .5, right = 1)
  
  # Base theme components
  base <- theme(
    text = element_text(color = col.txt, size = font.size, family = "Helvetica"),
    axis.text.x = element_text(color = col.txt, size = font.size),
    axis.text.y = element_text(color = col.txt, size = font.size),
    axis.title  = element_text(size = font.size),
    plot.title  = element_text(
      hjust = ttl.pos, face = ttl.face,
      size = font.size + 1, color = col.txt
    ),
    strip.text = element_text(face = facet.face, color = col.txt),
    legend.title = element_text(size = leg.ttl.size, face = "bold"),
    legend.text  = element_text(size = leg.size),
    legend.position = leg.pos,
    legend.key.height = unit(.4, "cm"),
    legend.key.width  = unit(.4, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.box = "vertical",
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(1,1,1,1),
    ...
  )
  
  # Preset themes
  preset <- switch(
    theme.style,
    minimal = theme_minimal(base_size = font.size),
    classic = theme_classic(base_size = font.size),
    bw      = theme_bw(base_size = font.size),
    test    = theme_test(base_size = font.size),
    void    = theme_void(base_size = font.size),
    dirty   = theme_minimal(base_size = font.size) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank()
      ),
    gray    = theme_gray(base_size = font.size) +
      theme(
        panel.background = element_rect(fill = "#EDEDED", color = NA),
        panel.grid.major = element_line(color = "#CCCCCC", linewidth = lw),
        panel.grid.minor = element_line(color = "#DDDDDD", linewidth = lw/2)
      )
  )
  
  th <- preset + base
  
  # Normalize panel borders for border-based themes
  if (theme.style %in% c("bw", "test", "gray")) {
    th <- th + theme(
      panel.border = element_rect(
        linewidth = lw,
        color = col.txt,
        fill = NA
      ),
      axis.line = element_blank()   # ⬅️ critical
    )
  }
  
  # Axis lines (classic-style)
  if (line && theme.style == "classic") {
    th <- th + theme(
      axis.line.x = element_line(color = col.txt, linewidth = lw),
      axis.line.y = element_line(color = col.txt, linewidth = lw)
    )
  } else {
    th <- th + theme(axis.line = element_blank())
  }
  
  # X-axis angle
  if (xy.val && xlab) {
    th <- th + theme(
      axis.text.x = element_text(angle = x.angle, hjust = hjust, vjust = vjust)
    )
  }
  
  # Label / tick / grid overrides
  if (!xy.lab) th <- th + theme(axis.text = element_blank(), axis.ticks = element_blank())
  if (!xlab)   th <- th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  if (!ylab)   th <- th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  if (!x.ttl)  th <- th + theme(axis.title.x = element_blank())
  if (!y.ttl)  th <- th + theme(axis.title.y = element_blank())
  
  if (!is.null(ticks)) {
    th <- th + if (ticks)
      theme(axis.ticks = element_line(color = col.txt, linewidth = lw))
    else
      theme(axis.ticks = element_blank())
  }
  
  if (!is.null(border)) {
    th <- th + if (border)
      theme(panel.border = element_rect(color = col.grid, fill = NA, linewidth = lw))
    else
      theme(panel.border = element_blank())
  }
  
  if (!is.null(grid.major)) {
    th <- th + if (grid.major)
      theme(panel.grid.major = element_line(color = col.grid, linewidth = lw))
    else
      theme(panel.grid.major = element_blank())
  }
  
  if (!is.null(grid.minor)) {
    th <- th + if (grid.minor)
      theme(panel.grid.minor = element_line(color = col.grid, linewidth = lw/2))
    else
      theme(panel.grid.minor = element_blank())
  }
  
  if (!facet.bg) {
    th <- th + theme(strip.background = element_blank())
  }
  
  if (!is.null(leg.ttl.text)) {
    th <- th + labs(color = leg.ttl.text)
  }
  
  # Void cleanup
  # Void cleanup
  if (theme.style == "void") {
    th <- th + theme(
      axis.text.x  = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks   = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line    = element_blank(),
      panel.grid   = element_blank(),
      panel.border = element_blank(),
      strip.text   = element_blank(),
      strip.background = element_blank()
    )
    line <- FALSE
    ticks <- FALSE
    border <- FALSE
    grid.major <- FALSE
    grid.minor <- FALSE
    facet.bg <- FALSE
  }
  
  # Dirty theme
  if (theme.style == "dirty") {
    line       <- TRUE
    ticks      <- FALSE
    border     <- FALSE
    grid.major <- FALSE
    grid.minor <- FALSE
    facet.bg   <- FALSE
  }
  
  th
}


#' Generate a Color Palette with Presets and Interpolation
#'
#' This function generates a vector of colors based on predefined palette presets
#' or a user-supplied set of base colors. Colors can be interpolated in multiple
#' color spaces and adjusted for saturation, lightness, and vividness.
#'
#' @param n Integer. Number of colors to generate.
#' @param preset Character. One of \code{"bright"}, \code{"pastel"}, \code{"warm"},
#'   \code{"cool"}, \code{"contrast"}, \code{"earth"}, \code{"base"}, or \code{"custom"}.
#'   If \code{"custom"} is chosen, \code{base_colors} must be supplied.
#' @param base_colors Character vector of HEX color codes. Required if
#'   \code{preset = "custom"}.
#' @param space Character. Color interpolation space. One of \code{"Lab"},
#'   \code{"rgb"}, or \code{"HCL"}.
#' @param oversample_factor Numeric. Factor by which to oversample colors before filtering.
#'   Larger values create smoother gradients. Default is \code{1.3}.
#' @param remove_gray Logical. If \code{TRUE} (default), removes grayish colors
#'   (low saturation).
#' @param reverse Logical. If \code{TRUE}, reverses the order of colors.
#' @param adjust_saturation Numeric multiplier for color saturation (chroma) in HCL space.
#'   Default is \code{1} (no change).
#' @param adjust_lightness Numeric multiplier for lightness in HCL space.
#'   Default is \code{1} (no change).
#'
#' @details
#' This function provides flexibility for both discrete and continuous color needs.
#' If the number of requested colors (\code{n}) is greater than the number of colors
#' in the base palette, the palette is interpolated in the chosen color space.
#'
#' The \code{remove_gray} option filters out low-chroma colors that appear grayish.
#'
#' @return A character vector of HEX color codes.
#'
#' @examples
#' # Get 5 bright colors
#' custom_palette(5, preset = "bright", space = "Lab")
#'
#' # Get 15 colors from the custom "base" preset
#' custom_palette(15, preset = "base", space = "HCL")
#'
#' # Use a fully custom palette
#' my_cols <- c("#123456", "#abcdef", "#ff0000")
#' custom_palette(10, preset = "custom", base_colors = my_cols)
#'
#' @import colorspace
#' @export
custom_palette <- function(n,
                           preset = c("base","bright", "pastel", "warm", "cool", "contrast", "earth", "custom"),
                           base_colors = NULL,
                           space = c("Lab", "rgb", "HCL"),
                           oversample_factor = 1.3,
                           remove_gray = TRUE,
                           reverse = FALSE,
                           adjust_saturation = 1,
                           adjust_lightness = 1) {
  
  library(colorspace)
  
  # Presets
  presets <- list(
    bright = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
               "#FFFF33", "#A65628", "#F781BF", "#999999"),
    pastel = c("#FDB462", "#B3DE69", "#BC80BD", "#CCEBC5", "#FFED6F",
               "#FB9A99", "#B2DF8A", "#CAB2D6", "#FFFFB3"),
    warm = c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#FDD49E",
             "#F4A582", "#D6604D", "#B2182B"),
    cool = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
             "#0571B0", "#74ADD1", "#ABD9E9", "#E0F3F8"),
    contrast = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                 "#66A61E", "#E6AB02", "#A6761D"),
    earth = c("#A6611A", "#DFC27D", "#80CDC1", "#018571",
              "#E5E5E5", "#F5F5F5", "#B2182B", "#D6604D"),
    # base = c(
    #   "#E41A1C", "#57A156", "#06A5FF", "#7289da", "#8D5B96", "#CB7647", "#F38E38",
    #   "#F781BE", "#CC95C8", "#B27E85", "#9A6242", "#5FA3C9", "#3766A4", "#204E75",
    #   "#1B9D77", "#86CC84", "#D3EC90", "#FBF583", "#E7C715", "#F0A957", "#F57994",
    #   "#E7298A", "#A90D55", "#52587E", "#17CDD3", "#8ECDE0", "#BC6298", "#AE2373",
    #   "#5E4EA1", "#7E8D86", "#507C51", "#1F5917", "#BEC603", "#C5DD3B", "#A8DA83",
    #   "#8DD3C7"
    # ),
    base = c(
      "#E41A1C", "#68618B", "#409388", "#57A156", "#8D5B96", "#CB7647", "#F38E38",
      "#F781BE", "#CC95C8", "#B27E85", "#9A6242", "#5FA3C9", "#3766A4", "#204E75",
      "#1B9D77", "#86CC84", "#D3EC90", "#FBF583", "#E7C715", "#F0A957", "#F57994",
      "#E7298A", "#A90D55", "#52587E", "#17CDD3", "#8ECDE0", "#BC6298", "#AE2373",
      "#5E4EA1", "#7E8D86", "#507C51", "#1F5917", "#BEC603", "#C5DD3B", "#A8DA83",
      "#8DD3C7"
    ),
    custom = NULL
  )
  
  # Pick preset
  preset <- match.arg(preset)
  if (preset != "custom") {
    base_colors <- presets[[preset]]
  } else if (is.null(base_colors)) {
    stop("For preset = 'custom', you must provide base_colors.")
  }
  
  # Match space argument
  space <- match.arg(space)
  
  # Short-circuit if enough colors
  if (n <= length(base_colors)) {
    cols <- base_colors[1:n]
    if (reverse) cols <- rev(cols)
    return(cols)
  }
  
  # Oversample
  extra_colors <- ceiling(n * oversample_factor)
  
  # Interpolation
  if (space %in% c("rgb", "Lab")) {
    extended_colors <- grDevices::colorRampPalette(base_colors, space = space)(extra_colors)
  } else if (space == "HCL") {
    extended_colors <- grDevices::colorRampPalette(base_colors, space = "hcl")(extra_colors)
  }
  
  # Remove grayish
  if (remove_gray) {
    is_grayish <- function(col) {
      rgb <- grDevices::col2rgb(col)
      sd(rgb) < 15
    }
    extended_colors <- Filter(function(c) !is_grayish(c), extended_colors)
  }
  
  # Adjust saturation / lightness
  if (adjust_saturation != 1 || adjust_lightness != 1) {
    hcl_vals <- colorspace::coords(
      methods::as(colorspace::hex2RGB(extended_colors), "polarLUV")
    )
    hcl_vals[, "C"] <- pmax(0, hcl_vals[, "C"] * adjust_saturation)
    hcl_vals[, "L"] <- pmax(0, pmin(100, hcl_vals[, "L"] * adjust_lightness))
    extended_colors <- colorspace::hex(colorspace::polarLUV(hcl_vals))
  }
  
  # Take first n
  final_colors <- extended_colors[1:n]
  if (reverse) final_colors <- rev(final_colors)
  
  return(final_colors)
}



#' Cell map visualization for Seurat objects
#'
#' @param object Seurat object
#' @param group.by Column(s) in metadata to group cells
#' @param reduction Dimensionality reduction to use (default "umap")
#' @param dims Dimensions to plot (1,2) or 3 for 3D
#' @param shuffle Logical, shuffle cells
#' @param raster Logical, rasterize plot
#' @param alpha Point transparency
#' @param repel Logical, repel labels
#' @param n.cells Logical, show number of cells in labels
#' @param label Logical, show cluster labels
#' @param label.size Cluster label size
#' @param label.face Cluster label font face
#' @param colors Named vector of colors
#' @param figplot Logical, minimal figure plot for figure panels
#' @param no.axes Logical, hide axes
#' @param plot.ttl Plot title
#' @param legend Logical, show legend
#' @param leg.ttl Legend title
#' @param item.size Legend item size
#' @param leg.pos Legend position
#' @param leg.just Legend justification
#' @param leg.dir Legend direction
#' @param leg.size Legend font size
#' @param leg.ncol Legend number of columns
#' @param item.border Logical, border around legend items
#' @param font.size Base font size
#' @param pt.size Point size
#' @param dark Logical, dark theme
#' @param total.cells Logical, include total cells in title
#' @param threeD Logical, 3D plot
#' @param theme Theme name
#' @param facet.bg Logical, add facet background
#' @param ... Additional arguments passed to DimPlot
#' @return ggplot or plotly object
#' @export
#' 
cellmap <- function(
    object, 
    group.by = NULL, 
    reduction = "umap", 
    dims = c(1,2), 
    shuffle = FALSE,
    raster = NULL,
    raster.dpi = c(512, 512), 
    alpha = 1, 
    repel = FALSE, 
    n.cells = TRUE,
    label = FALSE, 
    label.size = 3.5, 
    label.face = "plain", 
    cols = NULL, 
    figplot = FALSE, 
    no.axes = FALSE, 
    plot.ttl = NULL, 
    legend = TRUE,
    leg.ttl = NULL, 
    leg.ttl.size = font.size,
    item.size = 3.5, 
    leg.pos = "right", 
    leg.just = "center",
    leg.dir = "vertical", 
    leg.size = 10, 
    leg.ncol = NULL, 
    item.border = TRUE,
    font.size = 10, 
    pt.size = 0.5, 
    dark = FALSE, 
    total.cells = FALSE,
    threeD = FALSE, 
    theme.style = "classic", 
    facet.bg = FALSE,
    ...
) {
  
  if (!is.null(list(...)$theme)) {
    theme.style <- list(...)$theme
  }
  
  # Helper: prepare object and colors
  .prepare_object <- function(obj, group, cols=NULL){
    stopifnot(inherits(obj,"Seurat"))
    if(is.null(group)) group <- "ident"
    
    if(group=="ident"){
      obj@meta.data$ident <- Idents(obj)
    } else if(!group %in% colnames(obj@meta.data)){
      stop(paste("Grouping column", group, "not found."))
    }
    
    # Drop unused levels
    obj@meta.data[[group]] <- droplevels(factor(obj@meta.data[[group]]))
    values <- as.character(obj@meta.data[[group]])
    values[is.na(values)] <- "Unknown"
    
    obj@meta.data[[group]] <- factor(values)
    Idents(obj) <- obj@meta.data[[group]]
    
    levels_group <- levels(obj@meta.data[[group]])
    
    if(is.null(cols)){
      cols <- custom_palette(length(levels_group))
      names(cols) <- levels_group
      if("Unknown" %in% levels_group) cols["Unknown"] <- "gray70"
    } else {
      if(is.null(names(cols))) cols <- setNames(cols[seq_along(levels_group)], levels_group)
      missing <- setdiff(levels_group, names(cols))
      if(length(missing)) cols[missing] <- "gray70"
      cols <- cols[levels_group]
    }
    
    list(obj=obj, cols=cols, levels_group=levels_group)
  }
  
  # 3D plotting helper
  .plot_3D <- function(obj, dims, cols, alpha, pt.size, label, label.size, label.face, n.cells){
    emb <- obj@reductions[[reduction]]@cell.embeddings
    df <- data.frame(x=emb[, dims[1]], y=emb[, dims[2]], z=emb[, dims[3]], cluster=Idents(obj))
    hover_labels <- if(n.cells){
      tbl <- table(df$cluster)
      paste0(df$cluster, " (", tbl[as.character(df$cluster)], ")")
    } else as.character(df$cluster)
    
    p3d <- plotly::plot_ly(df, x=~x, y=~y, z=~z, color=~cluster, colors=cols,
                           type="scatter3d", mode="markers",
                           marker=list(size=pt.size, opacity=alpha, line=list(width=0)),
                           text=hover_labels, hoverinfo="text")
    if(label){
      centers <- df %>% dplyr::group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), z=median(z))
      p3d <- p3d %>% plotly::add_text(data=centers, x=~x, y=~y, z=~z, text=~cluster, textposition="top center")
    }
    p3d
  }
  
  if(length(group.by) > 1){
    plots <- lapply(group.by, function(g){
      cellmap(
        object = object,
        group.by = g,
        shuffle = shuffle,
        raster = raster,
        alpha = alpha,
        repel = repel,
        reduction = reduction,
        dims = dims,
        n.cells = n.cells,
        label = label,
        label.size = label.size,
        label.face = label.face,
        cols = cols,
        figplot = figplot,
        plot.ttl = g,
        legend = legend,
        leg.ttl = g,
        item.size = item.size,
        leg.pos = leg.pos,
        leg.just = leg.just,
        leg.dir = leg.dir,
        leg.ncol = leg.ncol,
        font.size = font.size,
        item.border = item.border,
        pt.size = pt.size,
        dark = dark,
        total.cells = total.cells,
        threeD = threeD,
        theme.style = theme.style,
        facet.bg = facet.bg,
        ...
      )
    })
    return(patchwork::wrap_plots(plots))
  }
  # Prepare object & colors
  prep <- .prepare_object(object, group.by, cols)
  object <- prep$obj
  cols <- prep$cols
  levels_group <- prep$levels_group
  if(is.null(leg.ncol)) leg.ncol <- if(length(levels_group) > 18) 2 else 1
  
  # Validate reduction
  if(!(reduction %in% names(object@reductions))){
    stop(paste0("Reduction '", reduction, "' not found. Available: ", paste(names(object@reductions), collapse=", ")))
  }
  emb <- object@reductions[[reduction]]@cell.embeddings
  if(max(dims) > ncol(emb)) stop("Selected dims exceed available dimensions in reduction.")
  
  # 3D plotting
  if(threeD || length(dims) == 3) return(.plot_3D(object, dims, cols, alpha, pt.size, label, label.size, label.face, n.cells))
  
  # 2D plotting
  plt <- Seurat::DimPlot(object, group.by=group.by, shuffle=shuffle, raster=raster, pt.size=pt.size,
                         repel=repel, alpha=alpha, reduction=reduction, dims=dims, raster.dpi=raster.dpi,...)
  
  # reset Seurat's forced theme_classic
  plt <- plt + ggplot2::theme_void()
  
  present_levels <- levels(droplevels(object@active.ident))
  
  if (n.cells) {
    cell.nb <- table(object@active.ident)[present_levels]
    clust.lab <- paste0(present_levels, " (", cell.nb, ")")
  } else {
    clust.lab <- present_levels
  }
  
  cols_use <- cols[present_levels]
  
  leg.ttl <- if(is.null(leg.ttl)) group.by else leg.ttl
  
  plt <- plt + ggplot2::scale_color_manual(
    breaks = present_levels,
    labels = clust.lab,
    values = cols_use
  )
  
  
  if (legend) {
    plt <- plt & ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = if(item.border)
          list(size = item.size, shape = 21, color = "black",
               stroke = 0.2, fill = unname(cols))
        else
          list(size = item.size),
        ncol = leg.ncol,
        title = leg.ttl,
        keyheight = grid::unit(0.25,"cm"),
        keywidth  = grid::unit(0.25,"cm")
      )
    )
  } else {
    plt <- plt & ggplot2::guides(color = "none")
  }
  
  # Plot title with total cells
  if(total.cells){
    plot.ttl <- paste0(plot.ttl, " (n=", format(ncol(object), big.mark=","), ")")
  }
  #if(!is.null(plot.ttl)) plt <- plt + labs(title = plot.ttl)
  if(!is.null(plot.ttl)) plt <- plt + labs(title = plot.ttl) else plt <- plt + labs(title = NULL)
  
  # Set default legend title size
  if (is.null(leg.ttl.size)) leg.ttl.size <- font.size
  
  # Apply plot_theme using do.call
  theme_args <- list(
    theme.style = theme.style,
    font.size = font.size,
    leg.size  = leg.size,
    leg.pos   = leg.pos,
    leg.dir   = leg.dir,
    leg.ttl   = leg.ttl,
    leg.ttl.size = leg.ttl.size,
    facet.bg  = facet.bg,
    mode      = if(dark) "dark" else "light"
  ) 
  
  if(figplot){
    # Warn if the user specified a non-classic theme
    if(!missing(theme.style) && theme.style != "classic"){
      warning(sprintf(
        "figplot = TRUE ignores custom themes (theme.style = '%s'). Use figplot = FALSE for full theming.",
        theme.style
      ))
    }
    
    # Apply plot_theme for figplot figure
    plt <- plt &
      do.call(plot_theme, c(
        theme_args,
        list(
          x.ttl = FALSE,
          ticks = FALSE,
          line = FALSE,
          border = FALSE,
          grid.major = FALSE,
          grid.minor = FALSE,
          panel.fill = "white"
        ),
        list(...)
      ))    
    text_col <- "black"
  } else {
    # Apply plot_theme normally
    plt <- plt & do.call(plot_theme, c(theme_args, list(...)))
  }
  
  
  # Add cluster labels if requested
  if(label){
    umap_data <- dplyr::tibble(x=emb[, dims[1]], y=emb[, dims[2]], cluster=as.character(object@active.ident)) %>%
      dplyr::group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), .groups="drop")
    plt <- plt + ggrepel::geom_text_repel(
      data=umap_data, aes(x, y, label=cluster),
      color = if(dark) "white" else "black",
      fontface = label.face,
      bg.color = if(dark) "#3A3A3A" else "grey95",
      bg.r = 0.1, size = label.size, seed = 42
    )
  }
  
  # figplot arrow axes (minimal figure)
  if(figplot){
    x.lab.reduc <- plt$labels$x %||% paste0(toupper(reduction), dims[1])
    y.lab.reduc <- plt$labels$y %||% paste0(toupper(reduction), dims[2])
    plt <- plt & Seurat::NoAxes()
    L <- 0.12
    axis.df <- data.frame(x0=c(0,0), y0=c(0,0), x1=c(L,0), y1=c(0,L))
    axis.plot <- ggplot2::ggplot(axis.df) +
      ggplot2::geom_segment(ggplot2::aes(x=x0, y=y0, xend=x1, yend=y1), linewidth=0.4, lineend="round") +
      ggplot2::xlab(x.lab.reduc) + ggplot2::ylab(y.lab.reduc) +
      ggplot2::coord_fixed() + ggplot2::theme_classic(base_size=font.size) +
      ggplot2::theme(plot.background=ggplot2::element_rect(fill="transparent", colour=NA),
                     panel.background=ggplot2::element_rect(fill="transparent", colour=NA),
                     axis.text=ggplot2::element_blank(),
                     axis.ticks=ggplot2::element_blank(),
                     axis.line=ggplot2::element_blank(),
                     panel.border=ggplot2::element_blank(),
                     axis.title=ggplot2::element_text(size=font.size, face="plain"),
                     plot.margin=ggplot2::margin(0,0,0,0))
    figure.layout <- c(patchwork::area(t=1,l=1,b=11,r=11), patchwork::area(t=10,l=1,b=11,r=2))
    return(plt + axis.plot + patchwork::plot_layout(design=figure.layout))
  }
  
  
  if(!legend) plt <- plt & Seurat::NoLegend()
  if(no.axes) plt <- plt & Seurat::NoAxes()
  
  plt
}




#' Cell Dot Plot (Enhanced Seurat DotPlot)
#'
#' A cleaner, more customizable wrapper around **Seurat::DotPlot**, providing
#' improved color handling, optional dot outlines, flexible axis formatting,
#' legend placement, and theme control. Useful for visualizing gene expression
#' patterns across clusters or metadata-defined groups.
#'
#' @param object A Seurat object.
#' @param features Character vector of features (genes or metadata fields) to plot.
#' @param group.by Column in `object@meta.data` used to group cells.
#'   Default: `"seurat_clusters"`.
#'
#' @param th.cols Color palette name from **RColorBrewer** used for the
#'   expression gradient. Default: `"Reds"`.
#' @param rev.th.cols Logical; reverse the gradient palette. Default: FALSE.
#'
#' @param dot.scale Numeric scale factor controlling the dot size range.
#'   Passed to `Seurat::DotPlot`. Default: 4.5.
#' @param dot.outline Logical; draw outlines around dots. Default: FALSE.
#'
#' @param x.angle Angle for x-axis labels (degrees). Default: 90.
#' @param vjust.x,hjust.x Vertical and horizontal justification for x labels.
#'
#' @param flip Logical; swap x and y axes using `coord_flip()`. Default: FALSE.
#'
#' @param font.size Base font size passed to internal theme helper. Default: 8.
#'
#' @param plot.title Optional plot title.
#'
#' @param leg.size Legend text size. Default: 8.
#' @param leg.pos Position of the legend (e.g., `"right"`, `"bottom"`).
#' @param leg.just Legend justification.
#' @param leg.hjust Logical; if TRUE, use a horizontal legend layout when possible.
#'
#' @param x.axis.pos Position of the x-axis (`"top"` or `"bottom"`).
#' @param theme ggplot2 theme name used by the internal theme helper.
#'
#' @param x.face,y.face Logical; italic styling for x and/or y-axis labels.
#' @param x.ttl,y.ttl Logical; italic styling for x and/or y-axis titles.
#'
#' @param ... Additional parameters passed to `Seurat::DotPlot()`.
#'
#' @details
#' This function enhances the standard Seurat dot plot by providing:
#' * Customizable Brewer color gradients  
#' * Optional dot outlines  
#' * Flexible axis label styling  
#' * Improved legend customization and ordering  
#' * Optional axis flipping  
#'
#' It retains all functionality of `Seurat::DotPlot` while adding cleaner,
#' publication-ready defaults.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' celldot(pbmc, features = c("MS4A1","CD3D"))
#'
#' celldot(pbmc, features = c("MS4A1","CD14"), th.cols = "Blues",
#'         dot.outline = TRUE, flip = TRUE)
#' }
#'
#' @export
#' 
celldot <- function(object, features, group.by="seurat_clusters", th.cols="Reds",
                    rev.th.cols=FALSE, dot.scale=4.5, x.angle=90, vjust.x=NULL,
                    hjust.x=NULL, flip=FALSE, font.size=8, plot.title=NULL,
                    leg.size=10, leg.pos="right", leg.just="bottom", leg.hjust=FALSE,
                    x.axis.pos="bottom", theme="classic", x.face=FALSE, y.face=FALSE,
                    x.ttl=FALSE, y.ttl=FALSE, dot.outline=FALSE, ...) {
  
  stopifnot(inherits(object,"Seurat"))
  object <- Seurat::SetIdent(object, value=group.by)
  features <- unique(features)
  
  pal <- RColorBrewer::brewer.pal(9, th.cols)
  if (rev.th.cols) pal <- rev(pal)
  
  outline_col <- if (dot.outline) "gray60" else NA
  outline_stroke <- if (dot.outline) 0.5 else 0
  
  plt <- suppressWarnings({
    suppressMessages({
      Seurat::DotPlot(object, features=features, dot.scale=dot.scale, ...) 
    })
  }) +
    scale_color_gradientn(colors=pal, oob=scales::squish) +
    geom_point(aes(size=pct.exp), shape=21, colour=outline_col, stroke=outline_stroke) +
    labs(title=plot.title, color="Average\nExpression", size="Percent\nExpressed") +
    plot_theme(theme=theme, font.size=font.size, x.angle=x.angle,
               x.hjust=hjust.x, x.vjust=vjust.x, xy.val=TRUE, x.lab=TRUE, y.lab=TRUE, ...) +
    theme(
      axis.text.x = if (x.face || (flip && y.face)) element_text(face="italic") else element_text(),
      axis.text.y = if (y.face || (flip && x.face)) element_text(face="italic") else element_text(),
      axis.title = element_blank(), 
      legend.spacing.y = unit(0.05, "cm"),
      legend.spacing.x = unit(0.05, "cm"),
      legend.box.spacing = unit(0.05, "cm"),
      legend.margin = margin(2,2,2,2)
    )
  
  if (flip) plt <- plt + coord_flip()
  if (is.list(features)) plt <- plt + theme(strip.text.x=element_text(angle=45))
  
  # Legend positioning
  # Legend positioning outside plot, bottom-right
  if (!is.null(leg.pos)) {
    if (leg.pos == "right") {
      plt <- plt + theme(
        legend.position = "right",
        legend.justification = c("right","bottom"),
        legend.box.just = "right",
        legend.box.margin = margin(0,0,0,0)
      )
    } else if (leg.pos == "left") {
      plt <- plt + theme(
        legend.position = "left",
        legend.justification = c("left","bottom"),
        legend.box.just = "left",
        legend.box.margin = margin(0,0,0,0)
      )
    } else if (leg.pos == "top") {
      plt <- plt + theme(
        legend.position = "top",
        legend.justification = c("right","top"),
        legend.box.just = "right"
      )
    } else if (leg.pos == "bottom") {
      plt <- plt + theme(
        legend.position = "bottom",
        legend.justification = c("right","bottom"),
        legend.box.just = "right"
      )
    }
  }
  
  # Keep original guides for color and size
  guide_color <- guide_colorbar(frame.colour="black", ticks.colour="black")
  guide_size <- guide_legend(override.aes=list(shape=21, colour=outline_col, fill="black"))
  guide_color$order <- 1
  guide_size$order  <- 2
  plt <- plt + guides(color=guide_color, size=guide_size)
  
  plt
}


#' Cell Feature Violin Plot with Statistics
#'
#' Plots expression or metadata features as violin plots for a Seurat object,
#' optionally adding median points, shared y-axis scaling, flipped axes, and
#' statistical comparisons (Wilcoxon for 2 groups, Kruskal-Wallis for >2 groups).
#'
#' @param obj A Seurat object.
#' @param features Character vector of feature names (genes or metadata columns) to plot.
#' @param ncol Number of columns in the output patchwork plot. Defaults to sqrt(#features / 1.5).
#' @param stack Logical; if TRUE, plots are stacked in a single column.
#' @param shared.y Logical; if TRUE, all violins share the same y-axis.
#' @param ttl.pos Position of subplot titles: "center", "left", or "right".
#' @param group.by Metadata column to group by. Defaults to "seurat_clusters".
#' @param split.by Optional metadata column to split violins by.
#' @param assay Assay to pull data from. Default is "RNA".
#' @param slot Slot to use for expression values. One of "data", "counts", or "scale.data".
#' @param log Logical; if TRUE, log-transform the expression values.
#' @param cols Optional named vector of colors for each group. If NULL, defaults are used.
#' @param med Logical; if TRUE, overlay median points on each violin.
#' @param med.size Size of median points if med = TRUE.
#' @param pt.size Size of jittered points. Set to 0 to hide points.
#' @param border.size Size of the violin border lines.
#' @param txtsize Base font size for titles and labels.
#' @param theme ggplot2 theme to use: "classic", "minimal", etc.
#' @param x.ang Rotation angle of x-axis labels.
#' @param leg.pos Position of legend: "none", "right", "left", etc.
#' @param title Optional overall title for the patchwork plot.
#' @param rm.subtitles Logical; if TRUE, removes individual subplot titles.
#' @param flip Logical; if TRUE, flips x and y axes.
#' @param auto.resize Logical; if TRUE, sets dynamic width/height attributes.
#' @param ylab.global Global y-axis label. Defaults to expression level.
#' @param xlab.global Global x-axis label. Defaults to blank.
#' @param pairwise Logical; if TRUE, perform pairwise comparisons between groups.
#' @param add.stats Logical; if TRUE, add p-values to plots.
#' @param show.pval Logical; if TRUE, show p-values above violins.
#' @param pval.label Character; label type for p-values, e.g., "p.signif" or "p.format".
#' @param ... Additional arguments passed to ggplot2 layers.
#'
#' @return A patchwork object containing the violin plots.
#' @examples
#' \dontrun{
#' cellvio(sub, features = c("MEG3","TP63","HES6"),
#'         group.by = "ann_level_2",
#'         pt.size = 0.1, ncol = 3, pairwise = TRUE,
#'         font.size = 10, show.pval = TRUE)
#' }
#' @export
#' 
cellvio <- function(
    obj, features,
    ncol = NULL,
    shared.y = FALSE,
    ttl.pos = c("center", "left", "right"),
    group.by = "seurat_clusters",
    split.by = NULL,
    stack = FALSE,
    assay = "RNA",
    slot = "data",
    log = FALSE,
    cols = NULL,
    med = FALSE,
    med.size = 1,
    pt.size = 0,
    border.size = 0.1,
    style = "classic",
    leg.pos = "none",
    x.ang = 45,
    title = NULL,
    rm.subttl = FALSE,
    flip = FALSE,
    auto.resize = TRUE,
    ylab.global = NULL,
    xlab.global = NULL,
    add.stats = FALSE,
    show.pval = FALSE,
    pairwise = FALSE,
    pval.label = "p.signif",
    txtsize = 10,
    ...
) {
  
  stopifnot(inherits(obj, "Seurat"))
  if (length(features) == 0) stop("features must be provided.")
  ttl.pos <- match.arg(ttl.pos)
  
  
  # determine ncol
  if (is.null(ncol)) {
    ncol <- if (stack) 1 else max(1, ceiling(sqrt(length(features) / 1.5)))
  }
  
  # Set identities if grouping
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(obj@meta.data))
      stop(paste(group.by, "not found in metadata."))
    Idents(obj) <- group.by
  }
  
  # Determine which features exist
  f.expr <- intersect(features, rownames(obj[[assay]]))
  f.meta <- intersect(features, colnames(obj@meta.data))
  features <- unique(c(f.expr, f.meta))
  if (!length(features)) stop("No features found in assay or metadata.")
  
  # Shared y-scale
  ymax <- NULL
  if (shared.y) {
    vals <- c(
      if (length(f.expr)) as.numeric(Seurat::GetAssayData(obj, assay, slot)[f.expr, ]),
      if (length(f.meta)) as.numeric(as.matrix(obj@meta.data[, f.meta, drop = FALSE]))
    )
    vals <- vals[is.finite(vals)]
    if (length(vals)) ymax <- max(vals)
  }
  
  # Colors
  if (is.null(cols)) {
    g <- tryCatch(unique(obj[[group.by]][, 1]), error = \(e) NULL)
    cols <- custom_palette(length(g))
    #cols <- scales::hue_pal()(if (is.null(g)) 8 else length(g))
    #cols <- ggpubr::get_palette("npg", length(g))
    names(cols) <- g
  }
  
  # Helper to build a single violin with stats
  vln <- function(f) {
    
    if (log && slot == "counts") {
      obj[[assay]]@data[f, ] <- log1p(obj[[assay]]@counts[f, ])
      slot_use <- "data"
    } else {
      slot_use <- slot
    }
    
    df <- Seurat::FetchData(obj, vars = c(group.by, f))
    names(df)[2] <- "value"
    df[[group.by]] <- factor(df[[group.by]]) # ensure factor
    
    # Base plot
    # p <- ggplot(df, aes_string(group.by, "value", fill = group.by)) +
    #   geom_violin(scale = "width", color = "black", size = border.size) +
    #   scale_fill_manual(values = cols)
    
    p <- suppressWarnings({
      suppressMessages({Seurat::VlnPlot(
        obj,
        features = f,
        group.by = group.by,
        split.by = split.by,
        assay = assay,
        slot = slot,
        pt.size = pt.size,
        cols = cols,
        ...
      ) + scale_y_continuous(
        expand = expansion(mult = c(0.05, 0.25))
      )
      })
    })
    
    # Points
    #if (pt.size > 0) p <- p + geom_jitter(width = 0.1, size = pt.size, alpha = 0.6)
    
    # Add statistics
    
    if (add.stats && show.pval) {
      if (nlevels(df[[group.by]]) > 1) {

        df$.grp <- df[[group.by]]
        if (!is.null(split.by)) {
          df$.grp <- interaction(df[[group.by]], obj@meta.data[[split.by]], drop = TRUE)
        }

        y_max <- max(df$value, na.rm = TRUE)
        y_step <- (ymax %||% y_max) * 0.08

        if (pairwise && nlevels(df$.grp) > 1) {

          cmp <- utils::combn(levels(df$.grp), 2, simplify = FALSE)
          stat_df <- ggpubr::compare_means(value ~ .grp, data = df, method = "wilcox.test", comparisons = cmp)
          stat_df$y.position <- y_max + seq_len(nrow(stat_df)) * y_step

          p <- p + ggpubr::stat_pvalue_manual(
            stat_df,
            label = pval.label,
            y.position = "y.position",
            tip.length = 0.02,
            size = txtsize * 0.25
          )

        } else {

          # Kruskal test if not pairwise
          stat_df <- ggpubr::compare_means(value ~ .grp, data = df, method = "kruskal.test")
          p <- p + annotate(
            "text",
            x = 1,
            y = y_max + y_step,
            label = paste0(signif(stat_df$p, 3)),
            hjust = 0,
            size = txtsize * 0.25
          )

        } 
      } 
    } 
    
    
    # Titles
    p <- if (!rm.subttl) p + labs(title = f) else p + labs(title = NULL)
    
    # Make feature titles bold + italic
    p <- p + theme(
      plot.title = element_text(face = "bold.italic")
    )
    
    .style_layers <- function(
    p,
    violin_lw = 0.15,
    point_size = NULL,
    jitter_width = NULL
    ) {
      for (i in seq_along(p$layers)) {
        
        layer <- p$layers[[i]]
        
        # Violin outline 
        if (inherits(layer$geom, "GeomViolin")) {
          layer$aes_params$linewidth <- violin_lw
        }
        
        # Points (Seurat uses GeomPoint + position_jitterdodge) 
        if (inherits(layer$geom, "GeomPoint")) {
          
          if (!is.null(point_size)) {
            layer$aes_params$size  <- point_size
            layer$aes_params$alpha <- 0.6
          }
          
          if (!is.null(jitter_width) &&
              inherits(layer$position, "PositionJitterdodge")) {
            layer$position$width <- jitter_width
          }
        }
        
        p$layers[[i]] <- layer
      }
      p
    }
    
    p <- .style_layers(
      p,
      violin_lw  = border.size,
      point_size = if (pt.size > 0) pt.size else NULL,
      jitter_width = 0.08
    )
    
    # Theme & formatting
    p <- p + plot_theme(style = style, txtsize = txtsize, x.ang = x.ang,
                        leg.pos = leg.pos, x.ttl = FALSE, ttl.pos = ttl.pos,...) +
      theme(
        plot.title = element_text(face = "bold.italic")
      )
    
    # Force final legend position
    if (!is.null(leg.pos)) {
      p <- p + theme(legend.position = leg.pos)
    }
    
    if (med) p <- p + stat_summary(fun = median, geom = "point", shape = 3, size = med.size)
    if (!is.null(ymax)) p <- p + ylim(0, ymax)
    if (flip) p <- p + coord_flip()
    p + ylab(NULL)
  }
  
  # number of cols
  if (is.null(ncol))
    ncol <- max(1, ceiling(sqrt(length(features) / 1.5)))
  
  plist <- lapply(features, vln)
  total <- length(plist)
  
  # Only show x-axis on bottom plots
  bottom <- sapply(1:ncol, \(i) max(seq(i, total, by = ncol)))
  for (i in seq_along(plist)) {
    if (!(i %in% bottom)) {
      plist[[i]] <- plist[[i]] +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    }
  }
  
  # Layout
  combo <- patchwork::wrap_plots(plist, ncol = ncol) +
    patchwork::plot_layout(guides = "collect")
  
  if (!is.null(title)) {
    combo <- combo +
      patchwork::plot_annotation(
        title = title,
        theme = theme(title = element_text(face = "bold"))
      )
  }
  
  # Global labels
  auto_y <- switch(slot,
                   data = "Expression level",
                   counts = "Raw counts",
                   scale.data = "Scaled expression",
                   "Expression level")
  
  ylab <- if (is.null(ylab.global)) auto_y else ylab.global
  xlab <- if (is.null(xlab.global)) "" else xlab.global
  
  plt <- cowplot::ggdraw(combo) +
    cowplot::draw_label(ylab, x = -0.01, y = 0.55, angle = 90, size = txtsize) +
    cowplot::draw_label(xlab, x = 0.5, y = 0.02, size = txtsize) +
    theme(plot.margin = margin(15, 15, 15, 15))
  
  # Auto resize attributes
  if (auto.resize) {
    ng <- length(unique(obj[[group.by]][, 1]))
    attr(plt, "dynamic_width") <- 6 + ng * 0.3
    attr(plt, "dynamic_height") <- 4 + length(features) * 0.25
  }
  
  plt
}


#' Plot Unique Gene Counts Per Cell (Robust Version)
#'
#' This function visualizes the number of selected genes expressed per cell.
#' It accepts either a Seurat object with a gene list or a precomputed per-cell table.
#' Missing or unexpressed genes are automatically handled with warnings.
#'
#' @param object Either a Seurat object or a data.frame with columns `Cell` and `Unique`.
#' @param gene.list Character vector of genes to evaluate (required if `object` is a Seurat object).
#' @param plot.type One of "bar", "hist", or "violin".
#' @param font.size Numeric font size.
#' @param theme Theme type passed to `plot_theme()`.
#' @param x.lab X-axis title.
#' @param y.lab Y-axis title.
#' @param ... Additional arguments passed to `plot_theme()`.
#'
#' @return A ggplot2 object.
#' @export
plot_unique_gene_counts <- function(
    object,
    gene.list = NULL,
    plot.type = c("bar", "hist", "violin"),
    font.size = 8,
    theme = "classic",
    x.lab = "Number of cells",
    y.lab = "Number of ORs",
    color = "#20679B",
    ...
) {
  
  # --- Get per-cell table ---
  if (inherits(object, "Seurat")) {
    
    if (is.null(gene.list)) stop("If `object` is a Seurat object, you must supply `gene.list`.")
    
    counts <- rownames(object[["RNA"]]@counts)
    genes.present <- intersect(gene.list, counts)
    
    if (length(genes.present) == 0) stop("None of the genes in `gene.list` exist in the object.")
    
    if (length(genes.present) < length(gene.list)) {
      warning(sprintf("Only %d/%d genes found in the object. Proceeding with available genes.", 
                      length(genes.present), length(gene.list)))
    }
    
    tbl <- get_unique_gene_table(object, genes.present)
    per.cell <- as.data.frame(tbl$per.cell)
    
  } else {
    per.cell <- as.data.frame(object)
    if (!all(c("Cell", "Unique") %in% colnames(per.cell))) {
      stop("`object` must be a Seurat object OR a data.frame with columns: Cell, Unique")
    }
  }
  
  # --- Filter non-expressing cells ---
  df <- per.cell[per.cell$Unique > 0, , drop = FALSE]
  
  if (nrow(df) == 0) stop("No cells express the selected genes.")
  
  df$Unique.factor <- factor(df$Unique)
  
  plot.type <- match.arg(plot.type)
  
  # --- Generate plot ---
  plt <- switch(
    plot.type,
    "bar" = ggplot(df, aes(x = Cell, y = Unique.factor)) +
      geom_bar(stat = "identity", colour = color) +
      labs(x = x.lab, y = y.lab) +
      plot_theme(theme = theme, legend.position = "none", font.size = font.size) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
    
    "hist" = ggplot(df, aes(x = Unique)) +
      geom_histogram(binwidth = 1, fill = color, colour = "black") +
      labs(x = y.lab, y = "Number of cells") +
      plot_theme(theme = theme, font.size = font.size),
    
    "violin" = ggplot(df, aes(x = "", y = Unique)) +
      geom_violin(trim = FALSE, fill = color) +
      geom_jitter(width = 0.1, alpha = 0.4) +
      labs(x = "", y = y.lab) +
      plot_theme(theme = theme, font.size = font.size, ...)
  )
  
  return(plt)
}


cellpct <- function(
    object,
    cell.col = "ann2",
    group.col = "sex",
    donor.col = "sample",
    xttl = NULL,
    yttl = "Cell type fraction per donor (log-shifted)",
    pseudo = 1e-3,
    jitter.size = 1.5,
    alpha = 0.85,
    plot_theme_fn = plot_theme,
    group.colors = NULL,
    donor.colors = NULL,
    font.size = 10,
    leg.pos = "right",
    x.angle = 45,
    theme.style = "classic",
    leg.ncol = 1,
    show.pval = TRUE,
    ...
) {
  
  library(dplyr)
  library(ggplot2)
  library(rlang)
  
  # 1. Compute donor-level fractions
  frac_df <- object@meta.data %>%
    group_by(!!sym(donor.col), !!sym(group.col), !!sym(cell.col)) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(!!sym(donor.col)) %>%
    mutate(
      total_cells = sum(n_cells),
      fraction = n_cells / total_cells,
      frac_log = log10(fraction + pseudo)
    ) %>%
    ungroup()
  
  # 2. Shift log values
  min_val <- min(frac_df$frac_log, na.rm = TRUE)
  frac_df <- frac_df %>%
    mutate(frac_shift = frac_log - min_val)
  
  # 3. Statistics per cell type
  stat_df <- frac_df %>%
    distinct(!!sym(donor.col), !!sym(group.col), !!sym(cell.col), fraction, frac_shift) %>%
    group_by(!!sym(cell.col)) %>%
    summarise(
      n_group1 = sum(!!sym(group.col) == levels(factor(!!sym(group.col)))[1]),
      n_group2 = sum(!!sym(group.col) == levels(factor(!!sym(group.col)))[2]),
      p_value = if (n_group1 >= 2 && n_group2 >= 2)
        wilcox.test(fraction ~ !!sym(group.col), exact = FALSE)$p.value
      else NA_real_,
      max_val = max(frac_shift),
      .groups = "drop"
    ) %>%
    mutate(
      sig = case_when(
        is.na(p_value)  ~ "ns",
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ "ns"
      ),
      y_pos = max_val * 1.08
    )
  
  # 3b. Merge stats into source data
  source_df <- frac_df %>%
    left_join(stat_df %>% dplyr::select(!!sym(cell.col), p_value, sig),by = cell.col)
  
  # 4. Donor colors
  if (!is.null(donor.colors)) {
    frac_df[[donor.col]] <- factor(frac_df[[donor.col]], levels = names(donor.colors))
  } else {
    donors <- unique(frac_df[[donor.col]])
    donor.colors <- setNames(scales::hue_pal()(length(donors)), donors)
  }
  
  # 5. Base plot
  p <- ggplot(frac_df, aes(x = !!sym(cell.col), y = frac_shift, fill = !!sym(group.col))) +
    stat_boxplot(geom = "errorbar", width = 0.8, position = position_dodge(0.8)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
    geom_jitter( aes(color = !!sym(donor.col)),size = jitter.size,alpha = alpha,
                 position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
    labs(x = xttl, y = yttl, fill = group.col) +
    plot_theme_fn(font.size = font.size,x.angle = x.angle,theme.style = theme.style,leg.pos = leg.pos,...)
  
  # 6. Add significance labels
  if (show.pval) {
    p <- p + geom_text(data = stat_df,aes(x = !!sym(cell.col), y = y_pos, label = sig),
                       inherit.aes = FALSE,size = 4)
  }
  
  # 7. Custom colors
  if (!is.null(group.colors)) p <- p + scale_fill_manual(values = group.colors)
  if (!is.null(donor.colors)) p <- p + scale_color_manual(values = donor.colors)
  
  # 8. Guides (number of columns)
  if (!is.null(leg.ncol)) {
    p <- p + guides(
      color = guide_legend(ncol = leg.ncol, override.aes = list(size = 3)),
      fill  = guide_legend(ncol = leg.ncol)
    )
  }
  
  # 9. Force legend position after all guides/scales
  if (!is.null(leg.pos)) {p <- p + theme(legend.position = leg.pos)}
  
  return(list(plot = p, source_data = source_df))
}



#' Plot Counts or Proportions by Group
#'
#' This function generates a ggplot2 bar, point, or box plot showing counts or
#' proportions of a variable (y) across groups (x). It supports Seurat objects
#' (uses `meta.data`), custom color palettes, stacking, coordinate flipping, and facetting.
#'
#' @param data Data frame or Seurat object (if Seurat, uses meta.data)
#' @param x Grouping variable (bare name)
#' @param y Identity variable (bare name)
#' @param plot.type Plot type: "bar", "point", or "box" (default "count")
#' @param prop Logical; if TRUE, plot proportions instead of counts
#' @param prop.multi Logical; compute proportions per group for multiple facets
#' @param stack Logical; stack bars (for bar plot) instead of dodging
#' @param coord.flip Logical; flip x and y axes
#' @param colors Vector of colors; if NULL, defaults to hue palette or RColorBrewer
#' @param use.brewer Logical; use RColorBrewer palette
#' @param brew.pal Brewer palette name (default "Set1")
#' @param raster Logical; not implemented (reserved)
#' @param theme ggplot2 theme type (default "classic")
#' @param font.size Base font size for text
#' @param x.angle Rotation angle for x-axis labels
#' @param ncol Number of columns for facet_wrap (if prop.multi)
#' @param legend Logical; whether to show legend
#' @param legend.title Legend title; if NULL, uses y variable
#' @param legend.text.size Legend text size
#' @param legend.position Legend position: "right", "bottom", etc.
#' @param legend.ncol Number of columns for legend
#' @param x.title Custom x-axis title
#' @param y.title Custom y-axis title
#' @param x.lab Logical; show x-axis labels
#' @param xy.lab Logical; show both x and y axis labels
#' @param show.contour Logical; add black border around bars (bar plot)
#'
#' @return ggplot object
#' @examples
#' \dontrun{
#' cellprop(df, x = dominance.bin, y = cluster, plot.type = "bar", prop = TRUE)
#' cellprop(sub, x = cell_group, y = dominance_bin, prop.multi = TRUE)
#' }
#' @export
cellprop <- function(
    data,
    x, y,
    plot.type = "bar",
    prop = FALSE,
    prop.multi = FALSE,
    percent.stack = FALSE,
    stack = FALSE,
    coord.flip = FALSE,
    x.reverse = FALSE,
    y.reverse = FALSE,
    colors = NULL,
    use.brewer = FALSE,
    brew.pal = "Set1",
    raster = FALSE,
    theme = "classic",
    font.size = 8,
    x.angle = 90,
    ncol = NULL,
    legend = TRUE,
    legend.title = NULL,
    legend.text.size = 8,
    legend.position = "right",
    legend.ncol = 1,
    x.title = NULL,
    y.title = NULL,
    x.lab = TRUE,
    xy.lab = TRUE,
    show.contour = TRUE,
    add.pval = FALSE,
    pval.test = "chisq",
    pval.size = 3,
    ...
) {
  require(ggplot2)
  require(dplyr)
  require(RColorBrewer)
  
  # Handle Seurat
  if (inherits(data, "Seurat")) data <- data@meta.data
  
  x_var <- rlang::enquo(x)
  y_var <- rlang::enquo(y)
  
  # Count table
  stat <- data %>%
    dplyr::select(!!x_var, !!y_var) %>%
    dplyr::rename(.group = !!x_var, .ident = !!y_var) %>%
    dplyr::group_by(.group, .ident) %>%
    dplyr::summarise(.n = n(), .groups = "drop")
  
  # proportions (per x-group)
  if (prop || prop.multi || percent.stack) {
    stat <- stat %>%
      dplyr::group_by(.group) %>%
      dplyr::mutate(.value = .n / sum(.n) * 100) %>%
      dplyr::ungroup()
  } else {
    stat$.value <- stat$.n
  }
  
  # reverse x-axis
  if (x.reverse) {
    stat$.group <- factor(stat$.group, levels = rev(sort(unique(stat$.group))))
  } else {
    stat$.group <- factor(stat$.group, levels = sort(unique(stat$.group)))
  }
  
  stat$.ident <- factor(stat$.ident)
  
  # colors
  n.colors <- length(unique(stat$.ident))
  if (is.null(colors)) {
    if (use.brewer) {
      colors <- colorRampPalette(brewer.pal(min(9, n.colors), brew.pal))(n.colors)
    } else {
      colors <- scales::hue_pal()(n.colors)
    }
  }
  
  # Build main plot
  p <- ggplot(stat, aes(x = .group, y = .value, fill = .ident))
  
  # stacked or dodged bars
  if (plot.type == "bar") {
    p <- p + geom_bar(
      stat = "identity",
      position = if (stack || percent.stack) "stack" else "dodge",
      color = if (show.contour) "black" else NA,
      linewidth = 0.2
    )
  }
  
  # point/box optional
  if (plot.type == "point") p <- p + geom_point(size = 2)
  if (plot.type == "box") p <- p + geom_boxplot()
  
  # percent stacked bar → fix y-scale to 100%
  if (percent.stack) {
    p <- p + scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0, 0.05)))
    y.title <- "Percent (%)"
  }
  
  # prop.multi facet mode
  if (prop.multi) {
    p <- p +
      facet_wrap(~.ident, scales = "free_y", ncol = ncol) +
      scale_y_continuous(labels = function(x) paste0(x, "%"))
  }
  
  # axis reversing
  if (coord.flip) p <- p + coord_flip()
  if (y.reverse && !percent.stack && !prop.multi) p <- p + scale_y_reverse()
  
  # labels + theme
  p <- p +
    scale_fill_manual(values = colors) +
    plot_theme(theme = theme, font.size = font.size, x.angle = x.angle, x.lab = x.lab, xy.lab = xy.lab, ...) +
    labs(
      x = if (!is.null(x.title)) x.title else rlang::as_name(x_var),
      y = if (!is.null(y.title)) y.title else if (prop || prop.multi || percent.stack) "Percent (%)" else "Count",
      fill = if (!is.null(legend.title)) legend.title else rlang::as_name(y_var)
    ) +
    theme(
      legend.title = element_text(size = legend.text.size),
      legend.text = element_text(size = legend.text.size),
      legend.position = if (legend) legend.position else "none"
    ) +
    guides(fill = guide_legend(ncol = legend.ncol))
  
  # ADD P-VALUES
  if (add.pval) {
    pval.df <- stat %>%
      tidyr::pivot_wider(names_from = .ident, values_from = .n, values_fill = 0) %>%
      dplyr::rowwise() %>%
      mutate(
        pval =
          if (pval.test == "chisq")
            chisq.test(c_across(!.group))$p.value
        else if (pval.test == "fisher")
          fisher.test(matrix(c_across(!.group), nrow = 1))$p.value
        else NA_real_
      ) %>%
      ungroup() %>%
      mutate(
        label = paste0("p=", signif(pval, 2)),
        y.position = max(stat$.value) * 1.05
      )
    
    p <- p +
      geom_text(
        data = pval.df,
        aes(x = .group, y = y.position, label = label),
        inherit.aes = FALSE,
        size = pval.size
      )
  }
  
  return(p)
}


##' Pseudobulk cell-type proportion plotting with statistics
#'
#' Generate pseudobulk cell-type proportion plots from single-cell data,
#' supporting Seurat objects or precomputed data frames. The function
#' computes donor-level proportions, summarizes data as mean ± SEM,
#' performs statistical testing, and returns both the plot and
#' Nature-style source data.
#'
#' @param input A \code{Seurat} object or a data frame containing
#'   donor-, stage-, and cell-level information.
#' @param donor.col Column name identifying biological replicates
#'   (e.g. patient, sample, donor).
#' @param stage.col Column name identifying experimental groups or stages.
#' @param cell.col Column name identifying cell types or clusters.
#' @param prop.col Optional column name containing precomputed proportions.
#'   If NULL, proportions are calculated automatically.
#' @param plot.type Plot type: \code{"bar"} (mean ± SEM with points)
#'   or \code{"box"} (distribution with jittered points).
#' @param min.lmm Minimum number of donors per group required to apply
#'   linear mixed models. Below this threshold, Kruskal–Wallis is used.
#' @param facet.by Variable used for faceting (typically \code{"cell"}).
#' @param cell.order Optional character vector specifying cell-type order.
#' @param palette Brewer palette name used when \code{stage.cols} is NULL.
#' @param donor.cols Optional named vector of colors for donors.
#' @param stage.cols Optional named vector of colors for stages/groups.
#' @param show.p Logical; whether to display p-values on the plot.
#' @param title Plot title.
#' @param x.lab X-axis label.
#' @param y.lab Y-axis label.
#' @param leg.pos Legend position passed to \code{plot_theme()}.
#' @param pt.size Point size for jittered donor-level points.
#' @param ncol Number of columns for facet wrapping.
#' @param ... Additional arguments passed directly to \code{plot_theme()}.
#'
#' @details
#' For Seurat objects, cell-type proportions are computed per donor
#' (pseudobulk). Summary statistics are calculated as mean ± SEM
#' across donors. Statistical testing is performed independently
#' for each cell type:
#'
#' \itemize{
#'   \item Linear mixed model (LMM): \code{prop ~ stage + (1 | donor)}
#'   \item Kruskal–Wallis test when donor numbers are insufficient
#' }
#'
#' The returned \code{source_data} table is suitable for direct inclusion
#' as Nature-style Source Data and includes raw proportions, donor counts,
#' summary statistics, and p-values.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{plot}{A \code{ggplot2} object.}
#'   \item{source_data}{A data frame containing raw and summarized data,
#'   statistical results, and donor counts.}
#' }
#'
#' @examples
#' res <- cellpbulk(
#'   input = seurat_obj,
#'   donor.col = "sample",
#'   stage.col = "condition",
#'   cell.col  = "celltype",
#'   plot.type = "bar",
#'   show.p = TRUE
#' )
#'
#' res$plot
#' head(res$source_data)
#'
#' @export
cellpbulk <- function(
    input,
    donor.col = "sample",
    stage.col = "stage",
    cell.col  = "ann2",
    prop.col  = NULL,
    plot.type = c("bar", "box"),
    min.lmm   = 3,
    facet.by  = "cell",
    cell.order = NULL,
    palette = "Set3",
    donor.cols = NULL,
    stage.cols = NULL,
    show.p = TRUE,
    title = NULL,
    x.lab = NULL,
    y.lab = "Proportion (%)",
    leg.pos = "right",
    pt.size = 2,
    label.size = 3,
    ncol = NULL,
    theme.style = "classic",
    facet.title.size = 10,
    ...
) {
  library(dplyr)
  library(ggplot2)
  library(lme4)
  
  plot.type <- match.arg(plot.type)
  
  # 1. Build pseudobulk dataframe
  if (inherits(input, "Seurat")) {
    meta <- input@meta.data %>% as.data.frame()
    
    df <- meta %>%
      mutate(
        donor = .data[[donor.col]],
        stage = .data[[stage.col]],
        cell  = .data[[cell.col]]
      ) %>%
      group_by(donor, stage, cell) %>%
      summarise(n_cells = n(), .groups = "drop") %>%
      group_by(donor, stage) %>%
      mutate(prop = n_cells / sum(n_cells)) %>%
      ungroup()
  } else {
    df <- input %>%
      rename(
        donor = all_of(donor.col),
        stage = all_of(stage.col),
        cell  = all_of(cell.col)
      )
    
    if (!is.null(prop.col)) {
      df <- df %>% rename(prop = all_of(prop.col))
    } else if (!"prop" %in% colnames(df)) {
      stop("prop.col must be provided or 'prop' column must exist")
    }
  }
  
  # 2. Factor ordering
  if (!is.null(cell.order)) {
    df$cell <- factor(df$cell, levels = intersect(cell.order, unique(df$cell)))
  } else {
    df$cell <- factor(df$cell)
  }
  
  # 3. Donor counts + Mean ± SEM
  donor_n <- df %>%
    group_by(cell, stage) %>%
    summarise(n_donor = n_distinct(donor), .groups = "drop")
  
  summary_stats <- df %>%
    group_by(cell, stage) %>%
    summarise(
      n = n_distinct(donor),
      mean = mean(prop, na.rm = TRUE),
      sd   = sd(prop, na.rm = TRUE),
      sem  = sd / sqrt(n),
      .groups = "drop"
    )
  
  # 4. Statistics (LMM / KW)
  stats <- lapply(levels(df$cell), function(ct) {
    tmp <- df %>% filter(cell == ct)
    n_ds <- donor_n %>% filter(cell == ct) %>% pull(n_donor)
    
    if (all(n_ds >= min.lmm)) {
      fit <- try(lmer(prop ~ stage + (1 | donor), data = tmp), silent = TRUE)
      p <- if (!inherits(fit, "try-error")) anova(fit)$`Pr(>F)`[1] else NA
      data.frame(cell = ct, method = "LMM", p.value = p)
    } else if (length(unique(tmp$stage)) > 1) {
      kw <- kruskal.test(prop ~ stage, data = tmp)
      data.frame(cell = ct, method = "KW", p.value = kw$p.value)
    } else {
      data.frame(cell = ct, method = NA, p.value = NA)
    }
  }) %>% bind_rows()
  
  # 5. Plot
  fill.scale <- if (!is.null(stage.cols)) {
    scale_fill_manual(values = stage.cols)
  } else {
    scale_fill_brewer(palette = palette)
  }
  
  if (plot.type == "bar") {
    p <- ggplot(df, aes(stage, prop, fill = stage)) +
      stat_summary(fun = mean, geom = "bar",
                   color = "black", linewidth = 0.25, width = 0.7) +
      stat_summary(fun.data = mean_se, geom = "errorbar",
                   width = 0.25, linewidth = 0.3) +
      geom_jitter(aes(color = donor),
                  color = "black",    
                  stroke = 0.3,         
                  shape = 21, size = pt.size,
                  width = 0.1, alpha = 0.8)
  } else {
    p <- ggplot(df, aes(stage, prop, fill = stage)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = donor),
                  color = "black",
                  shape = 21, size = pt.size,
                  width = 0.1, alpha = 0.8)
  }
  
  if (!is.null(donor.cols)) {
    p <- p + scale_color_manual(values = donor.cols)
  }
  
  p <- p +
    fill.scale +
    facet_wrap(as.formula(paste("~", facet.by)),scales = "free_y",ncol = ncol) +
    plot_theme(theme.style = theme.style, leg.pos = leg.pos, ...) +
    labs(title = title, x = x.lab, y = y.lab) +
    theme(legend.position = leg.pos, strip.text = element_text(size = facet.title.size))
  
  # 6. P-value annotation
  if (show.p) {
    annot <- df %>%
      group_by(cell) %>%
      summarise(y = max(prop, na.rm = TRUE), .groups = "drop") %>%
      left_join(stats, by = "cell") %>%
      mutate(label = ifelse(is.na(p.value), "",
                            paste0("p=", signif(p.value, 2))))
    
    annot$cell <- factor(annot$cell, levels = levels(df$cell))
    
    p <- p + geom_text(
      data = annot,
      aes(x = 2, y = y * 1.05, label = label),
      inherit.aes = FALSE,
      size = label.size
    )
  }
  
  # 7. Source data (Nature-ready)
  source_data <- df %>%
    left_join(donor_n, by = c("cell", "stage")) %>%
    left_join(summary_stats, by = c("cell", "stage")) %>%
    left_join(stats, by = "cell") %>%
    arrange(cell, stage, donor)
  
  return(list(
    plot = p,
    source_data = source_data
  ))
}




#' Summarize and Plot Gene Set Expression
#'
#' @description
#' Computes the total expression of a given gene set across cells, optionally summarized 
#' in pseudobulk per group. Produces a bar plot showing summed expression per group/fill variable.
#'
#' @param object A Seurat object containing gene expression data.
#' @param gene.list A character vector of gene names to include.
#' @param group.var Metadata column to group cells by (x-axis). Default: "ann2".
#' @param fill.var Metadata column to use for fill colors (stacked or dodged bars). Default: "PCW".
#' @param sample.var Metadata column for pseudobulk (optional). If provided, aggregates by sample.
#' @param pseudobulk.mode Character. If "sum", sums expression per sample; if "cpm", normalizes counts per million. Default: NULL (no pseudobulk).
#' @param group.levels Optional factor levels for group.var.
#' @param fill.levels Optional factor levels for fill.var.
#' @param fill.colors Optional named vector of colors for fill levels.
#' @param expr.threshold Minimum expression to retain a gene. Default: 0.
#' @param theme ggplot2 theme. Default: "minimal".
#' @param leg.pos Legend position. Default: "right".
#' @param leg.dir Legend direction. Default: "vertical".
#' @param leg.ttl Legend title. Default: "".
#' @param font.size Base font size. Default: 10.
#' @param x.angle X-axis text angle. Default: 45.
#' @param x.lab X-axis label. Default: NULL.
#' @param y.lab Y-axis label. Default: "Total transcripts".
#' @param plot.ttl Plot title. Default: NULL.
#' @param flip Logical. If TRUE, flips coordinates. Default: FALSE.
#' @param return.data Logical. If TRUE, returns a list with plot and summarized data. Default: FALSE.
#' @param ... Additional arguments passed to `plot_theme`.
#'
#' @return A ggplot object (or list with plot and summarized data if `return.data = TRUE`).
#' @export
#'
#' @examples
#' genesum(
#'   object = seurat_obj,
#'   gene.list = c("OR1", "OR2"),
#'   group.var = "ann2",
#'   fill.var = "stage",
#'   sample.var = "sample",
#'   pseudobulk.mode = "cpm",
#'   plot.ttl = "OR pseudobulk across development"
#' )
#' 
genesum <- function(
    object,
    gene.list,
    group.var = "ann2",
    fill.var = "stage",
    sample.var = NULL,
    mode = c("cell", "pseudobulk"),
    pb.norm = c("none", "cpm"),
    group.levels = NULL,
    fill.levels = NULL,
    fill.colors = NULL,
    expr.threshold = 0,
    theme = "classic",
    leg.pos = "right",
    leg.dir = "vertical",
    leg.ttl = "",
    font.size = 10,
    x.angle = 45,
    x.lab = NULL,
    y.lab = "Total transcripts",
    plot.ttl = NULL,
    flip = FALSE,
    return.data = FALSE,
    ...
) {
  require(ggplot2)
  require(dplyr)
  require(reshape2)
  
  mode <- match.arg(mode)
  pb.norm <- match.arg(pb.norm)
  
  ## 1. Filter genes (FAST, sparse-safe)
  counts <- Seurat::GetAssayData(object, assay = "RNA", slot = "counts")
  
  valid.genes <- intersect(gene.list, rownames(counts))
  if (!length(valid.genes)) stop("No valid genes found")
  
  counts <- counts[valid.genes, , drop = FALSE]
  
  # Remove low-expression genes
  keep.genes <- Matrix::rowSums(counts) > expr.threshold
  counts <- counts[keep.genes, , drop = FALSE]
  
  ## 2. Total expression per cell (FAST)
  cell_totals <- Matrix::colSums(counts)
  
  meta <- object@meta.data
  meta$total_expr <- cell_totals[colnames(object)]
  
  ## 3. Aggregation
  if (mode == "cell") {
    
    summary.df <- meta %>%
      group_by(.data[[group.var]], .data[[fill.var]]) %>%
      summarise(
        total = sum(total_expr),
        n     = n(),
        sd    = sd(total_expr),
        sem   = sd / sqrt(n),
        .groups = "drop"
      )
  }
  
  if (mode == "pseudobulk") {
    
    if (is.null(sample.var))
      stop("sample.var must be provided for pseudobulk mode")
    
    pb <- meta %>%
      group_by(.data[[sample.var]],
               .data[[group.var]],
               .data[[fill.var]]) %>%
      summarise(total_expr = sum(total_expr), .groups = "drop")
    
    if (pb.norm == "cpm") {
      libsize <- pb %>%
        group_by(.data[[sample.var]]) %>%
        summarise(lib = sum(total_expr), .groups = "drop")
      
      pb <- pb %>%
        left_join(libsize, by = sample.var) %>%
        mutate(total_expr = total_expr / lib * 1e6)
    }
    
    summary.df <- pb %>%
      group_by(.data[[group.var]], .data[[fill.var]]) %>%
      summarise(
        total = sum(total_expr),
        n     = n(),
        sd    = sd(total_expr),
        sem   = sd / sqrt(n),
        .groups = "drop"
      )
  }
  
  ## 4. Factors
  if (!is.null(group.levels))
    summary.df[[group.var]] <- factor(summary.df[[group.var]], levels = group.levels)
  
  if (!is.null(fill.levels))
    summary.df[[fill.var]] <- factor(summary.df[[fill.var]], levels = fill.levels)
  
  ## 5. Plot
  pd <- position_dodge(width = 0.7)  # increase width to create space between bars within fill.var
  bar_width <- 0.5                    # set bar width
  
  min_frac <- 0.05 
  
  p <- ggplot(summary.df, aes(x = .data[[group.var]], y = total, fill = .data[[fill.var]])) +
    geom_bar(stat = "identity", position = pd, color = "black", linewidth = 0.25, width = bar_width) +
    geom_errorbar(aes(ymin = total - pmax(sem, min_frac * total),
                      ymax = total + pmax(sem, min_frac * total)),
                  position = pd, width = 0.25, linewidth = 0.25) +
    guides(fill = guide_legend(title = leg.ttl)) +
    labs(x = x.lab, y = y.lab, title = plot.ttl) +
    plot_theme(theme = theme, leg.pos = leg.pos, leg.dir = leg.dir,
               font.size = font.size, x.angle = x.angle, leg.size = font.size, ...)
  
  if (!is.null(fill.colors))
    p <- p + scale_fill_manual(values = fill.colors)
  
  if (flip)
    p <- p + coord_flip()
  
  if (return.data)
    return(list(plot = p, data = summary.df))
  
  return(p)
}





#'#' Highlight selected clusters on UMAP/tSNE plots
#'
#' @description
#' Plots a dimensionality reduction (UMAP/tSNE) and highlights specified clusters.
#' Can highlight clusters individually (mode = "single") or
#' simultaneously (mode = "multi").
#'
#' @param object Seurat object
#' @param cluster.names Character vector of clusters to highlight.
#'   If NULL, all identities are used.
#' @param group.by Metadata column used for cluster identity.
#' @param reduction Dimensional reduction name (default: "umap")
#' @param mode "single" = one plot per cluster, "multi" = combined highlight
#' @param ncol Number of columns for "single" mode layout
#' @param custom.colors Named vector of colors corresponding to highlighted clusters.
#' @param pt.size Baseline point size for non-highlight cells
#' @param highlight.size Point size for highlighted cells
#' @param fontsize Base text size
#' @param show.cell.counts Logical, adds cell numbers to legend
#' @param background.color Color for all non-highlighted cells
#' @param label Whether to label clusters (multi mode only)
#' @param label.size Cluster label font size
#' @param plot.title Custom plot title (NULL removes the default title)
#' @param no_axes Remove axes (recommended for UMAP/tSNE)
#' @param legend.position Legend position ("none" by default)
#' @param item.size Size of legend point markers
#' @param legend.ncol Number of legend columns
#' @param legend.title Title for legend (NULL hides it)
#' @param legend.text.size Legend text size
#' @param ... Additional arguments to DimPlot
#'
#' @return A ggplot2 or patchwork object
#' @export
#'
cellmark <- function(
    object,
    cluster.names = NULL,
    group.by = "seurat_clusters",
    reduction = "umap",
    mode = c("single", "multi"),
    ncol = 3,
    custom.colors = NULL,
    pt.size = 0.1,
    highlight.size = 0.5,
    fontsize = 10,
    show.cell.counts = FALSE,
    background.color = "lightgray",
    label = FALSE,
    label.size = 4,
    plot.title = NULL,
    no.axes = TRUE,
    legend.position = "none",
    item.size = 3,
    legend.ncol = 1,
    legend.title = NULL,
    legend.text.size = 10,
    ...
) {
  require(Seurat)
  require(ggplot2)
  require(patchwork)
  require(dplyr)
  require(ggrepel)
  
  mode <- match.arg(mode)
  object <- SetIdent(object, value = group.by)
  clusters <- levels(Idents(object))
  
  if (is.null(cluster.names)) cluster.names <- clusters
  valid.clusters <- cluster.names[cluster.names %in% clusters]
  if (length(valid.clusters) == 0)
    stop("No valid clusters found in metadata.")
  
  if (is.null(custom.colors)) {
    pal <- custom_palette(n = length(valid.clusters))
    names(pal) <- valid.clusters
    custom.colors <- pal
  }
  
  
  # MODE = SINGLE
  if (mode == "single") {
    plots <- lapply(valid.clusters, function(cl) {
      cells <- WhichCells(object, idents = cl)
      cnt <- length(cells)
      title <- if (show.cell.counts) paste0(cl, " (n=", cnt, ")") else cl
      
      p <- suppressWarnings({
        suppressMessages({DimPlot(
          object,
          reduction = reduction,
          group.by = group.by,
          cells.highlight = cells,
          cols = background.color,
          cols.highlight = custom.colors[cl],
          pt.size = pt.size,
          sizes.highlight = highlight.size,
          order = TRUE,
          ...
        ) })})+
        ggtitle(title) +
        theme(
          plot.title = element_text(hjust = 0.5, size = fontsize, face = "bold"),
          text = element_text(size = fontsize)
        ) &
        NoLegend()
      
      if (no_axes) p <- p & NoAxes()
      return(p)
    })
    
    return(wrap_plots(plots, ncol = ncol))
  }
  
  # MODE = MULTI
  cells.to.highlight <- CellsByIdentities(object, idents = valid.clusters)
  
  plt <- suppressWarnings({
    suppressMessages({DimPlot(
      object,
      reduction = reduction,
      group.by = group.by,
      cells.highlight = cells.to.highlight,
      cols.highlight = custom.colors[valid.clusters],
      cols = background.color,
      pt.size = pt.size,
      sizes.highlight = highlight.size,
      order = TRUE,
      ...
    ) }) })
  
  # Ensure consistent color scale ONLY once
  plt <- plt +
    scale_color_manual(
      values = custom.colors[valid.clusters],
      breaks = valid.clusters,
      na.value = background.color
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          size = item.size,
          shape = 21,
          fill = unname(custom.colors[valid.clusters]),
          stroke = 0.3,
          color = "black"
        ),
        ncol = legend.ncol,
        title = legend.title
      )
    ) #+
  # theme(
  #   legend.position = legend.position,
  #   legend.text = element_text(size = legend.text.size),
  #   text = element_text(size = fontsize)
  # )
  
  plt <- plt +
    theme(
      legend.position = legend.position,
      legend.text = element_text(size = legend.text.size),
      legend.title = element_text(size = legend.text.size, face = "bold"),
      
      legend.spacing.y = unit(1, "pt"),
      legend.key.height = unit(4, "pt"),
      legend.key.width  = unit(4, "pt"),
      legend.margin = margin(1, 1, 1, 1),
      
      text = element_text(size = fontsize)
    )
  
  # Plot title rule 
  if (is.null(plot.title)) {
    plt <- plt + theme(plot.title = element_blank())
  } else {
    plt <- plt + ggtitle(plot.title) +
      theme(plot.title = element_text(
        hjust = 0.5, size = fontsize + 2, face = "bold"
      ))
  }
  
  # Labels on highlighted clusters
  if (label) {
    emb <- Embeddings(object[[reduction]])
    df <- data.frame(emb, cluster = Idents(object)) %>%
      filter(cluster %in% valid.clusters) %>%
      group_by(cluster) %>%
      summarise(x = median(umap_1), y = median(umap_2), .groups = "drop")
    
    plt <- plt +
      geom_text_repel(
        data = df,
        aes(x = x, y = y, label = cluster),
        size = label.size,
        fontface = "bold"
      )
  }
  
  # Cell count labels in legend
  if (show.cell.counts) {
    cell.nb <- table(Idents(object))
    lbl <- paste0(valid.clusters, " (",cell.nb[valid.clusters], ")")
    names(lbl) <- valid.clusters
    
    plt <- plt +
      scale_color_manual(
        values = custom.colors[valid.clusters],
        breaks = valid.clusters,
        labels = lbl,
        na.value = background.color
      )
  }
  
  if (no.axes) plt <- plt & NoAxes()
  return(plt)
}



 #' Cell Feature Plot for Seurat Objects
#'
#' A flexible wrapper around Seurat's `FeaturePlot` to visualize gene expression or metadata 
#' features in a Seurat object. Supports custom color palettes, viridis, and hotspot/rainbow palettes.
#'
#' @param object A `Seurat` object.
#' @param features Character vector of features (genes or metadata columns) to plot.
#' @param cols Optional character vector of colors for plotting.
#' @param theme.cols Character. Predefined theme color palette (default: `"Reds"`). Options include `"Reds"`, `"Blues"`, etc., or custom list palettes `"hotspot"` and `"rainbow"`.
#' @param viridis Logical. If TRUE, use viridis palette instead of RColorBrewer or custom palettes.
#' @param viridis.opt Character. Viridis palette option (default: `"D"`).
#' @param rev.cols Logical. Reverse the color palette (default: `FALSE`).
#' @param na.col Color for NA or below-cutoff expression values (default: `"lightgray"`).
#' @param order Logical. If TRUE, plot high-expression cells on top (default: `FALSE`).
#' @param pt.size Numeric. Point size. If NULL, automatically calculated based on number of cells.
#' @param font.size Numeric. Base font size for plot titles and axis labels (default: 10).
#' @param reduction Character. Dimensional reduction to use (default: first available in Seurat object).
#' @param na.cutoff Numeric. Minimum expression value for coloring; below this will be NA if palette requires (default: 1e-9).
#' @param raster Logical. If TRUE, rasterize points for faster plotting of large datasets.
#' @param raster.dpi Numeric vector of length 2. DPI for rasterization (default: c(512,512)).
#' @param split.by Character. Metadata column to split the plot.
#' @param ncol Numeric. Number of columns when combining multiple plots.
#' @param layer Character. Seurat assay slot to fetch data from (default: `"data"`).
#' @param label Logical. Whether to label clusters (default: FALSE).
#' @param axes Logical. Whether to show axes (default: TRUE).
#' @param combine Logical. Whether to return a single combined plot (default: TRUE).
#' @param blend Logical. Whether to blend exactly two features (default: FALSE).
#' @param merge.leg Logical. Whether to merge multiple legends into one (default: FALSE).
#' @param theme.style Character. ggplot2 theme to apply (default: `"classic"`).
#' @param ... Additional arguments passed to `Seurat::FeaturePlot`.
#'
#' @return A `ggplot` object (or `patchwork` object if multiple features).
#' @export
#'
#' @examples
#' # Single feature with default Reds palette
#' cellfeat(seurat_obj, features = "POMC")
#'
#' # Multiple features with viridis palette and merged legend
#' cellfeat(seurat_obj, features = c("POMC", "NPY"), viridis = TRUE, merge.leg = TRUE)
#'
#' # Blend two features
#' cellfeat(seurat_obj, features = c("POMC", "NPY"), blend = TRUE)
#' 
cellfeat <- function(
    object,
    features,
    cols = NULL,
    theme.cols = "Reds",
    viridis = FALSE,
    viridis.opt = "D",
    rev.cols = FALSE,
    na.col = "lightgray",
    order = FALSE,
    pt.size = NULL,
    font.size = 10,
    reduction = NULL,
    na.cutoff = 1e-9,
    raster = NULL,
    raster.dpi = c(512, 512),
    split.by = NULL,
    ncol = NULL,
    layer = "data",
    label = FALSE,
    axes = TRUE,
    combine = TRUE,
    blend = FALSE,
    merge.leg = FALSE,
    theme.style = "classic",
    ...
) {
  
  suppressWarnings(suppressMessages({
    
    # Checks
    if (!inherits(object, "Seurat"))
      stop("object must be a Seurat object")
    
    features <- intersect(
      unique(as.character(features)),
      c(rownames(object), colnames(object@meta.data))
    )
    
    if (length(features) == 0)
      stop("No valid features found")
    
    reduction <- reduction %||% SeuratObject::DefaultDimReduc(object)
    
    # Colors
    # if (is.null(cols)) {
    #   cols <- RColorBrewer::brewer.pal(9, theme.cols)
    # }
    # if (rev.cols) cols <- rev(cols)
    
    # Colors
    if (!is.null(cols)) {
      if (length(cols) < 2)
        stop("colors_use must contain at least 2 colors")
      cols <- cols
      if (rev.cols) cols <- rev(cols)
      
    } else if (isTRUE(viridis)) {
      
      if (!requireNamespace("viridis", quietly = TRUE))
        stop("Package 'viridis' is required when viridis = TRUE")
      
      cols <- viridis::viridis(9, option = viridis.opt)
      if (rev.cols) cols <- rev(cols)
      
    } else {
      
      pal <- list(
        hotspot = c("navy", "lightblue", "yellow", "orange", "red"),
        rainbow = c("blue", "green", "yellow", "orange", "red")
      )
      
      if (theme.cols %in% names(pal)) {
        cols <- pal[[theme.cols]]
      } else {
        cols <- RColorBrewer::brewer.pal(9, theme.cols)
      }
      
      if (rev.cols) cols <- rev(cols)
    }
    
    
    
    # Point size
    raster <- raster %||% (ncol(object) > 2e5)
    if (is.null(pt.size)) {
      pt.size <- if (raster) 1 else min(1583 / ncol(object), 1)
    }
    
    ## Global max for merged legend
    global_max <- NA
    if (merge.leg && is.null(split.by) && length(features) > 1 && !blend) {
      expr_data <- Seurat::FetchData(object, vars = features, layer = layer)
      global_max <- max(expr_data, na.rm = TRUE)
    }
    
    # Plot
    plt <- Seurat::FeaturePlot(
      object = object,
      features = features,
      reduction = reduction,
      order = order,
      pt.size = pt.size,
      raster = raster,
      raster.dpi = raster.dpi,
      split.by = split.by,
      combine = combine,
      blend = blend,
      label = label,
      ncol = ncol,
      ...
    )
    
    # Apply color scale (always)
    # if (!blend) {
    #   if (merge.leg && is.null(split.by) && length(features) > 1) {
    #     # Shared global scale
    #     plt <-  plt & scale_color_gradientn(
    #       colors = cols,
    #       limits = c(na.cutoff, global_max),
    #       na.value = na.col
    #     )
    #   } else {
    #     # Per-feature scale
    #     plt <- plt & scale_color_gradientn(
    #       colors = cols,
    #       limits = c(na.cutoff, NA),
    #       na.value = na.col,
    #       name = NULL
    #     )
    #   }
    # }
    
    # Determine palette behavior
    like <- isTRUE(viridis) || theme.cols %in% c("hotspot", "rainbow")
    
    # Brewer palettes should gray-out low expression
    if (!like && !blend) {
      
      plt <- plt & scale_color_gradientn(
        colors = cols,
        limits = c(na.cutoff, NA),
        oob = scales::censor,     # force < cutoff → NA
        na.value = na.col,        # show gray
        name = NULL
      )} else if (!blend) {
      # No gray background
      plt <- plt & scale_color_gradientn(
        colors = cols,
        limits = c(na.cutoff, NA),
        oob = scales::squish,     # clip instead of NA
        na.value = NA
      )
    }
    
    # Theme
    plt <- plt & plot_theme(theme.style = theme.style, font.size = font.size, ...) &
      theme(plot.title = element_text(
        hjust = 0.5, size = font.size + 2, face = "bold.italic"
      ))
    
    # ONLY add colorbar guide for non-blend plots
    if (!blend) {
      plt <- plt & guides(
        color = guide_colorbar(
          frame.colour = "black",
          ticks.colour = "black"
        )
      )
    }
    
    # Manage axes
    if (!axes) {
      plt <- plt & Seurat::NoAxes()
    }
    
    ## Merge legend 
    if (merge.leg && is.null(split.by) && length(features) > 1) {
      plt <- Seurat::CombinePlots(plots = plt, legend = "right", ncol = ncol)
    }
    
    return(plt)
  }))
}



#' Generate Expression Heatmaps for Stage or Celltype–Sample Modes
#'
#' This function generates heatmaps from a Seurat object using either:
#' 1) **Stage mode**: aggregated expression across developmental or experimental stages; or  
#' 2) **Celltype–Sample mode**: paired annotation of cell types and samples.  
#'
#' It supports z-score normalization, gene filtering, palette selection, optional
#' flipping of heatmaps, row/column label position control, and integrated saving
#' through a custom `save.fig()` function.
#'
#' @param object A Seurat object.
#' @param features Character vector of gene names to include.
#' @param mode Plot mode: `"stage"` (default) or `"celltype-sample"`.
#' @param assay Assay to use (default `"RNA"`).
#' @param slot Data slot to extract (default `"counts"`).
#' @param zlim Numeric vector defining z-score color scaling boundaries
#'   (default `c(-2, 0, 2)`).
#' @param fontsize Base font size for labels.
#' @param palette Color palette. Either `"viridis"` or manual colors.
#' @param manual.colors Colors used when `palette != "viridis"`.
#' @param flip.heatmap Logical. If TRUE, transpose the heatmap and switch label sides.
#' @param group.by Metadata column used for averaging in `"stage"` mode.
#' @param column.title Optional column title for the heatmap.
#' @param celltype.col Metadata column defining cell types (for `"celltype-sample"` mode).
#' @param sample.col Metadata column defining samples (for `"celltype-sample"` mode).
#' @param celltype.palette Named color vector for cell types.
#' @param sample.palette Named color vector for samples.
#' @param show.annotation.legend Logical; show annotation legends (default TRUE).
#' @param row.names.side Side for row names: `"left"` or `"right"` (default `"left"`).
#' @param column.names.side Side for column names: `"top"` or `"bottom"` (default `"top"`).
#' @param row.names.rot Rotation angle of row labels.
#' @param column.names.rot Rotation angle of column labels.
#' @param output.pdf File path to save heatmap as PDF. If NULL, no file is saved.
#' @param pdf.width PDF width (inches).
#' @param pdf.height PDF height (inches).
#' @param heatmap.legend.side Side for heatmap legend.
#' @param annotation.legend.side Side for annotation legend.
#'
#' @details  
#' The function uses Seurat's `AverageExpression()` to compute aggregated values,
#' applies per-gene z-scoring, and uses `ComplexHeatmap` for rendering.  
#'
#' If `flip.heatmap = TRUE`, the matrix is transposed and label sides are swapped
#' automatically to maintain readability.  
#'
#' If `output.pdf` is provided, the function saves the heatmap using the user’s
#' custom `save.fig()` helper, ensuring directory creation and figure consistency.
#'
#' @return Invisibly returns the ComplexHeatmap object used for plotting.
#'
#' @examples
#' \dontrun{
#' expression_heatmap(
#'   object,
#'   features = c("GATA3", "SOX2"),
#'   mode = "stage",
#'   group.by = "PCW.group",
#'   output.pdf = "plots/stage.heatmap.pdf"
#' )
#' }
#'
#' @import Seurat
#' @import dplyr
#' @import ComplexHeatmap
#' @import viridis
#' @import circlize
#' @import stringr
#' @import grid
#'
#' @export
#' 
expression_heatmap <- function(
    object,
    features,
    mode = c("stage", "celltype-sample"),
    assay = "RNA",
    slot = "counts",
    zlim = c(-2, 0, 2),
    fontsize = 12,
    palette = "viridis",
    manual.colors = c("purple","black","yellow"),
    flip.heatmap = FALSE,
    group.by = "PCW_group",
    column.title = NULL,
    celltype.col = "ann_level_2",
    sample.col = "sample",
    celltype.palette = NULL,
    sample.palette = NULL,
    show.annotation.legend = TRUE,
    row.names.side = "left",
    column.names.side = "top",
    row.names.rot = 0,
    column.names.rot = 90,
    output.pdf = NULL,
    pdf.width = 8,
    pdf.height = 6,
    heatmap.legend.side = "right",
    annotation.legend.side = "right"
){
  
  require(Seurat)
  require(ComplexHeatmap)
  require(circlize)
  require(viridis)
  require(stringr)
  require(grid)
  
  mode <- match.arg(mode)
  
  ## helpers 
  italic_labels <- function(x) {
    as.expression(lapply(x, function(g) bquote(italic(.(g)))))
  }
  
  make_labels <- function(mat, flip) {
    if (!flip) {
      list(row = italic_labels(rownames(mat)), col = NULL)
    } else {
      list(row = NULL, col = italic_labels(colnames(mat)))
    }
  }
  
  make_col_fun <- function() {
    if (palette == "viridis") {
      colorRamp2(seq(zlim[1], zlim[3], length.out = 256), viridis(256))
    } else {
      colorRamp2(zlim, manual.colors)
    }
  }
  
  ## gene filtering 
  genes <- intersect(features, rownames(object))
  if (length(genes) == 0) stop("None of the genes exist in Seurat object.")
  
  sub <- subset(object, features = genes)
  DefaultAssay(sub) <- assay
  
  ## compute average expression 
  if (mode == "stage") {
    avg <- AverageExpression(
      sub, group.by = group.by, assays = assay, slot = slot
    )[[assay]]
  } else {
    avg <- AverageExpression(
      sub, group.by = c(celltype.col, sample.col),
      assays = assay, slot = slot
    )[[assay]]
  }
  
  ## z-score & ordering
  avg <- avg[!rowSums(is.na(avg)), , drop = FALSE]
  avg <- t(scale(t(avg)))
  avg[is.na(avg)] <- 0
  avg <- avg[order(apply(avg, 1, which.max)), , drop = FALSE]
  
  ## annotations (celltype × sample)
  top.anno <- NULL
  if (mode == "celltype-sample") {
    
    parts <- str_split_fixed(colnames(avg), "_", 2)
    anno.df <- data.frame(
      CellTypes = factor(parts[,1], levels = levels(object[[celltype.col]][,1])),
      Samples   = factor(parts[,2], levels = levels(object[[sample.col]][,1]))
    )
    rownames(anno.df) <- colnames(avg)
    
    if (is.null(celltype.palette)) {
      celltype.palette <- setNames(
        viridis(length(levels(anno.df$CellTypes))),
        levels(anno.df$CellTypes)
      )
    }
    if (is.null(sample.palette)) {
      sample.palette <- setNames(
        viridis(length(levels(anno.df$Samples))),
        levels(anno.df$Samples)
      )
    }
    
    top.anno <- HeatmapAnnotation(
      df = anno.df,
      col = list(CellTypes = celltype.palette, Samples = sample.palette),
      show_legend = show.annotation.legend,
      annotation_name_gp = gpar(fontsize = fontsize, fontface = "bold")
    )
  }
  
  ## flip heatmap 
  if (flip.heatmap) {
    avg <- t(avg)
    row.names.side    <- ifelse(row.names.side == "left", "right", "left")
    column.names.side <- ifelse(column.names.side == "top", "bottom", "top")
  }
  
  labels <- make_labels(avg, flip.heatmap)
  
  ## build heatmap 
  ht_args <- list(
    matrix = avg,
    name = "Z-score",
    top_annotation = top.anno,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = make_col_fun(),
    border = TRUE,
    row_names_side = row.names.side,
    column_names_side = column.names.side,
    row_names_rot = row.names.rot,
    column_names_rot = column.names.rot,
    row_names_gp = gpar(fontsize = fontsize),
    column_names_gp = gpar(fontsize = fontsize),
    column_title = column.title,
    heatmap_legend_param = list(
      title = "Z-score",
      title_gp = gpar(fontface = "bold"),
      border = "black"
    )
  )
  
  if (!flip.heatmap) {
    ht_args$row_labels <- italic_labels(rownames(avg))
  } else {
    ht_args$column_labels <- italic_labels(colnames(avg))
  }
  
  ht <- do.call(Heatmap, ht_args)
  
  ## draw / save 
  if (!is.null(output.pdf)) {
    
    out.dir <- dirname(output.pdf)
    if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
    
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = heatmap.legend.side,
      annotation_legend_side = annotation.legend.side,
      merge_legend = TRUE
    )
    
    save_heatmap(
      ht = ComplexHeatmap::draw(
        ht,
        heatmap_legend_side = heatmap.legend.side,
        annotation_legend_side = annotation.legend.side,
        merge_legend = TRUE
      ),
      filename = output.pdf,
      formats = c("pdf", "png", "tiff"),
      width = pdf.width,
      height = pdf.height,
      dpi = 600
    )
    
    cat("Saved heatmap to:", output.pdf, "\n")
    
  } else {
    
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = heatmap.legend.side,
      annotation_legend_side = annotation.legend.side,
      merge_legend = TRUE
    )
  }
}


#'#' Plot a Transcription Factor Regulatory Network
#'
#' Visualizes a transcription factor (TF) regulatory network using selected source TFs and their top N targets based on absolute regulatory strength (mor score).
#' Supports visual differentiation of positive and negative regulation, with shape and color coding for source vs. target nodes.
#'
#' @param net A data frame representing the TF regulatory network with columns: \code{source}, \code{target}, and \code{mor} (mode of regulation).
#' @param selected_sources A character vector of TFs to be used as source nodes.
#' @param n_targets Integer specifying the number of top targets to display per source (default: 5).
#' @param label_size Numeric; font size of node labels (default: 2).
#' @param node_size Numeric; size of node points (default: 15).
#' @param repel Logical; whether to repel node labels to avoid overlap using \code{geom_node_text} (default: FALSE).
#' @param layout_type Layout type for the network graph; passed to \code{ggraph::ggraph()} (default: "fr").
#' @param positive_color Color for positively regulated targets (default: "darkgreen").
#' @param negative_color Color for negatively regulated targets (default: "#e59866").
#' @param plot_title Title of the plot (default: "TF Regulatory Network").
#'
#' @return A ggplot2 object representing the regulatory network.
#'
#' @import ggraph
#' @import igraph
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' plot_tf_network(net = tf_network_df,
#'                 selected_sources = c("SOX2", "FOXA1"),
#'                 n_targets = 10,
#'                 repel = TRUE)
#' }
plot_tf_network <- function(net, 
                            selected_sources,
                            n_targets = 5, 
                            label_size = 2, 
                            node_size = 15,
                            linewidth = 1.2,
                            repel = FALSE,
                            layout_type = "graphopt",   # deterministic
                            max_edges = 500,
                            positive_color = "darkred", 
                            negative_color = "navy",
                            plot_title = "TF Regulatory Network") {
  
  require(ggraph)
  require(igraph)
  require(dplyr)
  
  # Filter network
  filtered_net <- net %>%
    filter(source %in% selected_sources) %>%
    group_by(source) %>%
    top_n(n = n_targets, wt = abs(mor)) %>%
    ungroup()
  
  if (nrow(filtered_net) > max_edges) {
    filtered_net <- filtered_net %>%
      arrange(desc(abs(mor))) %>%
      slice(1:max_edges)
    warning(paste("Edge count exceeds", max_edges, "— truncating to top", max_edges, "edges."))
  }
  
  # Build graph
  g <- graph_from_data_frame(filtered_net, directed = TRUE)
  filtered_net$edge_color <- ifelse(filtered_net$mor > 0, "Positive", "Negative")
  
  node_df <- data.frame(name = V(g)$name) %>%
    mutate(node_type = ifelse(name %in% selected_sources, "Source", "Target"),
           shape = ifelse(node_type == "Source", 22, 21),
           edge_color = ifelse(node_type == "Source", "Source",
                               ifelse(name %in% filtered_net$target[filtered_net$mor > 0], "Positive", "Negative")),
           label_color = ifelse(node_type == "Source", "black", "white"))
  
  node_color_map <- c("Positive" = positive_color,
                      "Negative" = negative_color,
                      "Source" = "white")
  
  # FIXED LAYOUT → stable, reproducible!
  set.seed(42)   # <-- makes layout identical every run
  
  
  # Plot
  p <- ggraph(g, layout = layout_type) +
    geom_edge_link0(aes(edge_alpha = abs(mor), edge_colour = filtered_net$edge_color),
                    show.legend = TRUE,  linewidth = linewidth) +
    
    geom_node_point(aes(fill = factor(name, levels = node_df$name, labels = node_df$edge_color)),
                    size = node_size, shape = node_df$shape, show.legend = FALSE) +
    
    geom_node_text(aes(label = name, color = factor(name, levels = node_df$name, labels = node_df$label_color)),
                   size = label_size, fontface = "bold.italic", repel = repel) +
    
    
    # Legend 
    scale_edge_color_manual(
      name = "Regulation",
      values = c("Positive" = positive_color, "Negative" = negative_color),
      guide = guide_legend(
        override.aes = list(size = 6),
        title.position = "top",
        title.theme = element_text(size = 18, face = "bold"),
        label.theme = element_text(size = 16)
      )
    ) +
    
    scale_fill_manual(values = node_color_map, guide = "none") +
    scale_color_manual(values = c("black", "white"), guide = "none") +
    guides(edge_alpha = "none") +
    
    theme_void() +
    ggtitle(plot_title) +
    
    theme(
      plot.margin = unit(c(1,1,1,1), "cm"),
      legend.text  = element_text(size = 16),
      legend.title = element_text(size = 18, face = "bold"),
      legend.key.size = unit(1.4, "cm"),
      legend.key.height = unit(1.4, "cm"),
      legend.key.width = unit(1.4, "cm"),
      legend.position = "right"
    ) +
    
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
    scale_y_discrete(expand = expansion(mult = c(0.1, 0.1)))
  
  return(p)
}


plot_tf_activity <- function(
    seurat_object,
    tf_list,
    assay_activity = "tfsulm",
    group.by = "seurat_clusters",
    low_color = "navy",
    high_color = "darkred",
    zscore = TRUE,
    zscore_mode = c("per_tf", "global"),
    fontsize = 10,
    legend.text.size = 10,
    theme = "classic",
    ncol = 1,
    leg.pos = "bottom",
    return_data = TRUE
) {
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(patchwork)
  
  zscore_mode <- match.arg(zscore_mode)
  
  DefaultAssay(seurat_object) <- assay_activity
  found_activity <- intersect(tf_list, rownames(GetAssayData(seurat_object, slot = "data")))
  
  if (length(found_activity) == 0)
    stop("None of the TFs were found in the assay.")
  
  # Compute TF activity data (ACCUMULATE CORRECTLY)
  all_tf_data <- list()
  
  for (tf in found_activity) {
    
    df <- FetchData(seurat_object, vars = c(tf, group.by)) %>%
      group_by(.data[[group.by]]) %>%
      summarise(
        mean_val = mean(.data[[tf]], na.rm = TRUE),
        n_cells  = n(),
        .groups = "drop"
      ) %>%
      mutate(TF = tf)
    
    df$zval <- if (zscore) scale(df$mean_val)[, 1] else df$mean_val
    
    all_tf_data[[tf]] <- df
  }
  
  tf_activity_df <- bind_rows(all_tf_data)
  
  # Global z-score limits (if requested)
  if (zscore && zscore_mode == "global") {
    global_abs_lim <- max(abs(tf_activity_df$zval), na.rm = TRUE)
  } else {
    global_abs_lim <- NULL
  }
  
  
  # Plotting
  plot_list <- list()
  total_tf <- length(tf_list)
  
  for (i in seq_along(tf_list)) {
    
    tf <- tf_list[i]
    df <- tf_activity_df %>% filter(TF == tf)
    
    if (nrow(df) == 0) {
      p <- ggplot() + theme_void() + ggtitle(paste0(tf, " (not found)"))
    } else {
      
      show_x <- i == total_tf
      
      p <- ggplot(df, aes(x = .data[[group.by]], y = zval, fill = zval)) +
        geom_col(color = "black", linewidth = 0.2, width = 0.55) +
        scale_fill_gradient2(
          low = low_color,
          mid = "white",
          high = high_color,
          midpoint = 0,
          limits = if (!is.null(global_abs_lim)) c(-global_abs_lim, global_abs_lim) else NULL,
          name = if (zscore) "z-score" else "Mean activity"
        ) +
        plot_theme(font.size = fontsize, theme = theme, axis.title.x.show = FALSE) +
        theme(
          axis.text.x  = if (show_x) element_text(angle = 45, hjust = 1) else element_blank(),
          axis.ticks.x = if (show_x) element_line() else element_blank(),
          axis.title.x = if (show_x) element_text() else element_blank(),
          plot.title   = element_text(face = "bold.italic", hjust = 0.5)
        ) +
        ggtitle(tf) +
        ylab(if (zscore) "z-score" else "Mean TF activity") +
        xlab(if (show_x) group.by else "")
    }
    
    plot_list[[tf]] <- p
  }
  
  final_plot <- wrap_plots(
    plot_list,
    ncol = ncol,
    guides = "collect",
    axis_titles = "collect"
  ) &
    theme(
      legend.position = leg.pos,
      legend.title    = element_text(size = legend.text.size),
      legend.text     = element_text(size = legend.text.size)
    )
  
  # Return
  if (return_data) {
    return(list(
      plot = final_plot,
      data = tf_activity_df
    ))
  } else {
    return(final_plot)
  }
}


#' Plot dominance score results
#'
#' @param df Data frame returned by `dominance_score()`$df
#' @param mode Plot mode: "bins", "high_cells", "top_genes"
#' @param ident.col Grouping column for "bins" and "high_cells" (default "ident")
#' @param stage.col Optional column to facet/stratify plots by, e.g., "stage" or "sample"
#' @param bin.col Column name for dominance bins (default "adjusted_bin")
#' @param gene.col Column name for genes (default "top_gene")
#' @param min.adjust Minimum adjusted dominance score (used in high_cells & top_genes)
#' @param min.top.expr Minimum top gene expression
#' @param max.n.genes Maximum number of expressed genes
#' @param top.n Number of top genes to show (mode = "top_genes")
#' @param colors Optional named vector of colors
#' @param font.size Base font size
#' @param theme ggplot theme name
#' @param leg.pos Legend position
#' @param x.angle X-axis text angle
#'
#' @return ggplot object
#' @export
plot_domscore <- function(
    df,
    mode = c("bins", "high_cells", "top_genes"),
    ident.col = "ident",
    stage.col = NULL,       
    bin.col = "adjusted_bin",
    gene.col = "top_gene",
    min.adjust = 1,
    min.top.expr = 1,
    max.n.genes = 3,
    top.n = 20,
    colors = NULL,
    font.size = 10,
    theme = "classic",
    plot.ttl = NULL,
    leg.pos = "top",
    leg.dir = "horizontal",
    x.angle = 45,
    ncol = 2,
    flip = TRUE,
    facet.bg = FALSE,
    ...
) {
  require(dplyr)
  require(ggplot2)
  mode <- match.arg(mode)
  
  if (is.null(colors)) {
    groups <- unique(df[[ident.col]])
    colors <- setNames(
      grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(groups)),
      groups
    )
  }
  
  # MODE: dominance bins
  if (mode == "bins") {
    plot.df <- df %>%
      dplyr::filter(!is.na(.data[[bin.col]])) %>%
      {
        if (!is.null(stage.col) && stage.col %in% colnames(.)) {
          group_by(., .data[[stage.col]], .data[[bin.col]], .data[[ident.col]]) 
        } else {
          group_by(., .data[[bin.col]], .data[[ident.col]])
        }
      } %>%
      summarise(n = n(), .groups = "drop") %>%
      {
        if (!is.null(stage.col) && stage.col %in% colnames(.)) {
          group_by(., .data[[stage.col]], .data[[bin.col]])
        } else group_by(., .data[[bin.col]])
      } %>%
      mutate(perc = 100 * n / sum(n)) %>%
      ungroup()
    
    p <- ggplot(plot.df,
                aes(x = .data[[bin.col]], y = perc, fill = .data[[ident.col]])) +
      geom_col(color = "black", linewidth = 0.2, width=0.5) +
      scale_fill_manual(values = colors) +
      labs(x = "Dominance score", y = "Cell type proportion (%)", fill = "", title = plot.ttl) +
      plot_theme(theme = theme, font.size = font.size, leg.pos = leg.pos,leg.dir = leg.dir, x.angle = 0)
  }
  
  # MODE: high-dominance cells
  else if (mode == "high_cells") {
    plot.df <- df %>%
      filter(
        adjust_score > min.adjust,
        top_expr >= min.top.expr,
        n_genes_expr <= max.n.genes
      ) %>%
      {
        if (!is.null(stage.col) && stage.col %in% colnames(.)) {
          group_by(., .data[[stage.col]], .data[[ident.col]])
        } else group_by(., .data[[ident.col]])
      } %>%
      summarise(n_cells = n(), .groups = "drop")
    
    x_col <- if (!is.null(stage.col) && stage.col %in% colnames(df)) stage.col else ident.col
    
    p <- ggplot(plot.df,
                aes(x = factor(.data[[x_col]]),
                    y = n_cells,
                    fill = .data[[ident.col]])) +
      geom_col(color = "black", linewidth = 0.2, width = 0.5) +
      scale_fill_manual(values = colors) +
      labs(x = "", y = paste0("High-dominance cells\n(adjusted > ", min.adjust, ")"), fill = "", title = plot.ttl) +
      plot_theme(theme = theme, font.size = font.size, leg.pos = leg.pos, leg.dir = leg.dir, x.angle = x.angle)
  }
  
  
  # MODE: top dominant genes
  if (mode == "top_genes") {
    high.df <- df %>%
      filter(
        adjust_score > min.adjust,
        top_expr >= min.top.expr,
        n_genes_expr <= max.n.genes,
        !is.na(.data[[gene.col]])
      )
    
    high.df[[gene.col]] <- as.character(high.df[[gene.col]])
    
    if (!is.null(stage.col) && stage.col %in% colnames(high.df)) {
      plot.df <- high.df %>%
        group_by(.data[[stage.col]], .data[[gene.col]]) %>%
        summarise(n = n(), .groups = "drop") %>%
        arrange(.data[[stage.col]], desc(n)) %>%
        group_by(.data[[stage.col]]) %>%
        slice_head(n = top.n) %>%
        ungroup()
      
      plot.df$n <- as.integer(plot.df$n)
      
      p <- ggplot(plot.df, aes(x = reorder(.data[[gene.col]], -n), y = n)) +  # <- descending
        geom_col(fill = "steelblue", color = "black", linewidth = 0.2, width = 0.5) +
        facet_wrap(as.formula(paste0("~", stage.col)), ncol = ncol) +
        scale_y_continuous(breaks = scales::pretty_breaks()) +
        labs(x = "Top dominant OR genes", y = "Number of high-dominance cells", title = plot.ttl) +
        plot_theme(font.size = font.size, theme.style = theme, leg.pos = "none", x.angle = x.angle, facet.bg = facet.bg,...)
      
      if (flip) p <- p + coord_flip()
    } else {
      plot.df <- high.df %>%
        dplyr::count(!!rlang::sym(gene.col), sort = TRUE, name = "n") %>%
        slice_head(n = top.n)
      
      plot.df$n <- as.integer(plot.df$n)
      
      p <- ggplot(plot.df, aes(x = reorder(!!rlang::sym(gene.col), -n), y = n)) +  # <- descending
        geom_col(fill = "steelblue", color = "black", linewidth = 0.2, width = 0.5) +
        scale_y_continuous(breaks = scales::pretty_breaks()) +
        labs(x = "", y = "Number of high-dominance cells", title = plot.ttl) +
        plot_theme(font.size = font.size, theme.style = theme, leg.pos = "none", x.angle = x.angle,facet.bg = facet.bg,...)
      
      if (flip) p <- p + coord_flip()
    }
  }
  if (!is.null(colors)) p <- p + scale_fill_manual(values = colors)
  return(p)
}


#'------------------------------------------------------------------------------
#'------------------------------------------------------------------------------
#'
#' Plot Heatmap Along Pseudotime
#'
#' Generates a heatmap of gene expression along pseudotime for selected lineages, with
#' flexible annotation and legend options. Supports automatic legend type detection
#' (discrete vs continuous), smooth scaling, z-score transformation, and flexible
#' legend title positioning.
#'
#' @param heatmap_matrices List of heatmap matrices produced by `identify_temporal_genes` or similar. Each element should have `expr_raw`, `ordered_cells`, `top_genes`.
#' @param seurat_obj Seurat object containing cell metadata.
#' @param lineage_name Character. Name of the lineage to plot. Defaults to first element of `heatmap_matrices`.
#' @param pseudotime_vec Numeric vector of pseudotime values. Defaults to sequential ordering of `ordered_cells`.
#' @param gene_list Character vector of genes to plot. Defaults to `hm$top_genes`.
#' @param annotation_col Column in `seurat_obj@meta.data` to use for annotation. Default `"ann_level_2"`.
#' @param celltype_colors Named vector of colors for annotation clusters. If `NULL`, colors are auto-generated.
#' @param zscore_limits Numeric vector of length 3, c(min, mid, max), used for color scaling. Default c(-2,0,2).
#' @param pseudotime_palette Character. Palette name for pseudotime continuous scale. Default `"YlGnBu"`.
#' @param row_cluster_colors Not currently used. Placeholder for future row annotations.
#' @param font_size Font size for heatmap labels and legends. Default 10.
#' @param output_pdf Path to save heatmap. If `NULL`, heatmap is drawn on screen.
#' @param fig_name Filename for saving heatmap. Default uses `plot_title`.
#' @param pdf_width Width of saved PDF/PNG. Default 6.
#' @param pdf_height Height of saved PDF/PNG. Default 10.
#' @param plot_title Title of heatmap. Default `NULL`.
#' @param show_annotation_legend Logical. Show legend for top annotations. Default `TRUE`.
#' @param smooth_df Degrees of freedom for smoothing spline. Default 3.
#' @param legend_dir `"vertical"` or `"horizontal"`. Orientation of legend. Default `"vertical"`.
#' @param ht_legend_side Side of heatmap legend (`"right"` or `"left"`). Default `"right"`.
#' @param ann_legend_side Side of top annotation legend (`"top"` or `"right"`). Default `"right"`.
#' @param legend_width_cm Width of legends in cm. Default 3.
#' @param cont_legend_cm Height of continuous legends in cm. Default 3.
#' @param merge_legends Logical. Merge annotation and heatmap legends. Default `FALSE`.
#' @param center_legend_titles Logical. Center all heatmap legend titles. Default `TRUE`.
#' @param center_cells_title Logical. Center "Cells" annotation legend title. Default `TRUE`.
#' @param center_pt_title Logical. Center "Pseudotime" annotation legend title. Default `TRUE`.
#' @param cells_title_pos Custom position for "Cells" title (overrides centering). Default `NULL`.
#' @param pt_title_pos Custom position for "Pseudotime" title (overrides centering). Default `NULL`.
#' @param hm_title_pos Custom position for heatmap legend title. Default `NULL`.
#'
#' @return `ComplexHeatmap` object invisibly. If `output_pdf` is provided, the heatmap is saved as PNG.
#' @export
#'
#' @examples
#' plot_heatmap_pseudotime(
#'   heatmap_matrices = l1$heatmap_matrices,
#'   seurat_obj = sub,
#'   lineage_name = "Lineage1",
#'   pseudotime_vec = pt_vec,
#'   gene_list = gene_list,
#'   celltype_colors = cols,
#'   font_size = 8,
#'   pdf_width = 4,
#'   pdf_height = 6,
#'   output_pdf = file.path(dir.results),
#'   plot_title = "Neuronal lineage",
#'   legend_dir = "vertical",
#'   ht_legend_side = "right",
#'   ann_legend_side = "right",
#'   annotation_col = "ann2",
#'   merge_legends = TRUE,
#'   center_legend_titles = TRUE,
#'   center_cells_title = TRUE,
#'   center_pt_title = FALSE,
#'   cells_title_pos = "topcenter",
#'   pt_title_pos = "leftcenter",
#'   hm_title_pos = "topcenter"
#' )
plot_heatmap_pseudotime <- function(
    heatmap_matrices,
    seurat_obj,
    lineage_name = NULL,
    pseudotime_vec = NULL,
    gene_list = NULL,
    annotation_col = "ann_level_2",
    celltype_colors = NULL,
    zscore_limits = c(-2, 0, 2),
    pseudotime_palette = "YlGnBu",
    row_cluster_colors = NULL,
    font_size = 10,
    output_pdf = NULL,
    fig_name = paste(plot_title),
    pdf_width = 6,
    pdf_height = 10,
    plot_title = NULL,
    show_annotation_legend = TRUE,
    smooth_df = 3,
    legend_dir = "vertical",
    ht_legend_side = "right",
    ann_legend_side = "right",
    legend_width_cm = 3,
    cont_legend_cm = 3,
    merge_legends = FALSE,
    center_legend_titles = TRUE,
    center_cells_title = TRUE,
    center_pt_title = TRUE,
    cells_title_pos = NULL,
    pt_title_pos = NULL,
    hm_title_pos = NULL,
    file_path = file.path(outdir, fig_name)
){
  require(ComplexHeatmap)
  require(circlize)
  require(RColorBrewer)
  
  # --- Setup lineage and genes
  if(is.null(lineage_name)) lineage_name <- names(heatmap_matrices)[1]
  if(!lineage_name %in% names(heatmap_matrices)) stop("Lineage not found")
  hm <- heatmap_matrices[[lineage_name]]
  if(is.null(gene_list)) gene_list <- hm$top_genes
  cells_hm <- hm$ordered_cells
  expr_mat <- as.matrix(hm$expr_raw[gene_list, cells_hm, drop = FALSE])
  
  # --- Smooth
  if(!is.null(smooth_df) && smooth_df > 0){
    expr_mat <- t(apply(expr_mat, 1, function(x){
      if(sd(x) == 0) return(x)
      tryCatch(smooth.spline(x, df = smooth_df)$y, error = function(e) x)
    }))
  }
  
  # --- Z-score scaling
  row_means <- rowMeans(expr_mat, na.rm = TRUE)
  row_sds <- apply(expr_mat, 1, sd, na.rm = TRUE)
  row_sds[row_sds == 0] <- 1
  expr_mat_scaled <- (expr_mat - row_means) / row_sds
  
  # --- Order genes by pseudotime peak
  gene_peak <- apply(expr_mat_scaled, 1, which.max)
  expr_mat_scaled <- expr_mat_scaled[order(gene_peak), ]
  
  # --- Clusters and annotation
  clusters <- seurat_obj@meta.data[cells_hm, annotation_col]
  if(!is.factor(clusters)){
    cluster_order <- unique(clusters)
    clusters <- factor(clusters, levels = cluster_order)
  } else {
    cluster_order <- levels(clusters)
  }
  clusters[is.na(clusters)] <- "Unknown"
  
  if(is.null(celltype_colors)){
    celltype_colors <- setNames(
      colorRampPalette(brewer.pal(8,"Dark2"))(length(cluster_order)),
      cluster_order
    )
  }
  
  if(is.null(pseudotime_vec)) pseudotime_vec <- seq_along(cells_hm)
  pseudotime <- pseudotime_vec[cells_hm]
  
  pt_col_fun <- colorRamp2(
    seq(min(pseudotime), max(pseudotime), length.out = 100),
    viridis::viridis(100, option = "C")
  )
  
  # --- Helper to resolve title positions
  resolve_title_pos <- function(title_pos, center, is_discrete, legend_dir){
    if(is_discrete){
      if(!is.null(title_pos)) return(title_pos)
      if(legend_dir == "vertical") return("left") else return("top")
    } else {
      if(!is.null(title_pos)) return(title_pos)
      if(center){
        if(legend_dir == "vertical") return("leftcenter-rot") else return("topcenter")
      } else {
        if(legend_dir == "vertical") return("lefttop-rot") else return("topcenter")
      }
    }
  }
  
  is_cells_discrete <- is.factor(clusters)
  is_pt_discrete <- !is.numeric(pseudotime)
  cells_pos <- resolve_title_pos(cells_title_pos, center_cells_title, is_cells_discrete, legend_dir)
  pt_pos <- resolve_title_pos(pt_title_pos, center_pt_title, is_pt_discrete, legend_dir)
  
  top_annot <- HeatmapAnnotation(
    "Cells" = factor(clusters, levels = cluster_order),
    "Pseudotime" = pseudotime,
    col = list("Cells" = celltype_colors, "Pseudotime" = pt_col_fun),
    show_legend = show_annotation_legend,
    annotation_legend_param = list(
      "Cells" = list(
        border = "black",
        labels_gp = gpar(fontsize = font_size),
        title_gp = gpar(
          fontsize = font_size,
          fontface = "bold",
          just = if(center_cells_title) "center" else "left"
        ),
        title_position = cells_pos,
        legend_direction = legend_dir,
        legend_width = unit(legend_width_cm, "cm"),
        legend_height = unit(cont_legend_cm, "cm"),
        nrow = if(legend_dir == "horizontal") 1 else NULL
      ),
      "Pseudotime" = list(
        border = "black",
        labels_gp = gpar(fontsize = font_size),
        title_gp = gpar(
          fontsize = font_size,
          fontface = "bold",
          just = if(center_pt_title) "center" else "left"
        ),
        title_position = pt_pos,
        legend_direction = legend_dir,
        legend_width = unit(legend_width_cm, "cm"),
        legend_height = unit(cont_legend_cm, "cm")
      )
    ),
    annotation_name_gp = gpar(fontsize = font_size, fontface = "bold")
  )
  
  ht <- Heatmap(
    expr_mat_scaled,
    name = "Z-score",
    col = colorRamp2(seq(-2, 2, length = 11), rev(brewer.pal(11, "Spectral"))),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = font_size, fontface="italic"),
    top_annotation = top_annot,
    row_title_rot = 0,
    cluster_row_slices = FALSE,
    show_row_dend = FALSE,
    column_title = plot_title,
    column_title_gp = gpar(fontsize = font_size+2, fontface="bold"),
    heatmap_legend_param = list(
      title = "Z-score",
      legend_direction = legend_dir,
      title_position = resolve_title_pos(hm_title_pos, center_legend_titles, FALSE, legend_dir),
      legend_width = unit(legend_width_cm, "cm"),
      legend_height = unit(cont_legend_cm, "cm"),
      labels_gp = gpar(fontsize = font_size),
      title_gp = gpar(
        fontsize = font_size,
        fontface = "bold",
        just = if(center_legend_titles) "center" else "left"
      ),
      border = "black"
    )
  )
  
  if(!is.null(output_pdf)){
    ComplexHeatmap::draw(ht,
                         heatmap_legend_side = ht_legend_side,
                         annotation_legend_side = ann_legend_side,
                         merge_legends = merge_legends
    )
    save_heatmap(
      ht = ComplexHeatmap::draw(ht,
                                heatmap_legend_side = ht_legend_side,
                                annotation_legend_side = ann_legend_side,
                                merge_legends = TRUE
      ),
      filename = file_path,
      formats = c("png", "pdf", "tiff"),
      width = pdf_width,
      height = pdf_height
    )
  } else {
    ComplexHeatmap::draw(ht,
                         heatmap_legend_side = ht_legend_side,
                         annotation_legend_side = ann_legend_side,
                         merge_legends = TRUE
    )
  }
  
  invisible(ht)
}




#' Create Sankey Plot from Seurat Object
#'
#' This function generates a Sankey plot using the `networkD3` package from a Seurat object,
#' based on selected variables (metadata columns). It visualizes transitions between different 
#' levels or categories across the selected variables.
#'
#' @param object Seurat object containing the metadata.
#' @param selected_vars A character vector of metadata column names to visualize in the Sankey plot. 
#'                       Must contain at least two variables.
#' @param color_list A named list of colors where names correspond to the identity (e.g., cell type) 
#'                   and values correspond to the colors. If NULL, default colors are used.
#' @param sinksRight A boolean value indicating whether the sinks (final nodes) should be positioned on the right 
#'                   side of the plot. Defaults to FALSE.
#' @param fontSize An integer specifying the font size of node labels. Defaults to 13.
#' @param nodeWidth An integer specifying the width of the nodes. Defaults to 40.
#' @param nodePadding An integer specifying the padding between nodes. Defaults to 20.
#'
#' @return A Sankey plot object created using the `networkD3` package, which can be rendered in R Markdown or Jupyter Notebooks.
#'
#' @examples
#' \dontrun{
#' sankey_plot(seurat_obj, selected_vars = c("ann_level_1", "ann_level_2"))
#' }
#'
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @import networkD3
#' @import jsonlite
#' @export
#' 
sankey_plot <- function(seurat_obj, 
                        selected_vars = c("ann_level_1", "ann_level_2", "ann_level_3"),
                        plot.title = NULL,
                        custom_colors = list(),  
                        text.size = 14,
                        show_counts = TRUE,
                        show_percentages = TRUE,
                        show_labels = TRUE,
                        label.justify = "left", 
                        label.nudge = 0.1,
                        flow.alpha = 0.5,
                        x.label_position = "top",
                        custom_x.labels = NULL,  
                        show_x_axis = TRUE,
                        prevent_overlap = TRUE) {    
  
  require(Seurat)
  require(dplyr)
  require(ggplot2)
  require(ggsankey)
  require(viridis)  
  
  # Extract metadata
  metadata <- seurat_obj@meta.data
  
  # Ensure selected variables exist
  missing_cols <- setdiff(selected_vars, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Missing columns in metadata:", paste(missing_cols, collapse = ", ")))
  }
  
  # Convert data to long format
  df <- metadata %>%
    dplyr::select(all_of(selected_vars)) %>%
    make_long(!!!syms(selected_vars))  
  
  # Calculate total per selected variable level
  dagg <- df %>%
    group_by(x, node) %>%
    tally() %>%
    dplyr::mutate(pct = n / sum(n))  
  
  # Merge back into the main dataframe
  df <- left_join(df, dagg, by = c("x", "node"))
  
  # Unique ID for each node
  df$full_node_name <- paste(df$x, df$node, sep = "_")  
  
  # Build color map
  color_map <- c()
  for (var in selected_vars) {
    if (!is.null(custom_colors[[var]])) {
      color_map <- c(color_map, setNames(custom_colors[[var]], paste(var, names(custom_colors[[var]]), sep = "_")))
    }
  }
  
  # # Default to viridis if no custom colors
  # if (length(color_map) == 0) {
  #   viridis_colors <- viridis::viridis(length(unique(df$full_node_name)))
  #   color_map <- setNames(viridis_colors, unique(df$full_node_name))
  # }
  
  # Default to custom_palette if no custom colors
  if (length(color_map) == 0) {
    n_nodes <- length(unique(df$full_node_name))
    palette_colors <- custom_palette(n_nodes, preset = "base")  # choose any preset you like
    color_map <- setNames(palette_colors, unique(df$full_node_name))
  }
  
  # Check missing nodes
  unique_nodes <- unique(df$full_node_name)
  missing_nodes <- setdiff(unique_nodes, names(color_map))
  if (length(missing_nodes) > 0) {
    stop(paste("Error: Missing colors for nodes:", paste(missing_nodes, collapse = ", ")))
  }
  
  # Assign colors
  color_scale <- scale_fill_manual(values = color_map)
  
  #- ADAPTIVE LABEL COLORS-
  get_text_color <- function(hex_color) {
    rgb_vals <- col2rgb(hex_color) / 255
    brightness <- sum(rgb_vals * c(0.299, 0.587, 0.114))  
    if (brightness > 0.5) "black" else "white"
  }
  df$label_color <- sapply(df$full_node_name, function(name) get_text_color(color_map[name]))
  
  #- Labels-
  df$label <- df$node  
  if (show_labels) {
    if (show_counts & show_percentages) {
      df$label <- paste0(df$node, " n=", df$n, " (", round(df$pct * 100, 2), "%)")
    } else if (show_counts) {
      df$label <- paste0(df$node, " n=", df$n)
    } else if (show_percentages) {
      df$label <- paste0(df$node, " (", round(df$pct * 100, 2), "%)")
    }
  } else {
    df$label <- ""  
  }
  
  # Adjust justification
  hjust_value <- ifelse(label.justify == "left", 1, ifelse(label.justify == "center", 0.5, 0))
  nudge_x_value <- ifelse(label.justify == "left", label.nudge, ifelse(label.justify == "right", -label.nudge, 1))
  
  # Prevent overlap by ordering
  if (prevent_overlap) {
    df <- df %>%
      group_by(x) %>%
      arrange(desc(n)) %>%
      mutate(node = factor(node, levels = unique(node)))
  }
  
  #- Plot-
  pl <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, 
                       fill = full_node_name, label = label)) +
    geom_sankey(flow.alpha = flow.alpha, node.color = "black", show.legend = TRUE) +
    geom_sankey_label(aes(color = I(label_color)), size = 3, 
                      hjust = hjust_value, nudge_x = nudge_x_value) +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = margin(1, 1, 1, 1, "cm"),
          axis.text.x = element_text(color = "black", size = text.size, face = "bold"),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    scale_x_discrete(position = x.label_position) +
    color_scale +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.6))) +
    labs(title = plot.title, fill = "Nodes")
  
  # Custom x-axis labels
  if (!is.null(custom_x.labels) & length(custom_x.labels) == length(selected_vars)) {
    pl <- pl + scale_x_discrete(labels = custom_x.labels)
  }
  
  # Hide x-axis if needed
  if (!show_x_axis) {
    pl <- pl + theme(axis.text.x = element_blank())
  }
  
  return(pl)
}


#' Plot Cell Fractions Over Time
#'
#' This function creates various types of visualizations to display the distribution 
#' of cell types over time using Seurat metadata.
#'
#' @param seurat_obj A Seurat object containing cell metadata.
#' @param x.col Character. The column name in metadata representing the x-axis (e.g., "Age").
#' @param fill.col Character. The column name in metadata representing cell types (e.g., "CellType").
#' @param plot.type Character. Type of plot to generate. Choose from:
#'   - "stacked_area" (default)
#'   - "stacked_bar"
#'   - "grouped_bar"
#'   - "point_smooth"
#'   - "heatmap"
#'   - "ridgeline"
#'   - "line"
#' @param text.size Numeric. Text size for axis labels and titles (default: 10).
#' @param legend.title Character. Custom title for the legend (default: NULL).
#' @param legend.position Character. Position of legend ("right", "bottom", "top", "left") (default: "right").
#' @param legend.ncol Numeric. Number of columns in the legend (default: 1).
#' @param custom.colors Character vector. Custom colors for cell types (default: NULL).
#' @param x.title Character. Custom x-axis title (default: same as `x.col`).
#' 
#' @return A ggplot object.
#' @import ggplot2 dplyr RColorBrewer ggridges scales
#' @export
#' 
#' @examples
#' # Example usage:
#' plot_cell_fractions(seurat_obj, plot.type = "ridgeline")
#' plot_cell_fractions(seurat_obj, plot.type = "stacked_area", text.size = 12)
#' 
plot_cell_fractions <- function(seurat_obj, 
                                x.col = "Age", 
                                fill.col = "CellType",
                                plot.type = "stacked_area",
                                text.size = 12,
                                pt.size = 2,
                                line.size = 0.5,
                                legend.title = NULL,
                                legend.position = "right",
                                legend.ncol = 1,
                                custom.colors = NULL,
                                x.title = x.col
) {
  
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(ggridges)
  library(scales)
  
  # Ensure metadata contains the required columns
  metadata <- seurat_obj@meta.data %>%
    dplyr::select(all_of(c(x.col, fill.col))) %>%
    dplyr::group_by(!!sym(x.col), !!sym(fill.col)) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(!!sym(x.col)) %>%
    dplyr::mutate(fraction = count / sum(count)) %>%
    dplyr::ungroup()
  
  # Ensure x.col is numeric
  metadata[[x.col]] <- as.numeric(as.character(metadata[[x.col]]))
  
  # Define colors
  num_categories <- length(unique(metadata[[fill.col]]))
  if (is.null(custom.colors)) {
    custom.colors <- colorRampPalette(RColorBrewer::brewer.pal(min(9, num_categories), "Set3"))(num_categories)
  }
  
  # Create the plot based on plot type
  if (plot.type == "stacked_area") {
    p <- ggplot(metadata, aes(x = !!sym(x.col), y = fraction, fill = !!sym(fill.col))) +
      geom_area(position = "stack", alpha = 1, colour = "black", size = 0.2) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = custom.colors, name = legend.title) +
      labs(x = x.title, y = "Fraction of cells") +
      theme_classic()
  }
  
  else if (plot.type == "stacked_bar") {
    p <- ggplot(metadata, aes(x = !!sym(x.col), y = fraction, fill = !!sym(fill.col))) +
      geom_bar(stat = "identity", position = "fill", colour = "black", size = 0.2) +
      scale_y_continuous(labels = scales::percent) +
      scale_fill_manual(values = custom.colors, name = legend.title) +
      labs(x = x.title, y = "Fraction of cells") +
      theme_classic()
  }
  
  else if (plot.type == "grouped_bar") {
    p <- ggplot(metadata, aes(x = !!sym(x.col), y = count, fill = !!sym(fill.col))) +
      geom_bar(stat = "identity", position = "dodge", colour = "black", size = 0.2) +
      scale_fill_manual(values = custom.colors, name = legend.title) +
      labs(x = x.title, y = "Cell Count") +
      theme_classic()
  }
  
  else if (plot.type == "point_smooth") {
    p <- ggplot(metadata, aes(x = !!sym(x.col), y = fraction, color = !!sym(fill.col))) +
      geom_point(size = pt.size) +
      geom_smooth(method = "loess", se = FALSE, size = line.size) +
      scale_y_continuous(labels = scales::percent) +
      scale_color_manual(values = custom.colors, name = legend.title) +
      labs(x = x.title, y = "Fraction of cells") +
      theme_classic()
  }
  
  else if (plot.type == "heatmap") {
    p <- ggplot(metadata, aes(x = !!sym(x.col), y = !!sym(fill.col), fill = fraction)) +
      geom_tile() +
      viridis::scale_fill_viridis_c(option = "plasma") +
      labs(x = x.title, y = "Cell type", fill = "Fraction of cells") +
      theme_minimal()
  }
  
  else if (plot.type == "ridgeline") {
    p <- ggplot(metadata, aes(x = !!sym(x.col), y = !!sym(fill.col), height = fraction, fill = !!sym(fill.col))) +
      geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
      scale_fill_manual(values = custom.colors, name = legend.title) +
      labs(x = x.title, y = "Cell type", fill = "Fraction of cells") +
      theme_classic()
  } else if (plot.type == "line") {
    p <- ggplot(metadata, aes(x = !!sym(x.col), y = fraction, color = !!sym(fill.col), group = !!sym(fill.col))) +
      geom_line(size = line.size) + 
      geom_point(size = pt.size) +
      scale_y_continuous(labels = scales::percent) +
      scale_color_manual(values = custom.colors, name = legend.title) +
      labs(x = x.title, y = "Fraction of cells") +
      theme_classic()
  } else {
    stop("Invalid plot type. Choose from: stacked_area, stacked_bar, grouped_bar, point_smooth, heatmap, line, ridgeline.")
  }
  
  # Adjust theme settings
  p <- p +
    theme(
      plot.title = element_text(size = text.size + 2, face = "bold", hjust = 0.5),
      axis.title = element_text(size = text.size),
      axis.text = element_text(size = text.size, color = "black"),
      legend.key.size = unit(0.4, 'cm'),
      legend.text = element_text(size = text.size-1),
      legend.position = legend.position,
      legend.spacing.y = unit(0.05, 'cm'),  # Reduce space between legend items
      legend.margin = margin(t = -10, b = 0, unit = 'pt') # Reduce space between plot and legend
    ) +
    guides(fill = guide_legend(
      ncol = legend.ncol, 
      override.aes = list(color = "black", size = 0.5) 
    ))
  
  return(p)
}


#' Cell–cell (cluster) expression correlation heatmap
#'
#' Compute average gene expression per group (e.g. cluster, cell type)
#' and visualize pairwise Pearson correlations using a ComplexHeatmap.
#'
#' @param object A \code{Seurat} object.
#' @param group.by Metadata column used to define groups (default: \code{"seurat_clusters"}).
#' @param assay Assay to use (default: \code{"RNA"}).
#' @param slot Expression slot/layer to use (e.g. \code{"data"}, \code{"counts"}, \code{"scale.data"}).
#' @param plot.title Optional title displayed above the heatmap.
#' @param column_names_rot Rotation angle of column labels (default: 45).
#' @param fontsize Base font size for row/column labels.
#' @param pdf_width Width of saved figure (in inches).
#' @param pdf_height Height of saved figure (in inches).
#' @param fig_name File name (without extension) used when \code{save = TRUE}.
#' @param save Logical; whether to save the heatmap to disk.
#' @param dir.results Output directory for saved figures.
#'
#' @details
#' For each group defined by \code{group.by}, the function computes
#' mean gene expression across all cells in that group and calculates
#' a Pearson correlation matrix between groups.
#'
#' Compatible with Seurat v4 and v5.
#'
#' @return
#' A \code{ComplexHeatmap::Heatmap} object (invisibly).
#'
#' @examples
#' \dontrun{
#' cellcorr(
#'   object = seurat_obj,
#'   group.by = "ann0",
#'   plot.title = "Cluster correlation",
#'   column_names_rot = 90,
#'   save = TRUE,
#'   fig_name = "cluster_corr",
#'   dir.results = "figures"
#' )
#' }
#'
#' @import Seurat
#' @import ComplexHeatmap
#' @importFrom Matrix rowMeans
#' @importFrom grid unit gpar
#' @importFrom stats cor
#'
#' @export
cellcorr <- function(object, 
                     group.by = "seurat_clusters",
                     assay = "RNA", 
                     slot = "data",
                     plot.title = NULL,
                     column_names_rot = 45,
                     fontsize = 10,
                     pdf_width = 5,
                     pdf_height = 6,
                     fig_name = NULL,
                     save = FALSE,
                     dir.results = ".") {
  require(Seurat)
  require(ComplexHeatmap)
  require(RColorBrewer)
  
  # Set identity class based on the group.by column
  if (!group.by %in% colnames(object@meta.data)) {
    stop(paste("Column", group.by, "not found in metadata."))
  }
  
  Idents(object) <- object@meta.data[[group.by]]
  
  # Extract unique cluster identities
  uniq <- unique(Idents(object))
  
  # Extract expression data
  if (packageVersion("Seurat") >= "5.0.0") {
    expr_data <- object[[assay]][slot] # Seurat v5 fix
  } else {
    expr_data <- GetAssayData(object, slot = slot, assay = assay)  # Seurat v4
  }
  
  # Initialize an empty matrix with correct dimensions
  mnmat <- matrix(NA, nrow = nrow(expr_data), ncol = length(uniq))
  
  # Iterate over each cluster and compute the average expression
  for (i in seq_along(uniq)) {
    cluster_cells <- WhichCells(object, idents = uniq[i])
    
    # Check if the cluster contains valid cells
    valid_cells <- cluster_cells[cluster_cells %in% colnames(expr_data)]
    
    if (length(valid_cells) > 0) {
      mnmat[, i] <- Matrix::rowMeans(expr_data[, valid_cells, drop = FALSE])
    }
  }
  
  colnames(mnmat) <- as.vector(uniq)
  
  # Compute Pearson correlation matrix
  mat <- cor(mnmat, method = "pearson", use = "complete.obs")
  
  # Generate heatmap
  ht <- Heatmap(
    mat,
    name = 'Pearson correlation',
    #col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), # Red-Blue gradient
    col = viridis::viridis(100),
    border = "#8B8B8C",
    rect_gp = gpar(col = "#8B8B8C"),
    row_names_gp = gpar(fontsize = fontsize),  
    column_names_gp = gpar(fontsize = fontsize, fontface = "plain"),
    column_names_side = "bottom",
    column_names_rot = column_names_rot,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    show_column_dend = TRUE,
    show_row_dend = TRUE,
    heatmap_legend_param = list(
      title = 'Pearson correlation',
      title_position = 'leftcenter-rot',
      legend_height = unit(3.5, 'cm'),
      legend_width = unit(0.5, 'cm'),
      border = '#8B8B8C',
      #legend.direction = 'horizontal',
      labels_gp = gpar(fontsize = fontsize)
    )
  )
  
  if (save) {
    save_heatmap(
      ht = ComplexHeatmap::draw(
        ht,
        heatmap_legend_side = "right",
        annotation_legend_side = "right"
      ),
      filename = file.path(dir.results, fig_name),
      formats = c("pdf", "png", "tiff"),
      width = pdf_width,
      height = pdf_height
    )
  }
  
  # Draw heatmap with title
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  grid::grid.text(
    plot.title,
    x = 0.5, y = unit(1, "npc") - unit(2, "mm"),
    gp = gpar(fontsize = fontsize + 2, fontface = "bold")
  )
}



#' Plot Marker Positivity Percentage Before/After DecontX
#'
#' This function reproduces the logic of
#' \code{plotDecontXMarkerPercentage()} but allows you to group cells
#' using any \code{colData(sce)} column (e.g., \code{ann_level_2})
#' rather than being restricted to \code{decontX_clusters}.
#'
#' It computes the percentage of cells expressing one or more markers
#' (above a given threshold) and visualizes these percentages for one
#' or multiple assays (e.g. raw counts vs DecontX-cleaned counts).
#'
#' @param sce A \code{SingleCellExperiment} object containing counts and,
#' optionally, DecontX-cleaned counts.
#' @param markers A character vector of gene names to evaluate. Genes
#' not present in \code{rownames(sce)} are dropped automatically.
#' @param group_by Column name in \code{colData(sce)} to group cells by
#' (e.g., cell type annotation).
#' @param assays Character vector of assay names to compare, typically
#' \code{c("counts", "decontXcounts")}.
#' @param expr_threshold Numeric expression threshold. Cells with
#' \code{expr > expr_threshold} are considered positive for that marker.
#' Default = 1.
#' @param plot_type Either \code{"bar"} (default) or \code{"heatmap"}.
#' @details
#' For each marker and each group, the percentage of positive cells is:
#'
#' \deqn{
#' 100 * \frac{\sum(expr > threshold)}{N_{\text{cells in group}}}
#' }
#'
#' @return A list with:
#' \describe{
#'   \item{table}{A data.frame of marker positivity percentages.}
#'   \item{plot}{A ggplot2 object showing the result.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- plot_decontx_marker_percent(
#'   sce = clean$sce,
#'   markers = c("XIST","UTY","HBB"),
#'   group_by = "ann_level_2",
#'   assays = c("counts","decontXcounts"),
#'   plot_type = "bar"
#' )
#' res$plot
#' }
#'
#' @export
#' 
plot_decontx_marker_percent <- function(
    sce,
    markers,
    group_by = "ann2",
    assays = c("counts", "decontXcounts"),
    expr_threshold = 1,
    plot_type = c("bar","heatmap"),
    theme = "bw",
    fontsize = 10,
    x.angle = 45,
    x.lab = NULL,
    ncol = NULL,
    leg.pos = "right",
    leg.dir = "vertical",
    custom_colors = NULL,
    ...
){
  require(SingleCellExperiment)
  require(Matrix)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  
  plot_type <- match.arg(plot_type)
  
  if (!group_by %in% colnames(colData(sce))) {
    stop("group_by must be a column in colData(sce)")
  }
  
  # Check markers exist
  markers <- intersect(markers, rownames(sce))
  if (length(markers) == 0) stop("No markers found in SCE")
  
  # Build results table
  res_list <- list()
  
  for (assay_name in assays) {
    if (!assay_name %in% assayNames(sce))
      stop(paste("Assay", assay_name, "not found in sce"))
    
    mat <- assay(sce, assay_name)
    
    df <- lapply(markers, function(m){
      expr <- mat[m, ]
      is_pos <- as.numeric(expr > expr_threshold)
      
      data.frame(
        cell   = colnames(sce),
        group  = colData(sce)[[group_by]],
        marker = m,
        pos    = is_pos
      )
    }) %>% bind_rows()
    
    perc <- df %>%
      group_by(group, marker) %>%
      summarise(
        pct = mean(pos) * 100,
        .groups = "drop"
      ) %>%
      mutate(assay = assay_name)
    
    res_list[[assay_name]] <- perc
  }
  
  out <- bind_rows(res_list)
  
  # If custom_colors is NULL, assign default
  if (is.null(custom_colors)) {
    custom_colors <- c("counts" = "#E69F00", "decontXcounts" = "#56B4E9")
  }
  
  # PLOTS
  if (plot_type == "bar") {
    # BAR
    p <- ggplot(out, aes(x = group, y = pct, fill = assay)) +
      geom_bar(stat = "identity", position = "dodge",
               color = "black", linewidth = 0.2) +
      facet_wrap(~ marker, scales = "free_y", ncol = ncol) +
      plot_theme(theme = theme, font.size = fontsize, x.angle = x.angle, 
                 leg.pos = leg.pos, leg.dir = leg.dir,...) +
      labs(y = "% cells expressing", x = x.lab) +
      scale_fill_manual(values = custom_colors)
    
  } else if (plot_type == "heatmap") {
    # HEATMAP
    p <- ggplot(out, aes(x = marker, y = group, fill = pct)) +
      geom_tile() +
      facet_wrap(~ assay, ncol = ncol) +
      scale_fill_viridis_c() +
      plot_theme(theme = theme, font.size = fontsize, x.angle = x.angle, 
                 leg.pos = leg.pos, leg.dir = leg.dir,...) +
      labs(fill = "% positive", x = "Marker", y = group_by)
  }
  
  return(list(
    table = out,
    plot  = p
  ))
}



#' Plot Sex Chromosome QC
#'
#' Visualizes expression of X and Y chromosome marker genes per sample to assess sex
#' contamination and discordance in single-cell RNA-seq data. Optionally highlights
#' discordant cells flagged by `sex_contamination_qc()`.
#'
#' @param obj A \code{Seurat} object with RNA assay and metadata containing \code{sex_call} 
#' and optionally \code{sex_discordant} columns (as generated by `sex_contamination_qc()`).
#' @param genes Character vector of marker genes to plot. Default: c("XIST","UTY").
#' @param sample_col Character. Column in \code{obj@meta.data} specifying sample identity. Default: "sample".
#' @param highlight_discordant Logical. If TRUE, discordant cells are highlighted in black. Default: TRUE.
#' @param pt_size Numeric. Size of points in jitter overlay. Default: 0.5.
#' @param fontsize Numeric. Base font size for the plot. Default: 14.
#' @param palette Named character vector of colors for sex calls. Must include "female_like", "male_like", and "ambiguous". Default: c(female_like="#b15928", male_like="#6a3d9a","ambiguous"="gray70").
#'
#' @return A ggplot object displaying violin plots of normalized expression per sample for each gene.
#'
#' @details
#' This function generates violin plots for key sex chromosome markers (e.g., XIST, UTY)
#' across samples. Cells can be colored by sex call and discordant cells highlighted.
#' Useful for visual QC of sex assignment after sex contamination analysis.
#'
#' @examples
#' \dontrun{
#' qc_res <- sex_contamination_qc(obj)
#' p <- plot_sex_qc(qc_res$object)
#' print(p)
#' }
#'
#' @export
plot_sex_qc <- function(
    obj,
    genes = c("XIST","UTY"),
    sample_col = "sample",
    highlight_discordant = TRUE,
    pt_size = 0.5,
    fontsize = 14,
    palette = c("Female-like"="#b15928"," Male-like"="#6a3d9a","Ambiguous"="gray70")
){
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Ensure RNA assay
  DefaultAssay(obj) <- "RNA"
  expr <- GetAssayData(obj, slot="data")[genes, , drop=FALSE] %>% t() %>% as.data.frame()
  expr$cell <- rownames(expr)
  
  # Add metadata
  meta <- obj@meta.data %>% 
    dplyr::select(all_of(sample_col), sex_call, sex_discordant) %>%
    dplyr::mutate(cell = rownames(.))
  
  df <- expr %>% 
    tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expression") %>%
    dplyr::left_join(meta, by="cell")
  
  # Discordant coloring
  if(highlight_discordant & "sex_discordant" %in% colnames(df)){
    df$color_group <- ifelse(df$sex_discordant, "discordant", df$sex_call)
    palette <- c(palette, discordant="black")
  } else {
    df$color_group <- df$sex_call
  }
  
  p <- ggplot(df, aes(x=.data[[sample_col]], y=expression, fill=color_group)) +
    geom_violin(scale="width", trim=TRUE, adjust=1) +
    geom_jitter(width=0.1, size=pt_size, alpha=0.5) +
    facet_wrap(~gene, scales="free_y") +
    scale_fill_manual(values=palette, name="Sex Call") +
    plot_theme(theme = "classic", font.size=fontsize, x.angle = 45, legend.position = "top", legend.direction = "horizontal") +
    labs(x="Sample", y="Normalized expression", title="Sex Chromosome QC")
  
  return(p)
}




#' Generate QC plots for sex and erythroid contamination
#'
#' @param meta_df Metadata data frame from Seurat object (with sex_score, sex_call, eryth_sum)
#' @param tbl_summary Summary table per sample (n_cells, n_removed, median_eryth_sum)
#' @param theme ggplot theme type
#' @param font.size Base font size
#' @param x.angle Angle of x-axis labels
#' @param out_dir Optional directory to save plots
#'
#' @return List of ggplot objects
#' @export
plot_sex_eryth_qc <- function(meta_df, tbl_summary,
                              eryth_genes = c("HBB","HBA1","HBA2","ALAS2"),  # default markers
                              theme = "classic",
                              font.size = 8,
                              x.angle = 45,
                              out_dir = NULL) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  meta_df = obj@meta.data
  # Ensure eryth_sum exists
  # Ensure eryth_sum exists and is numeric
  if(!"eryth_sum" %in% colnames(meta_df)) {
    available_genes <- intersect(eryth_genes, colnames(meta_df))
    if(length(available_genes) > 0) {
      meta_df$eryth_sum <- rowSums(meta_df[, available_genes, drop=FALSE])
    } else {
      meta_df$eryth_sum <- NA
    }
  }
  # Convert to numeric in case it is factor/character
  meta_df$eryth_sum <- as.numeric(meta_df$eryth_sum)
  
  
  # 1) Sex classification fraction barplot
  df_sex_long <- tbl_summary %>%
    dplyr::select(sample, frac_female_like, frac_male_like, frac_ambiguous) %>%
    tidyr::pivot_longer(
      cols = starts_with("frac_"),
      names_to = "sex_class",
      values_to = "fraction"
    ) %>%
    mutate(sex_class = recode(sex_class,
                              frac_female_like="Female-like",
                              frac_male_like="Male-like",
                              frac_ambiguous="Ambiguous"))
  
  p1 <- ggplot(df_sex_long, aes(x=sample, y=fraction, fill=sex_class)) +
    geom_col(width=0.8, color="black", linewidth=0.2) +
    scale_fill_manual(values=c("Male-like" = "#6a3d9a", "Female-like" = "#b15928","Ambiguous" = "gray60"),
                      breaks = c("Female-like","Male-like","Ambiguous")) +
    plot_theme(theme = theme, x.angle=x.angle, font.size = font.size) +
    labs(#title = "Sex classification per sample",
      x = "", y = "Fraction of cells", fill = "Sex call")
  
  # 2) Sex score histogram
  p2 <- ggplot(meta_df, aes(x=sex_score, fill=sex_call)) +
    geom_histogram(bins=80, alpha=1) +
    scale_fill_manual(values=c("Male-like" = "#6a3d9a", "Female-like" = "#b15928","Ambiguous" = "gray60")) +
    plot_theme(theme = theme, font.size=font.size, x.angle=0, leg.pos = c(0.2, 0.85)) +
    labs(x="Log2(XIST / Y genes)", y="Number of cells", fill="Sex call")
  
  # 3) Erythroid gene sum histogram
  p3 <- ggplot(meta_df, aes(x=eryth_sum)) +
    geom_histogram(bins=80, fill="steelblue", alpha=1) +
    plot_theme(theme = theme, font.size = font.size) +
    labs(x="Erythroid gene sum", y="Number of cells")
  
  # 4) Median erythroid sum per sample
  p4 <- ggplot(tbl_summary, aes(x=sample, y=median_eryth_sum)) +
    geom_col(fill="tomato", color="black", linewidth=0.2) +
    plot_theme(theme = theme, font.size = font.size, x.angle=x.angle) +
    labs(x="Sample", y="Median erythroid sum")
  
  # 5) Cells before/after QC per sample (dodge bars)
  tbl2 <- tbl_summary %>%
    dplyr::mutate(n_before = n_cells,
                  n_after  = n_before - n_removed) %>%
    dplyr::select(sample, n_before, n_after) %>%
    tidyr::pivot_longer(cols = c(n_before, n_after),
                        names_to = "stage",
                        values_to = "n_cells") %>%
    dplyr::mutate(stage = recode(stage, n_before="Pre-sex-QC", n_after="Post-sex-QC"),
                  stage = factor(stage, levels = c("Pre-sex-QC","Post-sex-QC")))  # enforce order
  
  p5 <- ggplot(tbl2, aes(x=sample, y=n_cells, fill=stage)) +
    geom_col(position = position_dodge(width = 0.6), color="black", width=0.6, linewidth=0.2) +
    scale_fill_manual(values = c("Pre-sex-QC"="#151A1F","Post-sex-QC"="#8CB9E6"),
                      breaks = c("Pre-sex-QC","Post-sex-QC")) +
    plot_theme(theme = theme, font.size=font.size, x.angle=x.angle, leg.pos = c(0.2, 0.9)) +
    labs(x="Sample", y="Number of cells", fill="Stage")
  
  # Save plots if out_dir is provided
  if(!is.null(out_dir)) {
    if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    ggsave(file.path(out_dir, "Sex_class.pdf"), p1, width=5, height=4)
    ggsave(file.path(out_dir, "sex_score_hist.pdf"), p2, width=4, height=4)
    ggsave(file.path(out_dir, "eryth_sum_hist.pdf"), p3, width=4, height=4)
    ggsave(file.path(out_dir, "median_eryth_per_sample.pdf"), p4, width=4, height=4)
    ggsave(file.path(out_dir, "sex_cells_before_after_qc.pdf"), p5, width=4, height=4)
  }
  
  return(list(
    sex_class_fraction = p1,
    sex_score_hist = p2,
    eryth_sum_hist = p3,
    median_eryth = p4,
    cells_before_after = p5
  ))
}



#' Plot Slingshot Lineages on a Dimensional Reduction (UMAP/TSNE/PCA)
#'
#' @description
#' This function visualizes Slingshot lineages on top of a dimensionality 
#' reduction (e.g., UMAP). Cells can be colored by cluster identity or pseudotime.
#' The function supports plotting a subset of lineages, faceting by lineage,
#' customizing color palettes, and adding cluster/lineage labels.
#'
#' @param obj A Seurat object containing the specified reduction.
#' @param sds A \code{SlingshotDataSet} object containing lineage and 
#'   pseudotime information.
#' @param cluster_col Metadata column in \code{obj} used for cluster coloring.
#' @param reduction Dimensional reduction to use (default: \code{"umap"}).
#' @param dims Vector of two integers specifying which dimensions to plot.
#' @param pt.size Point size for cell embeddings.
#' @param cluster_palette Optional named vector of colors for clusters.
#' @param lineage_colors Optional vector of colors for individual lineages.
#'   If NULL, defaults will be generated.
#' @param base.size Base font size for the plot theme.
#' @param show_cluster_labels Logical; if TRUE, cluster text labels are shown.
#' @param show_lineage_labels Logical; if TRUE, labels per lineage are shown.
#' @param label.size Text size for cluster/lineage labels.
#' @param label.fontface Fontface for text labels (e.g., \code{"bold"}).
#' @param legend.position Position of ggplot legend (e.g., \code{"right"}).
#' @param plot.title Optional title for the plot.
#' @param fig.plot Logical; if TRUE, the plot is displayed. If FALSE, it is returned.
#' @param font.size Global font size for text elements.
#' @param color.by One of \code{"cluster"} or \code{"pseudotime"} indicating
#'   how cells should be colored.
#' @param lineage.select Optional vector specifying which lineages to plot
#'   (e.g., \code{c("Lineage1", "Lineage3")}).
#' @param facet.by.lineage Logical; if TRUE, each lineage is plotted in a facet.
#' @param pt.lim Optional numeric vector of length 2 specifying pseudotime limits
#'   for color scaling when \code{color.by = "pseudotime"}.
#'
#' @return A ggplot2 object. If \code{fig.plot = TRUE}, the plot is also displayed.
#'
#'#' @examples
#' \dontrun{
#' library(Seurat)
#' library(slingshot)
#'
#' # Example dataset
#' data("pbmc_small")
#' obj <- pbmc_small
#'
#' # Run PCA and UMAP
#' obj <- RunPCA(obj, verbose = FALSE)
#' obj <- RunUMAP(obj, dims = 1:10, verbose = FALSE)
#'
#' # Convert to SingleCellExperiment for Slingshot
#' sce <- as.SingleCellExperiment(obj)
#'
#' # Run Slingshot on PCA embedding
#' sce <- slingshot(
#'   sce,
#'   clusterLabels = "RNA_snn_res.1",
#'   reducedDim = "PCA"
#' )
#'
#' # Extract SlingshotDataSet
#' sds <- SlingshotDataSet(sce)
#'
#' # Basic usage
#' plot_lineages(
#'   obj = obj,
#'   sds = sds,
#'   plot.title = "Slingshot lineages on UMAP"
#' )
#'
#' # Color cells by pseudotime
#' plot_lineages(
#'   obj = obj,
#'   sds = sds,
#'   color.by = "pseudotime",
#'   lineage.select = "Lineage1"
#' )
#'
#' # Show cluster and lineage labels
#' plot_lineages(
#'   obj = obj,
#'   sds = sds,
#'   show_cluster_labels = TRUE,
#'   show_lineage_labels = TRUE
#' )
#'
#' # Facet each lineage separately
#' plot_lineages(
#'   obj = obj,
#'   sds = sds,
#'   facet.by.lineage = TRUE
#' )
#' }
#' @export
#' 
plot_lineages <- function(
    obj,
    sds,
    cluster_col = "ann_level_2",
    reduction = "umap",
    dims = c(1,2),
    pt.size = 0.5,
    cluster_palette = NULL,
    lineage_colors = NULL,
    base.size = 10,
    show_cluster_labels = FALSE,
    show_lineage_labels = FALSE,
    show_lineages = TRUE,
    lineage_rename = NULL,
    label.size = 4,
    theme = "classic",
    label.fontface = "bold",
    legend.position = "right",
    plot.title = NULL,
    fig.plot = TRUE,
    no.axes = FALSE,
    font.size = 8,
    alpha = 0.8,
    linewidth = 1,
    color.by = c("cluster", "pseudotime"),
    lineage.select = NULL,
    facet.by.lineage = FALSE,
    pt.lim = NULL
) {
  
  require(ggplot2)
  require(dplyr)
  require(mgcv)
  require(ggrepel)
  require(patchwork)
  require(viridis)
  
  color.by <- match.arg(color.by)
  
  # Base UMAP coordinates
  umap.df <- as.data.frame(Embeddings(obj, reduction))
  umap.df$cluster <- obj[[cluster_col, drop = TRUE]]
  colnames(umap.df)[dims[1]] <- "UMAP1"
  colnames(umap.df)[dims[2]] <- "UMAP2"
  
  # Handle pseudotime matrix
  if (!is.null(sds)) {
    pt.mat <- slingPseudotime(sds)
    lineage.names <- colnames(pt.mat)
    
    if (!is.null(lineage.select)) {
      pt.mat <- pt.mat[, lineage.select, drop = FALSE]
      lineage.names <- colnames(pt.mat)
    }
    
    # Short lineage labels
    lineage_label_short <- paste0("Lin", seq_len(ncol(pt.mat)))
    names(lineage_label_short) <- lineage.names
  }
  
  # Short lineage labels
  lineage_label_short <- paste0("Lin", seq_len(ncol(pt.mat)))
  names(lineage_label_short) <- lineage.names
  
  # --- NEW: Optional user-specified renaming ---
  if (!is.null(lineage_rename)) {
    # lineage_rename should be a *named vector*: c(old = "new")
    valid_rename <- intersect(names(lineage_rename), lineage.names)
    if (length(valid_rename) > 0) {
      lineage_label_short[valid_rename] <- lineage_rename[valid_rename]
    }
  }
  
  # Color points
  if (color.by == "cluster") {
    if (is.null(cluster_palette)) {
      cluster_palette <- RColorBrewer::brewer.pal(
        min(12, length(unique(umap.df$cluster))), "Set3"
      )
      names(cluster_palette) <- unique(umap.df$cluster)
    }
    umap.df$color_value <- umap.df$cluster
    color_scale <- scale_color_manual(values = cluster_palette)
    
  } else if (color.by == "pseudotime") {
    # If no specific lineage is selected, combine all pseudotimes
    if (is.null(lineage.select)) {
      # Average pseudotime across lineages (like in plot_slingshot_umap)
      umap.df$color_value <- rowMeans(pt.mat, na.rm = TRUE)
      pt.lim <- range(umap.df$color_value, na.rm = TRUE)
      
    } else {
      # Existing behavior for a single lineage
      umap.df$color_value <- pt.mat[, lineage.select]
      if (is.null(pt.lim)) pt.lim <- range(umap.df$color_value, na.rm = TRUE)
    }
    
    color_scale <- scale_color_viridis(
      option = "C",
      limits = pt.lim,
      na.value = "grey80",
      name = "Pseudotime",
      guide = guide_colorbar(
        direction = "vertical",
        title.position = "left",
        title.hjust = 0.5,
        frame.colour = "black",
        frame.linewidth = 0.2,
        title.theme = element_text(angle = 90, vjust = 0.5),
        barheight = 5,
        barwidth = 0.8
      )
    )
    
    legend.position <- "right"
  }
  
  # Faceting by lineage
  if (facet.by.lineage && !is.null(sds)) {
    facet.df <- umap.df
    facet.df$lineage <- NA
    for (i in seq_len(ncol(pt.mat))) {
      pts <- pt.mat[, i]
      facet.df$lineage[!is.na(pts)] <- lineage_label_short[lineage.names[i]]
    }
    p <- ggplot(facet.df, aes(UMAP1, UMAP2, color = color_value)) +
      geom_point(size = pt.size, alpha = alpha) +
      color_scale +
      theme_void(base_size = font.size) +
      facet_wrap(~ lineage, ncol = 2) +
      theme(
        strip.text = element_text(size = base.size, face = "bold"),
        legend.position = legend.position
      ) +
      labs(title = plot.title)
    return(p)
  }
  
  # Base plot (points)
  # plt <- ggplot(umap.df, aes(x = UMAP1, y = UMAP2, color = color_value)) +
  #   geom_point(size = pt.size, alpha = alpha) +
  #   color_scale +
  #   plot_theme(theme = theme) +
  #   labs(title = plot.title) +
  #   theme(legend.position = legend.position)
  
  plt <- ggplot(umap.df, aes(x = UMAP1, y = UMAP2, color = color_value)) +
    geom_point(size = pt.size, alpha = alpha) +
    color_scale +
    plot_theme(theme = theme) +
    labs(title = plot.title) +
    theme(
      legend.position = legend.position,
      plot.title = element_text(size = base.size + 2, face = "bold")
    )
  
  # Lineages
  if (!is.null(sds) && show_lineages) {
    if (is.null(lineage_colors)) lineage_colors <- viridis::viridis(ncol(pt.mat))
    
    for (i in seq_len(ncol(pt.mat))) {
      df_lineage <- umap.df
      df_lineage$pt <- pt.mat[, i]
      df_lineage <- df_lineage %>% filter(!is.na(pt))
      
      if (nrow(df_lineage) > 20) {
        g1 <- mgcv::gam(UMAP1 ~ s(pt, k = 20), data = df_lineage)
        g2 <- mgcv::gam(UMAP2 ~ s(pt, k = 20), data = df_lineage)
        grid <- seq(min(df_lineage$pt), max(df_lineage$pt), length.out = 200)
        
        curve_df <- data.frame(
          UMAP1 = predict(g1, newdata = data.frame(pt = grid)),
          UMAP2 = predict(g2, newdata = data.frame(pt = grid))
        )
        
        # Lineage curves behind points
        plt <- plt + geom_path(
          data = curve_df,
          aes(UMAP1, UMAP2),
          color = lineage_colors[i],
          linewidth = linewidth,
          inherit.aes = FALSE
        )
        
        # Lineage labels
        if (show_lineage_labels) {
          end.pt <- curve_df[nrow(curve_df), , drop = FALSE]
          end.pt$label <- lineage_label_short[lineage.names[i]]
          plt <- plt + ggrepel::geom_text_repel(
            data = end.pt,
            aes(UMAP1, UMAP2, label = label),
            color = lineage_colors[i],
            fontface = label.fontface,
            size = label.size,
            bg.color = "grey95",
            bg.r = 0.1,
            seed = 42
          )
        }
      }
    }
  }
  
  # Cluster labels on top
  if (show_cluster_labels && color.by == "cluster") {
    centroids <- umap.df %>% group_by(cluster) %>%
      summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))
    plt <- plt + ggrepel::geom_text_repel(
      data = centroids,
      aes(UMAP1, UMAP2, label = cluster),
      size = label.size,
      color = "black",
      fontface = label.fontface,
      bg.color = "grey95",
      bg.r = 0.1,
      seed = 42
    )
  }
  
  if (no.axes) plt <- plt & NoAxes()
  
  # Optional axis arrows
  if (fig.plot) {
    x.lab.reduc <- plt$labels$x %||% paste0(toupper(reduction), dims[1])
    y.lab.reduc <- plt$labels$y %||% paste0(toupper(reduction), dims[2])
    plt <- plt & NoAxes()
    
    L <- 0.12
    axis.df <- data.frame(x0 = c(0,0), y0 = c(0,0), x1 = c(L,0), y1 = c(0,L))
    axis.plot <- ggplot(axis.df) +
      geom_segment(aes(x=x0,y=y0,xend=x1,yend=y1),
                   linewidth=0.6,lineend="round") +
      xlab(x.lab.reduc) + ylab(y.lab.reduc) +
      coord_fixed() + theme_classic(base_size=font.size) +
      theme(plot.background=element_rect(fill="transparent",colour=NA),
            panel.background=element_rect(fill="transparent",colour=NA),
            axis.text=element_blank(), axis.ticks=element_blank(),
            axis.line=element_blank(), panel.border=element_blank(),
            axis.title=element_text(size=font.size, face="bold"),
            plot.margin=margin(0,0,0,0))
    
    figure.layout <- c(
      patchwork::area(t = 1, l = 1, b = 11, r = 11),
      patchwork::area(t = 10, l = 1, b = 11, r = 2)
    )
    
    return(plt + axis.plot + patchwork::plot_layout(design = figure.layout))
  }
  
  return(plt)
}



#' Gene Set Expression Plot
#'
#' @description
#' Visualize the expression of a set of genes across cell types and stages in a Seurat object.
#' Supports three modes:
#' \itemize{
#'   \item "standard": dodged bar plot with proportion ± SEM per cell type
#'   \item "stacked": stacked bar plot of gene expression bins per cell type
#'   \item "faceted": faceted bar plot by stage with optional mean ± SEM overlay
#' }
#' 
#' @param obj Seurat object containing RNA expression data.
#' @param geneset Character vector of gene names to plot.
#' @param group.by Column name in `obj@meta.data` for grouping cells (default: "ann2").
#' @param stage Column name in `obj@meta.data` for stage/facet grouping (default: "stage").
#' @param mode Plot mode: "standard", "stacked", or "faceted". Default is "standard".
#' @param assay Assay in Seurat object to use (default: "RNA").
#' @param slot Slot to extract data from (default: "data").
#' @param expr.thresh Numeric, threshold for considering a gene expressed (default: 0).
#' @param breaks Numeric vector of bin edges for counting genes expressed (default: c(-1,0,1,2,Inf)).
#' @param bin.labels Character vector of labels for expression bins (default: c("0","1","2","3+")).
#' @param fill.cols Named vector of colors for cell types (default: NULL, uses Set2 palette).
#' @param theme Plot theme, passed to `plot_theme` (default: "classic").
#' @param leg.pos Legend position (default: "right"). Can be numeric vector c(x,y).
#' @param x.angle Angle for x-axis labels (default: 0).
#' @param font.size Base font size (default: 10).
#' @param x.lab X-axis label (default: "Number of genes expressed").
#' @param ann.counts Logical, whether to annotate bar counts (default: FALSE).
#' @param show.sig Logical, whether to show significance stars (default: TRUE).
#' @param pval Logical, whether to annotate chi-square p-values (default: TRUE).
#' @param facet.ncol Number of columns for faceted plots (default: 2).
#' @param label.size Size of annotation text for counts (default: 3).
#' @param sig.size Size of significance stars (default: 4).
#' @param facet.label.size Size of facet labels (default: 10).
#' @param leg.size Legend text size (default: 10).
#' @param leg.ttl.size Legend title size (default: 10).
#' @param panel.fill Fill color for panel background (default: "white").
#' @param facet.bg Logical, whether to display facet background (default: TRUE).
#' @param verbose Logical, whether to print messages (default: TRUE).
#' @param ... Additional arguments passed to `plot_theme`.
#'
#' @return A list containing:
#' \item{plot}{ggplot object.}
#' \item{binned.data}{Data frame of counts and proportions per cell type/bin.}
#' \item{chi.test}{Chi-square test object for binned counts.}
#' \item{residuals}{Standardized residuals from chi-square test.}
#'
#' @details
#' The function calculates the number of genes expressed per cell, bins them according
#' to `breaks` and `bin.labels`, and computes proportions. In faceted mode, it also
#' calculates mean ± SEM for each cell type and overlays it on the bars.
#'
#' @examples
#' # Standard dodged bar plot
#' res_standard <- gsetplot(
#'   obj = seurat_obj,
#'   geneset = c("GeneA","GeneB","GeneC"),
#'   group.by = "celltype",
#'   mode = "standard"
#' )
#' res_standard$plot
#'
#' # Stacked bar plot
#' res_stacked <- gsetplot(
#'   obj = seurat_obj,
#'   geneset = c("GeneA","GeneB","GeneC"),
#'   group.by = "celltype",
#'   mode = "stacked"
#' )
#' res_stacked$plot
#'
#' # Faceted plot with mean ± SEM overlay
#' res_faceted <- gsetplot(
#'   obj = seurat_obj,
#'   geneset = c("GeneA","GeneB","GeneC"),
#'   group.by = "celltype",
#'   stage = "stage",
#'   mode = "faceted",
#'   facet.ncol = 3
#' )
#' res_faceted$plot
#'
#' @export
#' 
gsetplot <- function(
    obj,
    geneset,
    group.by = "ann2",
    stage = "stage",
    mode = c("standard", "stacked", "faceted"),
    assay = "RNA",
    slot = "data",
    expr.thresh = 0,
    breaks = c(-1, 0, 1, 2, Inf),
    bin.labels = c("0", "1", "2", "3+"),
    fill.cols = NULL,
    theme = "classic",
    leg.pos = "right",
    x.angle = 0,
    font.size = 10,
    x.lab = "Number of genes expressed",
    ann.counts = FALSE,
    show.sig = TRUE,
    pval = TRUE,
    facet.ncol = 2,
    label.size = 3,
    sig.size = 4,
    facet.label.size = 10,
    leg.size = 10,
    leg.ttl.size = 10,
    panel.fill = "white",
    facet.bg = TRUE,
    verbose = TRUE,
    ...
) {
  
  require(Seurat)
  require(dplyr)
  require(ggplot2)
  require(RColorBrewer)
  
  mode <- match.arg(mode)
  
  # 1. Extract genes 
  genes.available <- intersect(geneset, rownames(obj))
  if(length(genes.available)==0) stop("No genes from geneset found in object.")
  
  expr.mat <- GetAssayData(obj, assay=assay, leyer=slot)[genes.available,,drop=FALSE]
  if(!is.matrix(expr.mat)) expr.mat <- as.matrix(expr.mat)
  
  # 2. Count expressed genes per cell 
  n.expr <- colSums(expr.mat > expr.thresh)
  
  meta <- obj@meta.data
  if(!group.by %in% colnames(meta)) stop("group.by not found in meta.data")
  if(mode=="faceted" && !stage %in% colnames(meta)) stop("stage not found in meta.data")
  
  df <- data.frame(
    cell = colnames(expr.mat),
    cell.type = as.character(meta[[group.by]]),
    n.expr = n.expr,
    stringsAsFactors = FALSE
  )
  if(mode=="faceted") df$stage <- factor(meta[[stage]], levels = sort(unique(meta[[stage]])))
  
  # 3. Bin counts
  if(length(breaks)-1 != length(bin.labels)) stop("bin.labels must match breaks-1")
  df$expr.bin <- cut(df$n.expr, breaks=breaks, labels=bin.labels, right=TRUE, include.lowest=TRUE)
  
  # 4. Compute proportions + SEM per cell type (binomial SEM) 
  plot.df <- df %>%
    group_by(cell.type, expr.bin) %>%
    summarise(count=n(), .groups="drop") %>%
    group_by(cell.type) %>%
    mutate(prop = count / sum(count) * 100,
           sem = sqrt(prop*(100-prop)/sum(count))) %>%
    ungroup()
  
  # 5. Chi-square significance (per bin)
  residuals.df <- NULL
  
  # Compute per-bin 2xN chi-square tests
  bin_pvals <- sapply(levels(df$expr.bin), function(bin) {
    tab_bin <- table(df$cell.type, df$expr.bin == bin) # in-bin vs not-in-bin
    if(all(dim(tab_bin) == c(length(unique(df$cell.type)), 2))) {
      test <- tryCatch(chisq.test(tab_bin), error=function(e) NULL)
      if(!is.null(test)) return(test$p.value)
    }
    return(NA)
  })
  names(bin_pvals) <- levels(df$expr.bin)
  
  # Overall chi-square residuals (optional, for stars)
  tab <- table(df$cell.type, df$expr.bin)
  chi.res <- tryCatch(chisq.test(tab), error=function(e){warning(e); return(NULL)})
  if(!is.null(chi.res)){
    stdres <- chi.res$stdres
    residuals.df <- as.data.frame(as.table(stdres))
    colnames(residuals.df) <- c("cell.type","expr.bin","stdres")
    residuals.df$sig <- dplyr::case_when(
      abs(residuals.df$stdres)>=3.29 ~ "***",
      abs(residuals.df$stdres)>=2.58 ~ "**",
      abs(residuals.df$stdres)>=1.96 ~ "*",
      TRUE ~ ""
    )
    plot.df <- left_join(plot.df, residuals.df[,c("cell.type","expr.bin","sig")], by=c("cell.type","expr.bin"))
  } else {
    plot.df$sig <- ""
  }
  
  # Add per-bin p-value
  plot.df$chi.p.value <- bin_pvals[as.character(plot.df$expr.bin)]
  
  
  # 6. Default colors 
  if(is.null(fill.cols)){
    if(mode=="stacked") {
      levels_vec <- unique(plot.df$expr.bin)
      pal <- RColorBrewer::brewer.pal(max(3,min(8,length(levels_vec))), "Set3")
      fill.cols <- setNames(colorRampPalette(pal)(length(levels_vec)), levels_vec)
    } else {
      levels_vec <- unique(plot.df$cell.type)
      pal <- RColorBrewer::brewer.pal(max(3,min(8,length(levels_vec))), "Set2")
      fill.cols <- setNames(colorRampPalette(pal)(length(levels_vec)), levels_vec)
    }
  }
  
  # 7. Plot modes 
  if(mode=="standard"){
    p <- ggplot(plot.df, aes(x=expr.bin, y=prop, fill=cell.type)) +
      geom_col(position=position_dodge(width=0.6), color="black", linewidth=0.2, width=0.5) +
      geom_errorbar(aes(ymin=prop-sem, ymax=prop+sem), position=position_dodge(width=0.6), width=0.2, linewidth=0.25) +
      scale_fill_manual(values=fill.cols) +
      labs(x=x.lab, y="Proportion (%)", fill="") +
      plot_theme(theme=theme, font.size=font.size, x.angle=x.angle, leg.pos=leg.pos, ...) +
      theme(strip.text=element_text(size=facet.label.size),
            legend.text=element_text(size=leg.size),
            legend.title=element_text(size=leg.ttl.size),
            panel.background=element_rect(fill=panel.fill, colour=NA))
    if(ann.counts) p <- p + geom_text(aes(label=count), position=position_dodge(width=0.9), vjust=-0.4, size=label.size)
    if(show.sig) p <- p + geom_text(aes(label=sig), position=position_dodge(width=0.9), vjust=-1, size=sig.size)
  } else if(mode=="stacked"){
    plot.df <- plot.df %>% group_by(expr.bin) %>% mutate(prop = count/sum(count)*100) %>% ungroup()
    p <- ggplot(plot.df, aes(x=cell.type, y=prop, fill=expr.bin)) +
      geom_col(color="black", linewidth=0.2, width=0.5) +
      scale_fill_manual(values=fill.cols, name=x.lab) +
      labs(x="Cell type", y="Proportion (%)") +
      plot_theme(theme=theme, font.size=font.size, x.angle=x.angle, leg.pos=leg.pos, ...) +
      theme(strip.text=element_text(size=facet.label.size),
            legend.text=element_text(size=leg.size),
            legend.title=element_text(size=leg.ttl.size),
            panel.background=element_rect(fill=panel.fill, colour=NA))
  } else if(mode=="faceted"){
    # Binned prop for facets 
    plot.df <- df %>%
      group_by(.data[[stage]], cell.type, expr.bin) %>%
      summarise(count=n(), .groups="drop") %>%
      group_by(.data[[stage]], cell.type) %>%
      mutate(prop = count / sum(count) * 100) %>%
      ungroup()
    
    # Mean ± SEM per cell type 
    mean.df <- df %>%
      group_by(.data[[stage]], cell.type) %>%
      summarise(
        mean.expr = mean(n.expr),
        sem.expr = sd(n.expr)/sqrt(n()),
        .groups="drop"
      )
    
    p <- ggplot(plot.df, aes(x=expr.bin, y=prop, fill=cell.type)) +
      geom_col(position=position_dodge(width=0.7), color="black", linewidth=0.2, width=0.5) +
      facet_wrap(vars(.data[[stage]]), ncol=facet.ncol) +
      scale_fill_manual(values=fill.cols) +
      labs(x=x.lab, y="Proportion (%)", fill="") +
      plot_theme(theme=theme, font.size=font.size, leg.pos=leg.pos, facet.bg=facet.bg, ...) +
      theme(strip.text=element_text(size=facet.label.size),
            panel.background=element_rect(fill=panel.fill, colour=NA)) +
      # Add mean ± SEM on top of bars 
      # geom_errorbar(
      #   data=mean.df,
      #   aes(x=cell.type, ymin=mean.expr-sem.expr, ymax=mean.expr+sem.expr, color=cell.type),
      #   inherit.aes=FALSE,
      #   position=position_dodge(width=0.7),
      #   width=0.3,
      #   linewidth=0.6
      # ) +
      # geom_point(
      #   data=mean.df,
      #   aes(x=cell.type, y=mean.expr, color=cell.type),
      #   inherit.aes=FALSE,
      #   position=position_dodge(width=0.7),
      #   size=2
      # ) +
      scale_color_manual(values=fill.cols, guide="none")
  }
  
  return(list(plot=p, binned.data=plot.df, chi.test=chi.res, residuals=residuals.df))
}


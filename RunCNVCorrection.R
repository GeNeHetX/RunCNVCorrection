# RunCNVCorrection(): regress arm-level CNV scores out of a Seurat PCA. See README.md.

RunCNVCorrection <- function(object, cnv.pattern = "^cnv_fastcnv_",
                             cnv.exclude = "fraction|cluster|_pc[0-9]+$",
                             reduction = "pca", reduction.name = "pca.cnv",
                             reduction.key = "PCCNV_", verbose = TRUE) {

  # Checks
  if (!inherits(object, "Seurat")) stop("object must be a Seurat object")
  if (!reduction %in% names(object@reductions))
    stop("reduction '", reduction, "' not found: run RunPCA() first")
  cols <- grep(cnv.pattern, colnames(object@meta.data), value = TRUE)
  if (!is.null(cnv.exclude)) cols <- cols[!grepl(cnv.exclude, cols)]   # keep arm scores only
  if (length(cols) < 2) stop("no CNV arm scores found in the metadata: run a CNV method first")

  A <- Seurat::Embeddings(object, reduction)
  X <- as.matrix(object@meta.data[rownames(A), cols])
  storage.mode(X) <- "double"
  altered <- rowSums(!is.na(X)) > 0
  if (!any(altered)) stop("CNV scores are all missing")
  if (verbose) message(sprintf("RunCNVCorrection: %d arms, %d/%d cells with a CNV profile",
                               length(cols), sum(altered), length(altered)))

  # Cells without CNV profile = diploid (score 0); constant arms removed
  X[is.na(X)] <- 0
  X <- X[, apply(X, 2, var) > 0, drop = FALSE]

  # Residuals of each PC on all arm scores (least squares with intercept)
  R <- qr.resid(qr(cbind(1, X)), A)
  dimnames(R) <- list(rownames(A), paste0(reduction.key, seq_len(ncol(R))))
  object[[reduction.name]] <- Seurat::CreateDimReducObject(embeddings = R, key = reduction.key,
                                                           assay = Seurat::DefaultAssay(object))
  object
}

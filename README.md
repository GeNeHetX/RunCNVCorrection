# RunCNVCorrection

Regress arm-level copy number variation (CNV) scores out of the PCA of a Seurat object, so that malignant cells no longer cluster by patient.

> Pignolet C, Becht E, Nicolle R. *Regressing out copy number variation signal removes patient clustering of single-cell malignant transcriptomes.*

In multi-patient tumor scRNA-seq, each tumor carries its own copy number alterations, which shift the expression of whole chromosomal arms. Malignant cells therefore group by patient rather than by cell state. `RunCNVCorrection()` regresses each principal component on the per-cell CNV arm scores and stores the residuals as a new reduction. The correction acts mainly on malignant cells; the microenvironment, which carries almost no CNV, is barely moved.

## Requirements

- R and Seurat (≥ 4)
- A PCA (`RunPCA()`)
- Per-cell CNV scores per chromosomal arm in the metadata, e.g. from fastCNV (`cnv_fastcnv_1.p`, `cnv_fastcnv_1.q`, ..., `cnv_fastcnv_X.q`)

## Usage

```r
source("RunCNVCorrection.R")

RunCNVCorrection(
  object,
  cnv.pattern    = "^cnv_fastcnv_",
  cnv.exclude    = "fraction|cluster|_pc[0-9]+$",
  reduction      = "pca",
  reduction.name = "pca.cnv",
  reduction.key  = "PCCNV_",
  verbose        = TRUE
)
```

### Arguments

- **`object`**: a Seurat object with a PCA and per-cell CNV arm scores in `object@meta.data`.
- **`cnv.pattern`**: regular expression that selects the CNV columns of the metadata. The default matches fastCNV columns (`cnv_fastcnv_1.p`, `cnv_fastcnv_1.q`, ...). With another CNV method, give the prefix of its arm columns, e.g. `cnv.pattern = "^infercnv_arm_"`.
- **`cnv.exclude`**: regular expression of columns that match `cnv.pattern` but are not arm scores, such as the CNV burden (`fraction`), CNV clusters (`cluster`) or CNV principal components (`_pc1`, `_pc2`, ...). These columns are not used. Set to `NULL` to keep all matched columns.
- **`reduction`**: name of the PCA to correct. Use another name if your PCA is stored under a different reduction.
- **`reduction.name`**: name of the new reduction that holds the corrected PCA.
- **`reduction.key`**: prefix of the columns of the new reduction (`PCCNV_1`, `PCCNV_2`, ...), as `PC_` for the PCA.
- **`verbose`**: if `TRUE`, prints the number of CNV arms used and the number of cells with a CNV profile.

### Value

The same Seurat object with a new reduction (`"pca.cnv"` by default) that has the same number of components as the input PCA. Use it with `reduction = "pca.cnv"` in `FindNeighbors()` and `RunUMAP()`. Expression values and the original PCA are not modified.

## Example

Full Seurat workflow on tumor scRNA-seq from several patients. CNV arm scores are assumed to have been added to the metadata by fastCNV, run per patient (see the fastCNV documentation).

```r
library(Seurat)
source("RunCNVCorrection.R")

# 1. Seurat object and quality control
obj <- CreateSeuratObject(counts, meta.data = metadata, min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)

# 2. Normalization, variable genes, scaling, PCA
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)

# 3. CNV inference with fastCNV (per patient), which adds the arm scores to the metadata
#    ... see fastCNV documentation ...
head(grep("^cnv_fastcnv_", colnames(obj@meta.data), value = TRUE))
#> "cnv_fastcnv_1.p" "cnv_fastcnv_1.q" "cnv_fastcnv_2.p" ...

# 4. Uncorrected embedding (for comparison)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:20, reduction.name = "umap.raw")

# 5. CNV correction
obj <- RunCNVCorrection(obj)
#> RunCNVCorrection: 46 arms, 86013/86013 cells with a CNV profile

# 6. Clustering and UMAP on the corrected PCA
obj <- FindNeighbors(obj, reduction = "pca.cnv", dims = 1:20)
obj <- FindClusters(obj, resolution = 0.3)
obj <- RunUMAP(obj, reduction = "pca.cnv", dims = 1:20, reduction.name = "umap.cnv")

# 7. Before / after: malignant cells group by patient before correction, mix after
DimPlot(obj, reduction = "umap.raw", group.by = "patient_id") +
  DimPlot(obj, reduction = "umap.cnv", group.by = "patient_id")
DimPlot(obj, reduction = "umap.raw", group.by = "cell_type") +
  DimPlot(obj, reduction = "umap.cnv", group.by = "cell_type")

# 8. Markers on the corrected clusters (expression values are unchanged)
markers <- FindAllMarkers(obj, only.pos = TRUE)
```

## Checks

The function stops if the object is not a Seurat object, if the PCA is missing, or if no CNV arm scores are found. Cells without a CNV profile get scores of 0: they are treated as diploid and barely moved. Arms that do not vary across cells are removed.

## Method

Each principal component is regressed on all arm scores by ordinary least squares with an intercept, and the residuals form the corrected PCA. Scores are used per arm rather than per gene to keep a small number of covariates.

## Note

States that follow copy number alterations, or subclones that differ by them, may no longer form separate clusters in the corrected space. Keep the uncorrected PCA for clonal or CNA-driven analyses.

## License

[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) © 2026 Camille Pignolet. See [LICENSE](LICENSE).

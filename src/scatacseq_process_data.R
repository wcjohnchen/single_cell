# =================================================================== load packages
library(Signac)
library(Seurat)
library(SeuratDisk)

library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(biovizBase)

library(SingleCellExperiment)
library(scDblFinder)

library(singleCellTK)

# =================================================================== read matrix
counts <- ReadMtx(
  mtx = "./GSM5717477_MCF7_DMSO_filtered_peak_bc_matrix.mtx.gz", 
  features = "./GSM5717477_MCF7_DMSO_filtered_peak_bc_peaks2_concat.tsv.gz",
  cells = "./GSM5717477_MCF7_DMSO_filtered_peak_bc_barcodes.tsv.gz",
  feature.column = 1
)
dim(counts)

# =================================================================== perform doublet detection
sce <- SingleCellExperiment(assays = list(counts = counts))
sce3 <- scDblFinder(sce, aggregateFeatures=TRUE)
sce3$scDblFinder.score 
sce3$scDblFinder.class

table(sce3$scDblFinder.class)
sce4 <- subsetSCECols(sce3, colData = "scDblFinder.class == 'singlet'")

mat <- assay(sce4, "counts")
dim(mat)

# =================================================================== create chromatin assay
chrom_assay <- CreateChromatinAssay(
  counts = mat,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = './GSM5717477_MCF7_DMSO_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

# =================================================================== create seurat object
data <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks',
)

# =================================================================== annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
annotations
Annotation(data) <- annotations

# =================================================================== metrics, nucleosome signals, calculate ratio of fragments from monocucleosomes versus fragments from nucleosome-free
data <- NucleosomeSignal(object = data)

# =================================================================== calculate TSS enrichment
data <- TSSEnrichment(object = data, fast = FALSE)

VlnPlot(
  object = data,
  features = c('nucleosome_signal', 'TSS.enrichment'),
  pt.size = 0.1,
  ncol = 2
)

# =================================================================== filtering based on quantile percentages
hig_ns <- quantile(data[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)
low_ts <- quantile(data[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)
hig_ns
low_ts

data <- subset(
  x = data,
  subset = nucleosome_signal < hig_ns &
    TSS.enrichment > low_ts
)
data

# =================================================================== normalization and dimensional reduction
data <- RunTFIDF(data)

data <- FindTopFeatures(data, min.cutoff = 'q0')
data

data <- RunSVD(data)
DepthCor(data)

data <- RunUMAP(object = data, reduction = 'lsi', dims = 3:30)
data <- FindNeighbors(object = data, reduction = 'lsi', dims = 3:30)
data <- FindClusters(object = data, verbose = FALSE, algorithm = 3)
DimPlot(object = data, label = TRUE) + NoLegend()

# =================================================================== calculate gene activity
gene.activities <- GeneActivity(data)


data[['RNA']] <- CreateAssayObject(counts = gene.activities)

data <- NormalizeData(
  object = data,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(data$nCount_RNA)
)
data[['RNA']]

rownames(data[['RNA']])
rownames(data[[]])

DefaultAssay(data) <- 'RNA'

SaveH5Seurat(data, filename = "sc_atac_seq.h5Seurat")
Convert("sc_atac_seq.h5Seurat", dest = "h5ad")



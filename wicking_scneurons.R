# Code following 30/09/25 meeting   
# What transcriptome did this align to? In which version of Parse Bio Trailmaker?

## USER SELECTIONS ##
n.pcs <- 20 # Okay to be high
x.res <- 0.4 # I set this slightly lower to try see high level signatures
setwd("~/scneurons25")
## DONE ##

lib.dir <- file.path("~/scneurons25/software/R_LIBS")
data.dir <- "~/scneurons25/data"
out.dir <- "~/scneurons25/output"
dir.create(data.dir)
dir.create(out.dir)
# Set custom library path for installations
dir.create(lib.dir, recursive = TRUE)
.libPaths(c(lib.dir, .libPaths()))

install.packages("BiocManager", lib = lib.dir)
install.packages("devtools", lib = lib.dir)
install.packages("Seurat", lib = lib.dir) # Incl. Patchwork
install.packages("dplyr")
install.packages("HGNChelper")
install.packages("openxlsx")
BiocManager::install("DESeq2", lib = lib.dir)
# install.packages("ScPred", lib = lib.dir) # Needs fork to work with v5
devtools::install_github(repo="powellgenomicslab/scPred",
    ref="9f407b7436f40d44224a5976a94cc6815c6e837f", lib = lib.dir)
system("git clone https://github.com/IanevskiAleksandr/sc-type")
sctype.dir <- file.path(getwd(), "sc-type")

library(Seurat)
library(scPred)

data.dir <- file.path(getwd(), "data")
out.dir <- file.path(getwd(), "output")

# Data dir must have *.mtx, cell_metadata.csv, and all_genes.csv
main.data <- ReadParseBio(data.dir)

# TODO(TG): Check ParseBio filter defaults (probably just use the filtered mtx)
main.data <- CreateSeuratObject(counts = main.data, project = "scneurons25", min.cells = 3, min.features = 200)

print("Created Seurat object")
print(main.data)
print(head(main.data), 5)
print(main.data[1:3, 1:10])

print(sprintf("Number of cells: ", length(colnames(main.data))))
print(sprintf("Number of features: ", length(rownames(main.data))))


# Add QC metrics
main.data[["percent.mt"]] <- PercentageFeatureSet(main.data, pattern = "^MT-")
main.data[["percent.ribo"]] <- PercentageFeatureSet(main.data, pattern = "^RPS|^RPL")
main.data[["complexity"]] <- main.data$nFeature_RNA / main.data$nCount_RNA

pdf(file.path(out.dir, "qc_plots.pdf"), width=11.7, height=8.3)
VlnPlot(main.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "complexity"), ncol = 5)

plot1 <- FeatureScatter(main.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(main.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot3 <- FeatureScatter(main.data, feature1 = "nCount_RNA", feature2 = "percent.ribo")
print(plot1 + plot2 + plot3)
dev.off()

# Uses the defaults
main.data <- NormalizeData(main.data)
# https://doi.org/10.1016/j.cell.2019.05.031
main.data <- FindVariableFeatures(main.data, selection.method = "vst", nfeatures = 2000)


main.data <- ScaleData(main.data) # Scales using variable features by default
main.data <- RunPCA(main.data, features = VariableFeatures(object = main.data))
print(main.data[["pca"]], dims = 1:5, nfeatures = 10)

pdf(file.path(out.dir, "lineardim_plots.pdf"), width=8.3, height=11.7)
VizDimLoadings(main.data, dims = 1:5, nfeatures = 10, reduction = "pca", ncol=5)
DimPlot(main.data, reduction = "pca") + NoLegend()
DimHeatmap(main.data, dims = 1:5, cells = 250, balanced = TRUE)
dev.off()

pdf(file.path(out.dir, "elbow_plot.pdf"), width=14, height=7)
ElbowPlot(main.data, ndims = 30) # User decides how many PCs
dev.off()

print(sprintf("Select number of PCs on ElbowPlot. Script using: ", n.pcs))

main.data <- FindNeighbors(main.data, dims = 1:n.pcs)
main.data <- FindClusters(main.data, resolution = x.res) # User decides resolution

print("Number of cells per cluster: ")
print(table(Idents(main.data)))

main.data <- RunUMAP(main.data, dims = 1:n.pcs)

pdf(file.path(out.dir, "nonlineardim_plots.pdf"), width=11.7, height=8.3)
DimPlot(main.data, reduction = "umap", label = TRUE) + NoLegend()
dev.off()

saveRDS(main.data, file = file.path(out.dir, "scneurons25_main.RDS"))
q()

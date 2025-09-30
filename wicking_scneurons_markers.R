setwd("~/scneurons25")

lib.dir <- file.path("~/scneurons25/software/R_LIBS")
data.dir <- "~/scneurons25/data"
out.dir <- "~/scneurons25/output"
dir.create(data.dir)
dir.create(out.dir)
# Set custom library path
dir.create(lib.dir, recursive = TRUE)
.libPaths(c(lib.dir, .libPaths()))

library(BiocManager)
library(devtools)
library(Seurat)
library(dplyr)
library(HGNChelper)
library(openxlsx)

# See https://github.com/IanevskiAleksandr/sc-type
sctype.dir <- file.path(getwd(), "sc-type")

main.data <- readRDS(file = file.path(out.dir, "scneurons25_main.RDS"))

source(file.path(sctype.dir), "gene_sets_prepare.R"); source(file.path(sctype.dir), "sctype_score_.R")
gene.sets <- gene_sets_prepare(file.path(file.path(sctype.dir), "ScTypeDB_full.xlsx", "Brain"))

scale.data <- as.matrix(main.data[["RNA"]]$scale.data)

# TODO(TG): Curious what format is here
sctype.pred <- sctype_score(scRNAseqData = scale.data, scaled = TRUE, gs = gene.sets$gs_positive, gs2 = gene.sets$gs_negative)

# We could avoid various spaghetti code and use the provided wrapper for less control
sctype.results <- do.call("rbind", lapply(unique(main.data@meta.data$seurat_clusters), function(cl){
    sctype.cl.pred = sort(rowSums(sctype.pred[ ,rownames(main.data@meta.data[main.data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(sctype.cl.pred), scores = sctype.cl.pred, ncells = sum(main.data@meta.data$seurat_clusters==cl)), 10)
}))
sctype.scores <- do.call(rbind, lapply(split(sctype.results, sctype.results$cluster), function(df) {
    df[which.max(df$scores), ]
}))
# Set low ScType score clusters to "Unknown"
sctype.scores$type[as.numeric(as.character(sctype.scores$scores)) < sctype.scores$ncells/4] <- "Unknown"
print(sctype.scores[,1:3])

# Add ScType annotation to main meta data
main.data@meta.data$sctype.pred = ""
for(j in unique(sctype.scores$cluster)){
    sctype.tmp = sctype.scores[sctype.scores$cluster==j,]
    main.data@meta.data$sctype.pred[main.data@meta.data$seurat_clusters == j] = as.character(sctype.tmp$type[1])
}
rm(sctype.tmp)
# Make sctype.pred safe for use in file paths
main.data@meta.data$sctype.pred.format <- gsub("[^A-Za-z0-9_-]", "_", as.character(main.data@meta.data$sctype.pred))

pdf(file.path(out.dir, "umap_sctype_classification.pdf"), width=11.7, height=8.3)
DimPlot(main.data, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype.pred')        
dev.off()

print("Check numbers per annotation: ")
print("Seurat clusters")
print(table(main.data$seurat_clusters))
print("ScType clusters")
print(table(main.data$sctype.pred))

# We run find all markers with the Seurat clusters, though they now have ScType annotations
# Can check which will be used with: Idents(main.data)
# Slightly lower than defaults to keep sensitivity
main.markers <- FindAllMarkers(main.data, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(main.markers, file = file.path(out.dir, "scneurons25_seurat_allmarkers.csv"))

sig.markers <- main.markers[main.markers$p_val_adj < 0.05, ] # Save just significant
write.csv(main.markers, file = file.path(out.dir, "scneurons25_sig_seurat_allmarkers.csv"))



clusters <- unique(sig.markers$cluster)
for (clust in clusters) {
    clust.markers <- sig.markers[sig.markers$cluster == clust, ]
    top.genes <- head(clust.markers[order(clust.markers$p_val_adj), "gene"], 20)
    # Dot plot for each cluster
    pdf(file = file.path(out.dir, paste0("cluster_", clust, "_sig_allmarkers_20dotplot.pdf")), width = 11.7, height = 8.3)
    print(DotPlot(main.data, features = top.genes, group.by = "seurat_clusters") +
        ggtitle(paste("Cluster", clust, "Top 20 Markers (Padj < 0.05)")))
    dev.off()
}



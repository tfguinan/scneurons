Code following 30/09/25 meeting
What transcriptome did the data align to? In which version of Parse Bio Trailmaker?

Directory structure:

~/scneurons25
  ./data
  ./software
  ./output

QC filters for >3 cells and >200 features (use the Trailmaker mtx filter for their empty and doublet algorithm)
Uses 20 PCs and cluster resolution of 0.4, tests genes with minimum percent detection in one type >0.1 and log fold change >0.1

Signficance set at P < 0.05


Output Files:
scneurons25_main.RDS: Processed Seurat object for downstream analysis.
qc_plots.pdf: Quality control plots (violin and scatter plots).
lineardim_plots.pdf: PCA visualizations and heatmaps.
elbow_plot.pdf: Elbow plot for selecting principal components.
nonlineardim_plots.pdf: UMAP plots for cluster visualization.
umap_sctype_classification.pdf: UMAP with ScType annotations.
scneurons25_seurat_allmarkers.csv: All identified markers for clusters.
scneurons25_sig_seurat_allmarkers.csv: Significant markers for clusters.
cluster_X_sig_allmarkers_20dotplot.pdf: Dot plots of top 20 markers per cluster (one file per cluster).

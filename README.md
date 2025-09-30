## Code Following 30/09/25 Meeting

**Questions:**
- What transcriptome did the data align to?
- Which version of Parse Bio Trailmaker was used?

---

### Directory Structure

```
~/scneurons25/
├── data/
├── software/
└── output/
```

---

### QC Filters

- **Cells:** > 3 cells
- **Features:** > 200 features
- **Empty/Doublet Filtering:** Use Trailmaker mtx filter (empty and doublet algorithm)
- **Principal Components:** 20 PCs
- **Cluster Resolution:** 0.4
- **Gene Testing:** Minimum percent detection in one type > 0.1, log fold change > 0.1
- **Significance Threshold:** P < 0.05

---

### Output Files

| File Name                                 | Description                                         |
|--------------------------------------------|-----------------------------------------------------|
| `scneurons25_main.RDS`                     | Processed Seurat object for downstream analysis     |
| `qc_plots.pdf`                             | Quality control plots (violin and scatter plots)    |
| `lineardim_plots.pdf`                      | PCA visualizations and heatmaps                    |
| `elbow_plot.pdf`                           | Elbow plot for selecting principal components       |
| `nonlineardim_plots.pdf`                   | UMAP plots for cluster visualization                |
| `umap_sctype_classification.pdf`           | UMAP with ScType annotations                        |
| `scneurons25_seurat_allmarkers.csv`        | All identified markers for clusters                 |
| `scneurons25_sig_seurat_allmarkers.csv`    | Significant markers for clusters                    |
| `cluster_X_sig_allmarkers_20dotplot.pdf`   | Dot plots of top 20 markers per cluster (one per cluster) |

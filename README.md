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

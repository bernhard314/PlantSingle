# PlantSingle

Initial notes and toy example of Seurat based plant root cell analysis based on public data.
Download raw files from quantification and Meta data from

https://www.ebi.ac.uk/gxa/sc/experiments/E-GEOD-121619/downloads

unzip and rename files to matrix.mtx, genes.tsv and barcodes.tsv

Steps:
1. Load and normalize 
2. Add Meta
3. Split conditions and reintegrate to cluster based on Cell types 
4. Identify clusters based on annotation found in  www.ncbi.nlm.nih.gov/pubmed/31091459 (redo with original publication?)
5. Calculate cell type based differential expressions  


Raw data to www.ncbi.nlm.nih.gov/pubmed/31091459 is available on https://www.ebi.ac.uk/gxa/sc/experiments/E-CURD-4/downloads .


Work in progress. 


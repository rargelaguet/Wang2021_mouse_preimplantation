# TO-DO: FILL DIMENSIONALITY REDUCTION SLOTS

library(SingleCellExperiment)
library(data.table)
library(purrr)

sce <- readRDS("/Users/argelagr/data/Wang2021/original/20190909.RNA.mtx.sce.rds")

logcounts(sce) <- NULL

cols <- c("batch", "stage", "embryo_id", "sex", "cnv_info", 
          "WCG_sites", "GCH_sites", "gch_mean_rate", "wcg_mean_rate", 
          "chrM_gch_rate", "chrM_wcg_rate", "lambda_gch_rate", "lambda_wcg_rate", 
          "gene_num", "cluster", "qc.rna", "qc.dna", "qc.acc.extra", "cdr", 
          "DC1", "DC2", "DC3", "DC.order", "total_counts"
)
metadata.dt <- colData(sce)[,cols] %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","cell") %>%
  setnames(c("qc.rna","qc.dna","qc.acc.extra"),c("pass_rnaQC","pass_metQC","pass_accQC")) %>%
  setnames("cdr","num_genes")


metadata.dt %>%
  .[stage=="4cell",stage:="early_4cell"] %>% .[stage=="L4cell",stage:="late_4cell"] %>% 
  setnames("stage","celltype")

table(metadata.dt$celltype)

 # Rename genes 
gene_metadata.dt <- fread("/Users/argelagr/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt")
ens_ids <- intersect(gene_metadata.dt$ens_id, rownames(sce))

sce <- sce[ens_ids,]
gene_metadata.dt <- gene_metadata.dt %>%
  .[ens_id%in%rownames(sce)] %>% 
  setkey(ens_id) %>% .[rownames(sce)]
rownames(sce) <- gene_metadata.dt$symbol

# Save
fwrite(metadata.dt, "/Users/argelagr/data/Wang2021/cell_metadata.txt.gz", sep="\t", quote=F, na="NA")
saveRDS(sce, "/Users/argelagr/data/Wang2021/processed/rna/SingleCellExperiment.rds")

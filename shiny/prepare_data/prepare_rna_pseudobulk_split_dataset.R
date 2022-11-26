source(here::here("settings.R"))
source(here::here("utils.R"))
source(here::here("rna/pseudobulk/utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"shiny/cell_metadata.txt.gz")
io$rna.sce <- file.path(io$basedir,"shiny/rna_expression/SingleCellExperiment.rds")
io$genes <- file.path(io$basedir,"shiny/rna_expression/genes.txt")
io$outdir <- file.path(io$basedir,"shiny/rna_expression/pseudobulk"); dir.create(io$outdir, showWarnings = F)

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[,dataset_stage_lineage:=sprintf("%s-%s-%s",dataset,stage,lineage)]

dim(cell_metadata.dt)

######################################################
## Calculate pseudovulk stats and do some filtering ##
######################################################

pseudobulk_stats.dt <- cell_metadata.dt[,.N,by=c("dataset","stage","lineage","dataset_stage_lineage")]

# Filter each instance by minimum number of cells
pseudobulk_stats.dt <- pseudobulk_stats.dt[N>=10]

# Select celltypes that are measured in at least 5 WT samples
# celltypes.to.use <- pseudobulk_stats.dt[class=="WT",.N,by=c("celltype")] %>% .[N>=5,celltype] 
# pseudobulk_stats.dt <- pseudobulk_stats.dt[celltype%in%celltypes.to.use]

# For each class and celltype combination, require at least 3 samples
# tmp <- pseudobulk_stats.dt[,.N,by=c("celltype","class")] %>% .[N>=3] %>% .[,N:=NULL]
# pseudobulk_stats.dt <- pseudobulk_stats.dt %>% merge(tmp,by=c("celltype","class"))

# Print stats
# print(pseudobulk_stats.dt)

# Update metadata
cell_metadata.dt <- cell_metadata.dt[dataset_stage_lineage%in%pseudobulk_stats.dt$dataset_stage_lineage]

dim(cell_metadata.dt)

##############################
## Load RNA expression data ##
##############################

sce <- load_SingleCellExperiment(io$rna.sce, cells=cell_metadata.dt$id_rna)
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("id_rna") %>% DataFrame

sce$dataset_stage_lineage <- sprintf("%s-%s-%s",sce$dataset,sce$stage,sce$lineage)

################
## Pseudobulk ##
################

sce_pseudobulk <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = "dataset_stage_lineage",
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_pseudobulk) <- "counts"

##################
## Subset genes ##
##################

genes <- fread(io$genes, header = F)[[1]]
sce_pseudobulk <- sce_pseudobulk[genes,]

###################
## Normalisation ##
###################

logcounts(sce_pseudobulk) <- log2(1e6*(sweep(counts(sce_pseudobulk),2,colSums(counts(sce_pseudobulk)),"/"))+1)

# Remove counts assay
# assays(sce_pseudobulk)["counts"] <- NULL

##########
## Save ##
##########

# Save sample metadata
to.save <- pseudobulk_stats.dt %>% copy %>% setnames("dataset_stage_lineage","id")
fwrite(to.save, file.path(io$outdir,"sample_metadata.txt.gz"), na="NA", quote=F, sep="\t")
       
# Save expression matrix
saveRDS(round(logcounts(sce_pseudobulk),2), file.path(io$outdir,"rna_expr.rds"))

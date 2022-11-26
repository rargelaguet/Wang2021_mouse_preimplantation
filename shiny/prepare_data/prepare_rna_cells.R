library(HDF5Array)

source(here::here("settings.R"))
source(here::here("utils.R"))

################
## Define I/O ##
################

# I/O
io$metadata.cells <- file.path(io$basedir,"shiny/cell_metadata.txt.gz")
io$outdir <- file.path(io$basedir,"shiny/rna_expression"); dir.create(io$outdir, showWarnings = F)

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata.cells)

#########################
## Load RNA expression ##
#########################

sce.cells <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = cell_metadata.dt$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = TRUE
)

##################
## Rename genes ##
##################

# gene_metadata.dt <- fread(io$gene.metadata)
# ens_ids <- intersect(gene_metadata.dt$ens_id, rownames(sce.cells))
# 
# sce.cells <- sce.cells[ens_ids,]
# gene_metadata.dt <- gene_metadata.dt %>%
#   .[ens_id%in%rownames(sce.cells)] %>% 
#   setkey(ens_id) %>% .[rownames(sce.cells)]
# rownames(sce.cells) <- gene_metadata.dt$symbol

##################
## Filter genes ##
##################

genes <- rownames(sce.cells)[grep("^Rik|Rik$|^Mt-|^Rps-|^Rpl-|^Gm",rownames(sce.cells),invert=T)]
sce.cells <- sce.cells[genes,]

##########
## Save ##
##########

saveRDS(sce.cells, paste0(io$outdir,"/SingleCellExperiment.rds"))

# genes
write.table(genes, paste0(io$outdir,"/genes.txt"), row.names = F, col.names = F, quote=F)

# cells
cells <- colnames(sce.cells)
write.table(cells, paste0(io$outdir,"/rna_expr_cells.txt"), row.names = F, col.names = F, quote=F)
outfile = paste0(io$outdir,"/rna_expr_cells.hdf5")
if(file.exists(outfile)) { file.remove(outfile) }
writeHDF5Array(x = DelayedArray(round(logcounts(sce.cells),2)), file = outfile, name = "rna_expr_logcounts", verbose = TRUE)

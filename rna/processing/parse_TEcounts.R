source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$inputdir <- file.path(io$basedir,"TEcounts")
io$outdir <- file.path(io$basedir,"processed/rna"); dir.create(io$outdir, showWarnings = F)

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata) %>% .[!is.na(alias)]

opts$samples <- unique(sample_metadata.dt$alias)

###############
## Load data ##
###############

stopifnot(file.exists(file.path(io$inputdir,sprintf("%s.cntTable.gz",opts$samples))))
# opts$samples[!file.exists(file.path(io$inputdir,sprintf("%s.cntTable.gz",opts$samples)))]
# opts$samples <- opts$samples[file.exists(file.path(io$inputdir,sprintf("%s.cntTable.gz",opts$samples)))]

gene_expression_list <- list(); te_expression_list <- list()
          
# i <- "SRR7093810"
for (i in opts$samples) {
	tmp <- fread(file.path(io$inputdir,sprintf("%s.cntTable.gz",i))) %>%
	  setnames(c("feature",i))
	
	gene_expression_list[[i]] <- tmp[grepl("ENS",feature)] %>% matrix.please
	te_expression_list[[i]] <- tmp[!grepl("ENS",feature)] %>% matrix.please
	
}

###########
## Merge ##
###########

gene_expression.mtx <- do.call("cbind",gene_expression_list)
te_expression.mtx <- do.call("cbind",te_expression_list)

tmp <- rownames(te_expression.mtx) %>% strsplit(":")

########################
## Create TE metadata ##
########################

te_metadata.dt <- data.table(
  feature = rownames(te_expression.mtx),
  name = tmp %>% map_chr(1),
  family = tmp %>% map_chr(2),
  class = tmp %>% map_chr(3)
)

##########################
## Create gene metadata ##
##########################

gene_metadata.dt <- fread(io$gene_metadata, select=c("symbol","ens_id")) %>% .[!duplicated(symbol)]
genes.to.use <- intersect(rownames(gene_expression.mtx),gene_metadata.dt$ens_id)
gene_expression.mtx <- gene_expression.mtx[genes.to.use,]
gene_metadata.dt <- gene_metadata.dt[ens_id%in%genes.to.use] %>% setkey(ens_id) %>% .[genes.to.use]

###########################################
## Create SingleCellExperiment for genes ##
###########################################

coldata_to_se <- DataFrame(sample_metadata.dt[alias%in%opts$samples] %>% setkey(alias) %>% .[colnames(gene_expression.mtx)] %>% tibble::column_to_rownames("alias"))
stopifnot(colnames(gene_expression.mtx)==rownames(coldata_to_se))

rowdata_to_se <- DataFrame(gene_metadata.dt %>% tibble::column_to_rownames("ens_id"))
stopifnot(rownames(gene_expression.mtx)==rownames(rowdata_to_se))

genes.se <- SingleCellExperiment(
  assays = list(counts=gene_expression.mtx), 
  rowData = rowdata_to_se,
  colData = coldata_to_se
)

#########################################
## Create SingleCellExperiment for TEs ##
#########################################

coldata_to_se <- DataFrame(sample_metadata.dt[alias%in%opts$samples] %>% setkey(alias) %>% .[colnames(te_expression.mtx)] %>% tibble::column_to_rownames("alias"))
stopifnot(colnames(te_expression.mtx)==rownames(coldata_to_se))

rowdata_to_se <- DataFrame(te_metadata.dt %>% tibble::column_to_rownames("feature"))
stopifnot(rownames(te_expression.mtx)==rownames(rowdata_to_se))

te.se <- SingleCellExperiment(
  assays = list(counts=te_expression.mtx), 
  rowData = rowdata_to_se,
  colData = coldata_to_se
)

############
## Filter ##
############

te_metadata.dt <- te_metadata.dt[!class%in%c("Unknown","LTR?","LINE?","SINE?","Satellite","RC?","DNA?")]

te.se <- te.se[te_metadata.dt$feature,]
te_expression.mtx <- te_expression.mtx[te_metadata.dt$feature,]

##########
## Save ##
##########

write.table(gene_expression.mtx, file.path(io$outdir,"gene_counts_split.txt.gz"), sep = ",", quote = F, na="NA", row.names = T, col.names = T)
write.table(te_expression.mtx, file.path(io$outdir,"TE_counts_split.txt.gz"), sep = ",", quote = F, na="NA", row.names = T, col.names = T)
saveRDS(te.se, file.path(io$outdir,"TE_SingleCellExperiment_split.rds"))
saveRDS(genes.se, file.path(io$outdir,"gene_SingleCellExperiment_split.rds"))

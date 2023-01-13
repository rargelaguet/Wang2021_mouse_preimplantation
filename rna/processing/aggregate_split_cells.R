source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

io$metadata <- file.path(io$basedir,"cell_metadata_split.txt")
io$infile <- file.path(io$basedir,"processed/rna/gene_SingleCellExperiment_split.rds")
io$outdir <- file.path(io$basedir,"processed/rna")

##########################
## Load sample metadata ##
##########################

sample_metadata.dt <- fread(io$metadata) %>% .[!is.na(alias)]

sample_metadata_aggregated.dt <- sample_metadata.dt[,c("stage","celltype","sample_before_split")] %>%
  setnames("sample_before_split","sample") %>% .[,alias:=sample] %>%
  unique

###############
## Load data ##
###############

sce <- readRDS(io$infile)[,sample_metadata.dt$alias]

colData(sce) <- sample_metadata.dt %>% copy %>%
  .[,sample:=NULL] %>% setnames("sample_before_split","sample") %>%
  tibble::column_to_rownames("alias") %>% DataFrame

###############
## Aggregate ##
###############

sce_aggr <- pseudobulk_sce_fn(
  x = sce,
  assay = "counts",
  by = "sample",
  fun = "sum",
  scale = FALSE # Should pseudo-bulks be scaled with the effective library size & multiplied by 1M?
)

assayNames(sce_aggr) <- "counts"

# colData
stopifnot(sort(sample_metadata_aggregated.dt$alias)==sort(colnames(sce_aggr)))

colData(sce_aggr) <- sample_metadata_aggregated.dt %>%
  tibble::column_to_rownames("alias") %>% DataFrame

head(colData(sce_aggr))

# rowData
rowData(sce_aggr) <- rowData(sce)

##########
## Save ##
##########

fwrite(sample_metadata_aggregated.dt, paste0(io$basedir,"/cell_metadata_aggregated.txt"), quote=F, na="NA", sep="\t")
saveRDS(sce_aggr, file.path(io$outdir,"gene_SingleCellExperiment.rds"))



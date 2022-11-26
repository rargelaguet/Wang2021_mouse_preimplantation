here::i_am("rna/dimensionality_reduction/dimensionality_reduction_sce.R")

source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

## START TEST ##
# args$stages <- "E4.5"
# args$features <- 2500
# args$npcs <- 30
# args$colour_by <- c("plate", "celltype", "celltype2", "celltype3", "stage", "nFeature_RNA")
# args$colour_by <- c("celltype")
# opts$vars_to_regress <- NULL # c("nFeature_RNA")
# opts$batch_correction <- NULL
# args$remove_ExE_cells <- FALSE
# args$n_neighbors <- 25
# args$min_dist <- 0.5
## END TEST ##

# I/O
io$outdir <- paste0(io$basedir,"/results/rna/dimensionality_reduction"); dir.create(io$outdir, showWarnings = F)

# Options
opts$celltypes <- c(
  "zygote",
  "2cell", 
  "early_4cell", 
  "late_4cell", 
  "8cell", 
  "16cell", 
  "ICM", 
  "TE" 
)
opts$seed <- 42

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE & celltype%in%opts$celltypes]

table(cell_metadata.dt$celltype)

###################
## Sanity checks ##
###################

if (length(opts$batch_correction)>0) {
  stopifnot(opts$batch_correction%in%colnames(cell_metadata.dt))
  if (length(unique(cell_metadata.dt[[opts$batch_correction]]))==1) {
    message(sprintf("There is a single level for %s, no batch correction applied",opts$batch_correction))
    opts$batch_correction <- NULL
  } else {
    library(batchelor)
  }
}

if (length(opts$vars_to_regress)>0) {
  stopifnot(opts$vars_to_regress%in%colnames(cell_metadata.dt))
}

###############
## Load data ##
###############

# Load RNA expression data as SingleCellExperiment object
sce <- load_SingleCellExperiment(io$rna.sce, cells=cell_metadata.dt$cell, normalise = TRUE)

# Add sample metadata as colData
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

#######################
## Feature selection ##
#######################

decomp <- modelGeneVar(sce)
# if (length(opts$batch_correction)>0) {
#   decomp <- modelGeneVar(sce, block=colData(sce)[[opts$batch_correction]])
# } else {
#   decomp <- modelGeneVar(sce)
# }
decomp <- decomp[decomp$mean > 0.01,]
hvgs <- decomp[order(decomp$FDR),] %>% head(n=2500) %>% rownames

# Subset SingleCellExperiment
sce_filt <- sce[hvgs,]

############################
## Regress out covariates ##
############################

if (length(opts$vars_to_regress)>0) {
  print(sprintf("Regressing out variables: %s", paste(opts$vars_to_regress,collapse=" ")))
  logcounts(sce_filt) <- RegressOutMatrix(
    mtx = logcounts(sce_filt),
    covariates = colData(sce_filt)[,opts$vars_to_regress,drop=F]
  )
}

############################
## PCA + Batch correction ##
############################

# outfile <- sprintf("%s/%s_pca_features%d_pcs%d.txt.gz",io$outdir, paste(args$samples,collapse="-"), args$features, args$npcs)
sce_filt <- runPCA(sce_filt, ncomponents = 15, ntop=2500)

# if (length(opts$batch_correction)>0) {
#   suppressPackageStartupMessages(library(batchelor))
#   print(sprintf("Applying MNN batch correction for variable: %s", opts$batch_correction))
#   outfile <- sprintf("%s/%s_pca_features%d_pcs%d_batchcorrectionby%s.txt.gz",io$outdir, paste(args$samples,collapse="-"), args$features, args$npcs,paste(opts$batch_correction,collapse="-"))
#   pca <- multiBatchPCA(sce_filt, batch = colData(sce_filt)[[opts$batch_correction]], d = args$npcs)
#   pca.corrected <- reducedMNN(pca)$corrected
#   colnames(pca.corrected) <- paste0("PC",1:ncol(pca.corrected))
#   reducedDim(sce_filt, "PCA") <- pca.corrected
# } else {
#   # outfile <- sprintf("%s/%s_pca_features%d_pcs%d.txt.gz",io$outdir args$features, args$npcs)
#   sce_filt <- runPCA(sce_filt, ncomponents = args$npcs, ntop=args$features)
# }

# Save PCA coordinates
# pca.dt <- reducedDim(sce_filt,"PCA") %>% round(3) %>% as.data.table(keep.rownames = T) %>% setnames("rn","cell")
# fwrite(pca.dt, sprintf("%s/pca_features%d_pcs%d.txt.gz",io$outdir, args$features, args$npcs))

# Plot PCA
plotPCA(sce_filt, colour_by="celltype", point_size=3, point_alpha=1)

###########
## t-SNE ##
###########

sce_filt <- runTSNE(sce_filt, dimred="PCA", perplexity=5)
plotTSNE(sce_filt, colour_by="celltype", point_size=3, point_alpha=1)

##########
## UMAP ##
##########

# Run
set.seed(opts$seed)
sce_filt <- runUMAP(sce_filt, dimred="PCA", n_neighbors = 15, min_dist = 0.3)

# Fetch UMAP coordinates
umap.dt <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>%
  .[,cell:=colnames(sce_filt)] %>%
  setnames(c("UMAP1","UMAP2","cell"))

# Save UMAP coordinates
# fwrite(umap.dt, sprintf("%s/umap_features%d_pcs%d_neigh%d_dist%s.txt.gz",io$outdir, args$features, args$npcs, args$n_neighbors, args$min_dist))

##########
## Plot ##
##########

to.plot <- reducedDim(sce_filt,"UMAP") %>% as.data.table %>% 
  .[,cell:=colnames(sce_filt)] %>%
  merge(cell_metadata.dt, by="cell")

ggplot(to.plot, aes_string(x="V1", y="V2", fill="celltype")) +
  geom_point(size=3, shape=21, stroke=0.15) +
  theme_classic() +
  ggplot_theme_NoAxes() + theme(
    legend.position  ="right",
    legend.title=element_blank()
  )

# pdf(outfile, width=9, height=5)
# print(p)
# dev.off()


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(argparse))

#########
## I/O ##
#########

io <- list()
if (Sys.info()[['nodename']]=="rargelaguet.local") {
  io$basedir <- "/Users/rargelaguet/data/Wang2021_mouse_preimplantation"
  io$gene.metadata <- io$gene_metadata <- "/Users/rargelaguet/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
  io$mm10.genome <- "/Users/rargelaguet/data/mm10_sequence/mm10.genome"
} else if (Sys.info()[['nodename']]=="BI2404M") {
  io$basedir <- "/Users/argelagr/data/Wang2021_mouse_preimplantation"
  io$gene.metadata <- io$gene_metadata <- "/Users/argelagr/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
  io$mm10.genome <- "/Users/argelagr/data/mm10_sequence/mm10.genome"
} else if (grepl("pebble|headstone", Sys.info()['nodename'])) {
    if (grepl("argelag", Sys.info()['effective_user'])) {
      io$basedir <- "/bi/group/reik/ricard/data/Wang2021_mouse_preimplantation"
      io$atlas.basedir <- "/bi/group/reik/ricard/data/pijuansala2019_gastrulation10x"
      io$gene.metadata <- io$gene_metadata <- "/bi/group/reik/ricard/data/ensembl/mouse/v87/BioMart/all_genes/Mmusculus_genes_BioMart.87.txt"
      io$mm10.genome <- "/bi/group/reik/ricard/data/mm10_sequence/mm10.genome"
    }
} else {
  stop("Computer not recognised")
}

io$metadata <- paste0(io$basedir,"/cell_metadata.txt.gz")

# Methylation
io$met_data_raw <- paste0(io$basedir,"/processed/met/cpg_level")
io$met_data_pseudobulk_raw <- paste0(io$basedir,"/processed/met/cpg_level/pseudobulk")
io$met_data_parsed <- paste0(io$basedir,"/processed/met/feature_level")

# Accessibility
io$acc_data_raw <- paste0(io$basedir,"/processed/acc/gpc_level")
io$acc_data_pseudobulk_raw <- paste0(io$basedir,"/processed/acc/gpc_level/pseudobulk")
io$acc_data_parsed <- paste0(io$basedir,"/processed/acc/feature_level")

# RNA
io$rna.sce <- paste0(io$basedir,"/processed/rna/SingleCellExperiment.rds")

# Other
io$features.dir <- paste0(io$basedir,"/features/genomic_contexts")

#############
## Options ##
#############

opts <- list()

opts$celltypes <- opts$stages <- c(
  "zygote",
  "2cell", 
  "early_4cell", 
  "late_4cell", 
  "8cell", 
  "16cell", 
  "ICM", 
  "TE" 
)

# opts$celltype.colors = c(
#   "zygote" = "",
#   "2cell" = "",
#   "early_4cell" = "",
#   "late_4cell" = "",
#   "8cell" = "",
#   "16cell" = "",
#   "ICM" = "",
#   "TE" = "" 
# )


opts$chr <- paste0("chr",c(1:19,"X","Y"))

# opts$stages <- c("E3.5", "E4.5", "E5.5", "E6.5", "E7.5", "E8.5")
# opts$stage.colors <- viridis::viridis(n=length(opts$stages))
# names(opts$stage.colors) <- rev(opts$stages)

##########################
## Load sample metadata ##
##########################

# factor.cols <- c("id_rna","id_met","id_acc","stage","lineage","lab","plate","embryo")

# sample_metadata <- fread(io$metadata) %>% 
#   .[,stage_lineage:=paste(stage,lineage10x_2,sep="_")]
  # %>% .[,(factor.cols):=lapply(.SD, as.factor),.SDcols=(factor.cols)] %>% droplevels

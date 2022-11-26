source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# io$metadata <- file.path(io$basedir,"results/mapping/sample_metadata_after_mapping.txt.gz")
io$outfile <- paste0(io$basedir,"/shiny/cell_metadata.txt.gz")

########################
## Load cell metadata ##
########################

cell_metadata.dt <- fread(io$metadata) %>%
  .[pass_rnaQC==TRUE] %>%
  .[,c("cell","batch","celltype","embryo_id","sex","num_genes")]

table(cell_metadata.dt$celltype)

##########
## Save ##
##########

fwrite(cell_metadata.dt, io$outfile, sep="\t", na="NA", quote=F)

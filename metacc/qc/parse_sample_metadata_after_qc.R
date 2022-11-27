here::here("metacc/qc/qc.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))


######################
## Define arguments ##
######################

p <- ArgumentParser(description='')
p$add_argument('--metadata_met',  type="character",              help='Cell metadata file')
p$add_argument('--metadata_acc',  type="character",              help='Cell metadata file')
p$add_argument('--outfile',  type="character",              help='Output file')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define settings ##
#####################

## START TEST ##
# args <- list()
# args$metadata_met <- file.path(io$basedir,"results/met/qc/sample_metadata_after_met_qc.txt.gz")
# args$metadata_acc <- file.path(io$basedir,"results/acc/qc/sample_metadata_after_acc_qc.txt.gz")
# args$outfile  <- file.path(io$basedir,"results/metacc/qc/sample_metadata_after_metacc_qc.txt.gz")
## END TEST ##

###################
## Load metadata ##
###################

cell_metadata_met.dt <- fread(args$metadata_met)
cell_metadata_acc.dt <- fread(args$metadata_acc)

###########
## Merge ##
###########

cell_metadata.dt <- merge(cell_metadata_met.dt[,c("cell","pass_metQC","nCG","met_rate")], cell_metadata_acc.dt, by="cell")

print(head(cell_metadata.dt))

##########
## Save ##
##########

fwrite(cell_metadata.dt, args$outfile, sep="\t", na = "NA", quote=F)

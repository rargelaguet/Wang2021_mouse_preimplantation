# here::here("metacc/quantify_feature_level/quantify_feature_level.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# I/O
io$indir <- file.path(io$basedir,"original/met")
io$outdir <- file.path(io$basedir,"processed/met/cpg_level")
io$metadata <- file.path(io$basedir,"results/met/qc/sample_metadata_after_met_qc.txt.gz")

# Options
opts$context <- "CG"

# Sanity checks
stopifnot(opts$context %in% c("CG","GC"))

# chr10   3000490 -       0       1       0       TCG
# chr10   3000894 +       0       1       0       ACG
# chr10   3000895 -       1       0       1       TCG
# chr10   3000924 +       0       1       0       TCG
# chr10   3000925 -       0       1       0       ACG

##################
## Define cells ##
##################

# Define cells
if (opts$context=="CG") {
  cells <- fread(io$metadata) %>% .[!is.na(id_met),id_met]
} else if (opts$context=="GC") {
  cells <- fread(io$metadata) %>% .[!is.na(id_acc),id_acc]
}

cells <- head(cells,n=3)

##################
## Process data ##
##################

for (i in cells) {
  print(i)

  data.dt <- fread(sprintf("%s/%s.gz",io$indir,i), sep="\t", verbose=F, showProgress=F) %>%
    setnames(c("chr","pos","strand","met_reads","nonmet_reads","rate","context"))
    .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
    .[chr%in%opts$chr] %>%
    [,rate:=round(rate,2)] %>%
    # .[rate!=0.5] %>% .[rate>0.5,rate:=1] %>% .[rate<0.5,rate:=0] %>%
    .[,c("chr","pos","rate")]


  # Sanity check
  # stopifnot(all(dat_sample$rate %in% c(0,1)))

      
  # Save results
  fwrite(data.dt, file.path(io$outdir,sprintf("%s.tsv.gz",i)), quote=FALSE, sep="\t", col.names=FALSE)
}

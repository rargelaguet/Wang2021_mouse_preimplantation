# here::here("metacc/quantify_feature_level/quantify_feature_level.R")

# Load default settings
source(here::here("settings.R"))
source(here::here("utils.R"))

#####################
## Define settings ##
#####################

# io$indir <- file.path(io$basedir,"original/met")
# io$outdir <- file.path(io$basedir,"processed/met/cpg_level"); dir.create(io$outdir, showWarnings=F)
# opts$context <- "CG"

io$indir <- file.path(io$basedir,"original/acc")
io$outdir <- file.path(io$basedir,"processed/acc/gpc_level"); dir.create(io$outdir, showWarnings=F)
opts$context <- "GC"

##################
## Define cells ##
##################

cell_metadata.dt <- fread(io$metadata)
if (opts$context=="CG") {
  cells <- cell_metadata.dt %>% .[!is.na(id_met),id_met]
} else if (opts$context=="GC") {
  cells <- cell_metadata.dt %>% .[!is.na(id_acc),id_acc]
}

# cells <- head(cells,n=3)
stopifnot(length(setdiff(cells,gsub(".gz","",list.files(io$indir))))==0)
# cell_metadata.dt[!id_acc%in%gsub(".gz","",list.files(io$indir))]

cells <- setdiff(list.files(io$indir) %>% gsub(".gz","",.), list.files(io$outdir) %>% gsub(".tsv.gz","",.))

##################
## Process data ##
##################

for (i in cells) {
  print(i)

  data.dt <- fread(sprintf("%s/%s.gz",io$indir,i), sep="\t", verbose=F, showProgress=F) %>%
    setnames(c("chr","pos","strand","met_reads","nonmet_reads","rate","context")) %>%
    .[,chr:=ifelse(grepl("chr",chr),chr,paste0("chr",chr))] %>%
    .[chr%in%opts$chr] %>%
    .[,rate:=round(rate,2)] %>%
    # .[rate!=0.5] %>% .[rate>0.5,rate:=1] %>% .[rate<0.5,rate:=0] %>%
    .[,c("chr","pos","rate")]


  # Sanity check
  # stopifnot(all(dat_sample$rate %in% c(0,1)))

      
  # Save results
  fwrite(data.dt, file.path(io$outdir,sprintf("%s.tsv.gz",i)), quote=FALSE, sep="\t", na="NA", col.names=TRUE)
}

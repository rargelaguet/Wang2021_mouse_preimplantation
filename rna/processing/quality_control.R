here::i_am("rna/processing/quality_control.R")

source(here::here("settings.R"))

#####################
## Define settings ##
#####################

# I/O
io$metadata <- file.path(io$basedir,"cell_metadata_aggregated.txt")
io$outdir <- paste0(io$basedir,"/results/rna/qc"); dir.create(io$outdir, showWarnings=F, recursive=T)

# Options
opts$min_reads <- 500000
opts$max_mit_percentage <- 15
opts$max_rib_percentage <- 15

###################
## Load metadata ##
###################

sample_metadata.dt <- fread(io$metadata) %>% .[!is.na(alias)]

# sample_metadata.dt <- sample_metadata.dt[!alias%in%c("GSM4055925_Zygote1","GSM4056174_8-cell")]

##########################
## Load expression data ##
##########################

sce <- readRDS(io$gene.sce)[,sample_metadata.dt$alias]
rownames(sce) <- rowData(sce)$symbol

########################
## Calculate QC stats ##
########################

qc_stats.dt <- data.table(
  alias = colnames(sce),
  number_genes = colSums(assay(sce,"counts")>0),
  number_reads = colSums(assay(sce,"counts")),
  mit_percentage = 100*colSums(assay(sce[grep("^mt-",rownames(sce), value=T),],"counts")) / colSums(assay(sce,"counts")),
  rib_percentage = 100*colSums(assay(sce[grep("^Rp[l|s]",rownames(sce), value=T),],"counts")) / colSums(assay(sce,"counts")) %>% round(2)
) %>% .[,c("mit_percentage","rib_percentage"):=list(round(mit_percentage,2),round(rib_percentage,2))]

sample_metadata.dt <- sample_metadata.dt %>% merge(qc_stats.dt,by="alias",all.x=T)

##############
## QC calls ##
##############

sample_metadata.dt %>%
  .[,pass_rnaQC:=number_reads>=opts$min_reads & mit_percentage<=opts$max_mit_percentage & rib_percentage<=opts$max_rib_percentage]

table(sample_metadata.dt$pass_rnaQC)

#####################
## Plot QC metrics ##
#####################

to.plot <- sample_metadata.dt %>% 
  # .[pass_rnaQC==TRUE] %>%
    melt(id.vars=c("sample","pass_rnaQC"), measure.vars=c("number_reads","mit_percentage","rib_percentage"))

facet.labels <- c("number_reads" = "Num. of reads", "mit_percentage" = "Mitochondrial %", "rib_percentage" = "Ribosomal %")
    
## Box plot 

p <- ggplot(to.plot, aes(x=sample, y=value, fill=pass_rnaQC)) +
    geom_bar(stat="identity") +
    facet_wrap(~variable, scales="free_y", nrow=1, labeller = as_labeller(facet.labels)) +
    # scale_fill_manual(values=opts$day.colors) +
    guides(x = guide_axis(angle = 90)) +
    theme_classic() +
    theme(
        axis.text.y = element_text(colour="black",size=rel(1)),
        # axis.text.x = element_text(colour="black",size=rel(0.45)),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
    )

pdf(sprintf("%s/qc_metrics_boxplot.pdf",io$outdir), width=11, height=6)
print(p)
dev.off()

## histogram 

tmp <- data.table(
    variable = c("number_reads", "mit_percentage", "rib_percentage"),
    value = c(opts$min_reads, opts$max_mit_percentage, opts$max_rib_percentage)
)
p <- gghistogram(to.plot, x="value", bins=20) +
    geom_vline(aes(xintercept=value), linetype="dashed", data=tmp) +
    facet_wrap(~variable, scales="free", nrow=1, labeller = as_labeller(facet.labels)) +
    theme(
        axis.text =  element_text(size=rel(0.8)),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.75))
    )
    
pdf(sprintf("%s/qc_metrics_histogram.pdf",io$outdir), width=13, height=6)
print(p)
dev.off()

##########
## Save ##
##########

fwrite(qc_stats.dt, paste0(io$outdir,"/qc_stats.txt"), quote=F, na="NA", sep="\t")
fwrite(sample_metadata.dt, paste0(io$outdir,"/cell_metadata_after_qc.txt"), quote=F, na="NA", sep="\t")


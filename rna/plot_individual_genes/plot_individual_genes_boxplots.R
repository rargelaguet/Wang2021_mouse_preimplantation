here::i_am("rna/plot_individual_genes/plot_individual_genes_boxplots.R")

#####################
## Define settings ##
#####################

source(here::here("settings.R"))
source(here::here("utils.R"))

## I/O ##

# io$metadata <- file.path(io$basedir,"results/rna/celltype_assignment/cell_metadata.dt_after_celltype_rename.txt.gz")
io$outdir <- paste0(io$basedir,"/results/rna/individual_genes"); dir.create(io$outdir, showWarnings = F)

## Define options ##

opts$celltypes = c(
  "zygote",
  "2cell", 
  "early_4cell", 
  "late_4cell", 
  "8cell", 
  "16cell", 
  "ICM", 
  "TE" 
)

##########################
## Load sample metadata ##
##########################

cell_metadata.dt <- fread(io$metadata) %>% 
  .[pass_rnaQC==T & celltype%in%opts$celltypes] %>%
  .[,celltype:=factor(celltype,levels=opts$celltypes)]

table(cell_metadata.dt$celltype)

###############
## Load data ##
###############

sce <- load_SingleCellExperiment(
  file = io$rna.sce, 
  cells = cell_metadata.dt$cell, 
  normalise = TRUE, 
  remove_non_expressed_genes = TRUE
)

# Add sample metadata as colData
colData(sce) <- cell_metadata.dt %>% tibble::column_to_rownames("cell") %>% DataFrame

#############################
## Plot one gene at a time ##
#############################

# genes.to.plot <- c("Eomes","Lefty1","Lefty2","Nodal")
genes.to.plot <- grep("Tet", rownames(sce),value=T)

for (i in genes.to.plot) {
  
    to.plot <- data.table(
      cell = colnames(sce),
      expr = logcounts(sce)[i,]
    ) %>% merge(cell_metadata.dt, by="cell")
    
    p <- ggplot(to.plot, aes(x=celltype, y=expr, fill=celltype)) +
      geom_violin(scale = "width", alpha=0.8) +
      geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.8) +
      geom_jitter(shape=21, size=1.5, alpha=0.5, width=0.05) +
      # scale_fill_manual(values=opts$stage.colors, drop=F) +
      stat_summary(fun.data = give.n, geom = "text", size=2.5) +
      # facet_wrap(~stage, nrow=1, scales="free_x") +
      theme_classic() +
      labs(title=i, x="",y=sprintf("%s expression",i)) +
      theme(
        strip.text = element_text(size=rel(0.85)),
        # plot.title = element_text(hjust = 0.5, size=rel(1.1), color="black"),
        plot.title = element_blank(),
        axis.text.x = element_text(colour="black",size=rel(1.25)),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour="black",size=rel(1.0)),
        axis.title.y = element_text(colour="black",size=rel(1.0)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(0.85))
      )
    
      pdf(sprintf("%s/%s_expr.pdf",io$outdir,i), width=7.5, height=5)
      print(p)
      dev.off()
}


##########################################
## Plot multiple genes at the same time ##
##########################################

# genes.to.plot <- c("Eomes","Lefty1","Lefty2","Nodal")
# genes.to.plot <- grep("Tet", rownames(sce),value=T)
genes.to.plot <- c("Dnmt1","Dnmt3a","Dnmt3b")

to.plot <- logcounts(sce[genes.to.plot]) %>% as.data.table(keep.rownames = T) %>%
  setnames("rn","gene") %>%
  melt(id.vars=("gene"), variable.name="cell", value.name = "expr") %>% 
  merge(cell_metadata.dt, by="cell")
  
p <- ggplot(to.plot, aes(x=celltype, y=expr, fill=celltype)) +
  geom_violin(alpha=0.30) +
  geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.8) +
  # geom_jitter(shape=21, size=0.5, alpha=0.4, width=0.03, stroke=0.05) +
  ggrastr::geom_jitter_rast(shape=21, size=0.5, alpha=0.4, width=0.03, stroke=0.05) +
  # scale_fill_manual(values=opts$stage.colors, drop=F) +
  facet_wrap(~gene, nrow=1, scales="free_x") +
  theme_classic() +
  labs(x="", y="Gene expression") +
  theme(
    strip.text = element_text(size=rel(1.25)),
    strip.background = element_blank(),
    axis.text.x = element_text(colour="black",size=rel(1)),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour="black",size=rel(1.0)),
    axis.title.y = element_text(colour="black",size=rel(1.0)),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size=rel(0.85))
  )
  
# pdf(sprintf("%s/%s_expr.pdf",io$outdir,i), width=7.5, height=4)
print(p)
# dev.off()


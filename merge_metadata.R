
foo <- fread("/Users/argelagr/data/Wang2021/original/20190826.updated.sample.metadata.tsv") %>%
  setnames("cellid","cell")
bar <- fread("/Users/argelagr/data/Wang2021/cell_metadata.txt.gz") %>% .[,id_met:=NULL]
  
foo[,id_met:=paste0(DNA.info,".ACG.TCG")]
foo[,id_acc:=paste0(DNA.info,".GCA.GCT.GCC")]

colnames(foo)
colnames(bar)

bar <- bar %>% merge(foo[,c("cell","id_met","id_acc")], by="cell")

fwrite(bar,"/Users/argelagr/data/Wang2021/cell_metadata.txt.gz", quote = F, sep="\t", na="NA")



cell_metadata.dt <- fread(file.path(io$basedir,"backup/cell_metadata.txt.gz"))

foo <- cell_metadata.dt[,c("cell", "id_met", "id_acc", "embryo_id", "batch", "celltype", "cluster", "cnv_info", "sex", "pass_rnaQC", "num_genes", "DC1", "DC2", "DC3", "DC.order")] %>%
  setnames("num_genes","nFeature_RNA") %>%
  setnames("embryo_id","embryo") %>%
  .[,c("DC1","DC2","DC3"):=list(round(DC1,3),round(DC2,3),round(DC3,3))] 

bar <- cell_metadata.dt[,c("cell","WCG_sites", "GCH_sites", "gch_mean_rate", "wcg_mean_rate", "chrM_gch_rate", "chrM_wcg_rate", "lambda_gch_rate", "lambda_wcg_rate", "gene_num")]

table(baz$celltype)

fwrite(foo,file.path(io$basedir,"cell_metadata.txt.gz"), quote = F, sep="\t", na="NA")
fwrite(bar,file.path(io$basedir,"backup/cell_metadata_met_acc_stats.txt.gz"), quote = F, sep="\t", na="NA")
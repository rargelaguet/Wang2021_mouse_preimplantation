
foo <- fread("/Users/argelagr/data/Wang2021/original/20190826.updated.sample.metadata.tsv") %>%
  setnames("cellid","cell")
bar <- fread("/Users/argelagr/data/Wang2021/cell_metadata.txt.gz") %>% .[,id_met:=NULL]
  
foo[,id_met:=paste0(DNA.info,".ACG.TCG")]
foo[,id_acc:=paste0(DNA.info,".GCA.GCT.GCC")]

colnames(foo)
colnames(bar)

bar <- bar %>% merge(foo[,c("cell","id_met","id_acc")], by="cell")

fwrite(bar,"/Users/argelagr/data/Wang2021/cell_metadata.txt.gz", quote = F, sep="\t", na="NA")

baz <- fread("/Users/argelagr/data/Wang2021/cell_metadata.txt")

table(baz$celltype)

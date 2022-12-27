
foo <- fread("/Users/argelagr/data/Wang2021_mouse_preimplantation/results/met/stats/sample_metadata_after_met_stats.txt.gz") 
bar <- fread("/Users/argelagr/data/Wang2021_mouse_preimplantation/cell_metadata.txt.gz") %>% .[,c("cell","embryo","batch","cnv_info","sex")]
  
bar <- bar %>% merge(foo, by="cell")

fwrite(bar, file.path(io$basedir,"results/met/stats/sample_metadata_after_met_stats.txt.gz"), quote = F, sep="\t", na="NA")

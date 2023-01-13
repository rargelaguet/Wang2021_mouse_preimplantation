source(here::here("settings.R"))

foo <- fread("/Users/argelagr/data/Wang2021_mouse_preimplantation/fastq/SraRunTable.txt") %>%
  .[,c("Run","GEO_Accession (exp)","develomental_stage")] %>% 
  setnames(c("sra","geo","stage"))

bar <- fread("/Users/argelagr/data/Wang2021_mouse_preimplantation/fastq/foo.txt")


foobar <- merge(foo,bar,by="geo")
rm(foo,bar)

fwrite(foobar, "/Users/argelagr/data/Wang2021_mouse_preimplantation/fastq/sra_sample_metadata.txt", sep="\t", quote=F, col.names = T)


# cell_metadata.dt <- fread(io$metadata) %>%
#   .[,id:=strsplit(id_met,"\\.") %>% map_chr(1)] %>%
#   .[,c("cell","id")] %>%
#   merge(geo2id,by="id",all.x=T)

# bar <- bar %>% merge(cell_metadata.dt,by="geo",all.x=T)

# fwrite(bar, "/Users/argelagr/data/Wang2021_mouse_preimplantation/fastq/sra_sample_metadata.txt", sep="\t", quote=F, col.names = F)

# baz <- fread("/Users/argelagr/data/Wang2021_mouse_preimplantation/fastq/sra_sample_metadata.txt")



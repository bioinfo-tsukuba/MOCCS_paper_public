totalization <- readRDS(url("https://figshare.com/ndownloader/files/34065686","rb")) #MOCCSout_hg38_all_qval_annotated.rds
#totalization <- readRDS("~/MOCCS_paper_public/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
ID_hard <- readRDS(paste0(annotation_path, "hg38_hard_filter_ID.rds"))

tmp2 <- totalization %>% filter(ID %in% ID_hard & Antigen_class == "TFs and others") %>% 
  select(ID, Antigen_class, Antigen, Cell_type_class, Cell_type) %>% distinct()
dim(tmp2)

print("Antigen")
length(unique(tmp2$Antigen))

print("Cell type class")
length(unique(tmp2$Cell_type_class))

print("Cell type")
length(unique(tmp2$Cell_type))
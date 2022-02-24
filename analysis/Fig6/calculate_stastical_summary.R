# SLE, Cell type class
target_phenotype <- "SLE"
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig6/"
top_x <- 100

library(tidyverse)
target_df <- readRDS(paste0(annotation_path, "result_output_binded_all/", target_phenotype, "_peak_rand_binded_all_qval.rds"))

totalization_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds"
totalization <- readRDS(totalization_path)
annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID") %>% filter(q_value < 0.05)
df1 <- df_phenotype_binded_all_selected_annotated %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% select(ID, p_value, Cell_type_class, Cell_type, Antigen, abs_dMOCCS2score) %>% distinct()

ID_snp_for_join <- df_phenotype_binded_all_selected_annotated %>% unite("ID_snp_for_join", c(ID, position)) %>% .$ID_snp_for_join %>% as.character()
df_phenotype_binded_all_selected_annotated2 <- df_phenotype_binded_all_selected_annotated %>% mutate(ID_snp_join = ID_snp_for_join)
df2 <- df_phenotype_binded_all_selected_annotated2 %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% group_by(ID, position, Antigen, Cell_type_class, Cell_type) %>% summarise(max_abs_dMOCCS2score = max(abs_dMOCCS2score))
df3 <- df2  %>% unite("ID_position", c(ID, position)) %>% arrange(desc(max_abs_dMOCCS2score))
df4 <- df3[1:top_x,]

target_CLC <- "Blood"
target_df <- df4 %>% filter(Cell_type_class == target_CLC)
target_CLC_snp_N <- nrow(target_df)

# CD, TFs
target_phenotype <- "CD"
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig6/"
top_x <- 50

library(tidyverse)
target_df <- readRDS(paste0(annotation_path, "result_output_binded_all/", target_phenotype, "_peak_rand_binded_all_qval.rds"))

totalization_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds"
totalization <- readRDS(totalization_path)
annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID") %>% filter(q_value < 0.05)
df1 <- df_phenotype_binded_all_selected_annotated %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% select(ID, p_value, Cell_type_class, Cell_type, Antigen, abs_dMOCCS2score) %>% distinct()

ID_snp_for_join <- df_phenotype_binded_all_selected_annotated %>% unite("ID_snp_for_join", c(ID, position)) %>% .$ID_snp_for_join %>% as.character()
df_phenotype_binded_all_selected_annotated2 <- df_phenotype_binded_all_selected_annotated %>% mutate(ID_snp_join = ID_snp_for_join)
df2 <- df_phenotype_binded_all_selected_annotated2 %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% group_by(ID, position, Antigen, Cell_type_class, Cell_type) %>% summarise(max_abs_dMOCCS2score = max(abs_dMOCCS2score))
df3 <- df2  %>% unite("ID_position", c(ID, position)) %>% arrange(desc(max_abs_dMOCCS2score))
df4 <- df3[1:top_x,]

target_TF <- "SPI1"
target_df <- df4 %>% filter(Antigen == target_TF)
target_TF_snp_N <- nrow(target_df)
target_TF_snp_N

target_TF <- "FOS"
target_df <- df4 %>% filter(Antigen == target_TF)
target_TF_snp_N <- nrow(target_df)
target_TF_snp_N



# rs17293632
target_phenotype <- "CD"
target_rs <- "rs17293632"
target_position <- "chr15_67150258"
top_x <- 10
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig6/"
target_df <- readRDS(paste0(annotation_path, "result_output_binded_all/", target_phenotype, "_peak_rand_binded_all_qval.rds"))

totalization_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds"
totalization <- readRDS(totalization_path)
annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID") %>% filter(q_value < 0.05)
df1 <- df_phenotype_binded_all_selected_annotated %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% select(ID, p_value, Cell_type_class, Cell_type, Antigen, abs_dMOCCS2score) %>% distinct()

ID_snp_for_join <- df_phenotype_binded_all_selected_annotated %>% 
  filter(position == target_position) %>% 
  unite("ID_position", c(ID, position)) %>% 
  .$ID_position %>% 
  as.character()
df_phenotype_binded_all_selected_annotated2 <- df_phenotype_binded_all_selected_annotated %>% 
  filter(position == target_position) %>%
  mutate(ID_position = ID_snp_for_join) %>% 
  drop_na(dMOCCS2score) 
df2 <- df_phenotype_binded_all_selected_annotated2 %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% 
  select(ID_position, position, abs_dMOCCS2score, ID, Antigen, Cell_type_class, Cell_type) %>%
  distinct() %>%
  group_by(ID, Antigen, Cell_type_class, position) %>% 
  summarise(max_abs_dMOCCS2score = max(abs_dMOCCS2score))
df3 <- df2 %>% unite("ID_position", c(ID, position)) %>% arrange(desc(max_abs_dMOCCS2score))
df4 <- df3[1:top_x,]

target_TF <- "FOS"
target_df <- df4 %>% filter(Antigen == target_TF)
target_TF_snp_N <- nrow(target_df)
target_TF_snp_N


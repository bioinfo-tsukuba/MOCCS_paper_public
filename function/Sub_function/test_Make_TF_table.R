system("wget https://figshare.com/ndownloader/files/34692295 -P ~/MOCCS_paper_public/data/Fig2/obj")
system("mv ~/MOCCS_paper_public/data/Fig2/obj/34692295 ~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_raw <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_fam <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
all_ns <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/all_ns.rds")
df_p_1 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_1.rds")
df_p_2 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_2.rds")
df_p_3 <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3.rds")
df_p_3_gp <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_p_3_gp.rds")

df1 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_1.txt", col_names = F)
df2 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_2.txt", col_names = F)
df3 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_3.txt", col_names = F)
df4 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_4.txt", col_names = F)
df5 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_5.txt", col_names = F)
df6 <- read_tsv("~/MOCCS_paper_public/data/TableS1/Supplementary_table_S1_parts_6.txt", col_names = F)

TF_list <- c(df1$X1, df3$X1, df5$X1)
length(TF_list)
pval_list <- c(df2$X1, df4$X1, df6$X1)
length(pval_list)
df7 <- tibble(TF = TF_list, pval = pval_list)
dim(df7)
head(df7)

ctc_num_list <- c()
for (target_ant in TF_list) {
  target_flag <- df_p_3_gp$ID1_Antigen == target_ant & df_p_3_gp$ID2_Antigen == target_ant
  target_df <- df_p_3_gp[target_flag, ]
  
  tmp1 <- target_df %>% select(ID1, ID1_Cell_type_class)
  tmp2 <- target_df %>% select(ID2, ID2_Cell_type_class)
  colnames(tmp1) <- c("ID", "Cell_type_class")
  colnames(tmp2) <- c("ID", "Cell_type_class")
  tmp3 <- rbind(tmp1, tmp2) %>% distinct()
  tgt_ctc_num <- tmp3$Cell_type_class %>% unique() %>% length()
  ctc_num_list <- c(ctc_num_list, tgt_ctc_num)
}

df8 <- df7 %>% mutate(ctc_num = ctc_num_list)
head(df8)

df9 <- df8 %>% mutate(ctc_depedent = ifelse(pval < 0.05, "Y", "N"))

# Add TF family


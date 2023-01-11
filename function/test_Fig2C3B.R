system("wget https://figshare.com/ndownloader/files/34692295 -P ~/MOCCS_paper_public/data/Fig2/obj")
system("mv ~/MOCCS_paper_public/data/Fig2/obj/34692295 ~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_raw <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_fam <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
all_ns <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/all_ns.rds")

# Fig. 2C and 3B: UMAP
source("~/MOCCS_paper_public/function/sub_function/Umap_df.R")
UMAP_df(df_raw, df_fam, all_ns)

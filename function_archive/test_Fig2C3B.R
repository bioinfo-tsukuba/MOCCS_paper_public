system("wget https://figshare.com/ndownloader/files/34692295 -P ~/MOCCS_paper_public/data/Fig2/obj")
system("mv ~/MOCCS_paper_public/data/Fig2/obj/34692295 ~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_raw <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_raw.rds")
df_fam <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/df_fam.rds")
all_ns <- readRDS("~/MOCCS_paper_public/data/Fig2/obj/all_ns.rds")

# Fig. 2C and 3B: UMAP
source("~/MOCCS_paper_public/function/sub_function/Umap_df.R")
UMAP_df(df_raw, df_fam, all_ns)

# UMAP_df.Rの中身
df_umap_2_Antigen <- df_umap_2[df_umap_2$Plot == "Antigen", ]
desc_list <- df_umap_2_Antigen %>% group_by(Antigen) %>% summarise(n = n()) %>% arrange(desc(n)) %>% .$Antigen
top_antigen <- desc_list[1:15] %>% as.character()
tmp1 <-  df_umap_2_Antigen %>% filter(!Annotation %in% top_antigen) 
tmp2 <- tmp1 %>% select(-Annotation) %>% mutate(Annotation = rep("others", nrow(tmp1)))
tmp3 <- df_umap_2_Antigen %>% filter(Annotation %in% top_antigen) 
df_umap_2_Antigen_new <- rbind(tmp2, tmp3)


library("colorspace")
p_umap_1_new <- ggplot2::ggplot(df_umap_2_Antigen_new[df_umap_2_Antigen_new$Plot == "Antigen", ], aes(x = UMAP1, y = UMAP2, color = Annotation)) +
  geom_point() +
  theme(aspect.ratio = 1.0) +
  theme_bw() +
  #theme(legend.position = "none") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(
    #values = sequential_hcl(30, h = c(260, 60), c = 60, l = c(40, 95), power = 1)
    values = c("#ff4b00","#fff100","#03af7a", "#005aff","#4dc4ff","#ff8082","#f6aa00","#990099","#804000", "#c8c8cb",
    "#ffff80", "#d8f255", "#bfe4ff", "#ffca80", "#77d9a8", "#c9ace6")
    #values = c("red", "orange", "yellow", "green", "lightblue", "blue", "purple",  "gray", "pink","red4", "black", "brown")
    #na.value = "gray",
    ## Add label for NA (here, the label is `none`)
    ## Replace `...` with your labels for non-NA levels below
    #labels = c(..., "none") 
  )+
  facet_wrap(~ Plot)

p_umap_1_new 
  
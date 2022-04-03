##################
### Figure1 ###
#################
## Fig. 1B --------
source("~/MOCCS_paper_public/function/Fig1B_plot.R")
target_ID_Fig1B <- "SRX1156473" #GATA3
# /Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1にhg38 target_ID_Fig1BのMOCCS outputを移してから以下を実行
Fig1B <- Fig1B_plot(target_ID_Fig1B)
plot(Fig1B)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig1/Fig1B.pdf"), Fig1B, width = 7, height = 7)


## Fig. 1C --------
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1C_plot.R")
target_TF <- "JUN"
load <- "local" #local or figshare
filter <- "all" #soft or hard or all
Fig1C_plot(target_TF, load, filter)

## Fig. 1D --------
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1D_plot.R")
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/"
Fig1D <- Fig1D_plot(annotation_path)
plot(Fig1D[[1]])
plot(Fig1D[[2]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig1/Fig1D_CL.pdf"), Fig1D[[1]], width = 7, height = 7)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig1/Fig1D_TF.pdf"), Fig1D[[2]], width = 7, height = 7)

## Fig. 1E --------
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1E_plot.R")
target_TF <- "CTCF"
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/"
Fig1E <- Fig1E_plot(target_TF, annotation_path)
plot(Fig1E)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig1/Fig1E_", target_TF,".pdf"), Fig1E, width = 7, height = 7)


## Fig. 1F --------
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1F_plot.R")
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/"
Fig1F <- Fig1F_plot(annotation_path)
plot(Fig1F)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig1/Fig1F.pdf"), Fig1F, width = 7, height = 7)


## Fig. 1G --------
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1G_plot.R")
target_TF <- "CTCF"
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/"
Fig1G <- Fig1G_plot(target_TF, annotation_path)
plot(Fig1G)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig1/Fig1G.pdf"), Fig1G, width = 7, height = 7)



##################
### Figure2-3 ###
#################
source("~/MOCCS_paper_public/function/Fig2_3_plot.R")
Fig2_3_plot(calc_opt = FALSE)

##################
### Figure4 ###
#################
## Fig. 4B --------
library(patchwork)
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig4B_plot.R")
simulation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig4/"
Fig4B <- Fig4B_plot(simulation_path)
Fig4B_plot <- Fig4B[[1]] / Fig4B[[2]] / Fig4B[[3]]
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig4/Fig4B_all.pdf"), Fig4B_plot, width = 14, height = 14)


## Fig. 4C --------
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig4C_plot.R")
#target_ID1 <- "SRX150600" # JUN K562
#target_ID2 <- "SRX186614" # JUN K562

target_ID1 <- "SRX150600" # JUN K562
target_ID2 <- "SRX150358" # JUN HUVEC

#target_ID1 <- "SRX150481" # EP300 K562
#target_ID2 <- "SRX150573" # EP300 K562

#target_ID1 <- "SRX150481" # EP300 K562
#target_ID2 <- "SRX2343934" # EP300 Cardiovascular Large pulmonary artery endothelial cells

target_TF <- "JUN"
#target_TF <- "EP300"

Fig4C <- Fig4C_plot(target_ID1, target_ID2, target_TF)
plot(Fig4C)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig4/Fig4C_", target_TF, "_", target_ID1, "_", target_ID2, ".pdf"), Fig4C, width = 7, height = 7)


## Fig. 4D --------
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig4D_plot.R")
#target_ID1 <- "SRX150546" #JUN
#target_ID2 <- "SRX029087" #FOS

target_ID1 <- "SRX150546" #JUN
target_ID2 <- "SRX190276" #CTCF

target_CT <- "K-562"
Fig4D <- Fig4D_plot(target_ID1, target_ID2, target_CT)
plot(Fig4D)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig4/Fig4D_", target_CT, "_", target_ID1, "_", target_ID2, ".pdf"), Fig4D, width = 7, height = 7)




##################
### Figure5 ###
#################
## Fig. 5C --------
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig5C_plot.R")
target_TF <- "HOXB13"
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig5/"
Fig5C <- Fig5C_plot(target_TF, annotation_path)
plot(Fig5C)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig5/Fig5C_", target_TF,".pdf"), Fig5C, width = 7, height = 7)


## Fig. 5E --------
### scatter plot and bar plot
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig5E_plot.R")
path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig5/"
target_TF <- "FOXA1"
Fig5E <- Fig5E_plot(target_TF, path)
plot(Fig5E[[1]])
plot(Fig5E[[2]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig5/Fig5E_scatter_", target_TF,".pdf"), Fig5E[[1]], width = 7, height = 7)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig5/Fig5E_bar_plot.pdf"), Fig5E[[2]], width = 7, height = 7)



##################
### Figure6 ###
#################
## Fig. 6B --------
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig6B_plot.R")
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig6/"
Fig6B <- Fig6B_plot(annotation_path)
plot(Fig6B)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig6/Fig6B.pdf"), Fig6B, width = 7, height = 7)


## Fig. 6C --------
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS_paper_public/function/Fig6C_plot.R")
target_phenotype <- "CD"
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig6/"
threshold <- 150
Fig6C <- Fig6C_plot(target_phenotype, annotation_path, threshold)
plot(Fig6C[[1]])
plot(Fig6C[[2]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig6/Fig6C_", target_phenotype,"_CL.pdf"), Fig6C[[1]], width = 7, height = 7)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig6/Fig6C_", target_phenotype,"_TF.pdf"), Fig6C[[2]], width = 7, height = 7)


## Fig. 6D --------
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS_paper_public/function/Fig6D_plot.R")
target_phenotype <- "CD"

target_rs <- "rs17293632"
target_position <- "chr15_67150258"

#target_rs <- "rs1057233"
#target_position <- "chr11_47354897"

#target_rs <- "rs10769258"
#target_position <- "chr11_47369488"

#target_rs <- "rs4752829"
#target_position <- "chr11_47375103"

annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig6/"
threshold <- 100
Fig6D <- Fig6D_plot(target_phenotype, target_rs, target_position ,annotation_path, threshold)
plot(Fig6D[[1]])
plot(Fig6D[[2]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig6/Fig6D_", target_phenotype,"_CL.pdf"), Fig6D[[1]], width = 7, height = 7)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig6/Fig6D_", target_phenotype,"_TF.pdf"), Fig6D[[2]], width = 7, height = 7)

### IBDの場合
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig6D_plot_2.R")
target_phenotype <- "IBD"

target_rs <- "rs17293632"
target_position <- "chr15_67150258"

annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig6/"
threshold <- 100
Fig6D <- Fig6D_plot(target_phenotype, target_rs, target_position ,annotation_path, threshold)
plot(Fig6D[[1]])
plot(Fig6D[[2]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig6/Fig6D_", target_phenotype,"_CL.pdf"), Fig6D[[1]], width = 14, height = 14)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig6/Fig6D_", target_phenotype,"_TF.pdf"), Fig6D[[2]], width = 14, height = 14)







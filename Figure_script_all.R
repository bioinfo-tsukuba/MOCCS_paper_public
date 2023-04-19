##################
### Figure1 ###
#################
options(timeout=1000)

## Fig. 1B --------
source("~/MOCCS_paper_public/function/Fig1B_plot_ver2.R")
Fig1B <- Fig1B_plot_ver2()
plot(Fig1B[[1]])
plot(Fig1B[[2]])
ggsave("~/MOCCS_paper_public/plot/Fig1/Fig1B_TF_ver2.pdf", Fig1B[[1]])
ggsave("~/MOCCS_paper_public/plot/Fig1/Fig1B_CTC_ver2.pdf", Fig1B[[2]])

## Fig. 1C --------
source("~/MOCCS_paper_public/function/Fig1C_plot.R")
target_ID_Fig1C <- "SRX1156473" #GATA3
# /Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1にhg38 target_ID_Fig1BのMOCCS outputを移してから以下を実行
Fig1C <- Fig1C_plot(target_ID_Fig1C)
plot(Fig1C)
ggsave(paste0("~/MOCCS_paper_public/plot/Fig1/Fig1C.pdf"), Fig1C, width = 7, height = 7)


F## Fig. 1D ------------
source("~/MOCCS_paper_public/function/Fig1D_plot_v2.R")
TF_list <- c("FOXA1", "SPI1", "CTCF")
load <- "local" #local or figshare
filter <- "hard" #soft or hard or all
p <- Fig1D_plot(TF_list, load, filter)
ggsave("~/MOCCS_paper_public/plot/Fig1/AUC_nega_3TF.pdf", p, width = 5, height = 7)


## Fig.1E------
source("~/MOCCS_paper_public/function/Fig1E_plot_v2.R")
TF_list <- c("FOXA1", "SPI1", "CTCF")
annotation_path <- "~/MOCCS_paper_public/data/Fig1/"
Fig1E <- Fig1E_plot(TF_list, annotation_path)
plot(Fig1E[[1]])
ggsave(paste0("~/MOCCS_paper_public/plot/Fig1/Fig1E_ver2.pdf"), Fig1E[[1]], width =7, height = 7)

# all TF p plot
source("~/MOCCS_paper_public/function/Fig1E_plot_v2.R")
TF_list_all <- readRDS("~/MOCCS_paper_public/data/Fig1/TF_all.rds")
annotation_path <- "~/MOCCS_paper_public/data/Fig1/"
Fig1E <- Fig1E_plot(TF_list_all, annotation_path)
plot(Fig1E[[2]])
ggsave(paste0("~/MOCCS_paper_public/plot/Fig1/Fig1E_ver2_allTF.pdf"), Fig1E[[2]], width =7, height = 7)


## Fig. 1F ------
source("~/MOCCS_paper_public/function/exec_Fig1_large_heatmap.R")
exec_Large_Heatmap()


##################
### Figure1G-2 ###
#################
source("~/MOCCS_paper_public/function/Fig1G_2_plot_v2.R")
Fig1G_2_plot(calc_opt = FALSE)

##################
### Figure3 ###
#################
## Fig. 3B --------
source("~/MOCCS_paper_public/function/Fig3B_plot.R")
simu_name <- "a01-02_N12000_varW5_stranded_sameA_func5_m90l45"
Fig3B_plot <- Fig3B_plot(simu_name)
ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3B.pdf"), Fig3B_plot, width = 7, height = 7)


## Fig. 3C --------
library(patchwork)
source("~/MOCCS_paper_public/function/Fig3C_plot.R")
simulation_path <- "~/MOCCS_paper_public/data/Fig4/"
Fig3C <- Fig3C_plot(simulation_path)
Fig3C_plot <- Fig3C[[1]] / Fig3C[[2]] / Fig3C[[3]]
ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3C_all.pdf"), Fig3C_plot, width = 14, height = 14)
ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3C_sensi.pdf"), Fig3C[[1]], width = 14, height = 14)
ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3C_spe.pdf"), Fig3C[[2]], width = 14, height = 14)
ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3C_FDR.pdf"), Fig3C[[3]], width = 14, height = 14)

## Fig. 3D --------
source("~/MOCCS_paper_public/function/Fig3D_plot.R")
target_ID1 <- "SRX150600" # JUN K562
target_ID2 <- "SRX186614" # JUN K562

#target_ID1 <- "SRX150600" # JUN K562
#target_ID2 <- "SRX150358" # JUN HUVEC

target_TF <- "JUN"

Fig3D <- Fig3D_plot(target_ID1, target_ID2, target_TF)
plot(Fig3D)
ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3D_", target_TF, "_", target_ID1, "_", target_ID2, ".pdf"), Fig3D, width = 7, height = 7)

## Fig. 3E --------
source("~/MOCCS_paper_public/function/Fig3E_plot.R")
target_ID1 <- "SRX150546" #JUN
target_ID2 <- "SRX029087" #FOS

#target_ID1 <- "SRX150546" #JUN
#target_ID2 <- "SRX190276" #CTCF

target_CT <- "K-562"
Fig3E <- Fig3E_plot(target_ID1, target_ID2, target_CT)
plot(Fig3E)
ggsave(paste0("~/MOCCS_paper_public/plot/Fig3/Fig3E_", target_CT, "_", target_ID1, "_", target_ID2, ".pdf"), Fig3E, width = 7, height = 7)




##################
### Figure4 ###
#################
## Fig. 4C --------
source("~/MOCCS_paper_public/function/Fig4C_plot.R")
target_TF <- "GATA3"
#annotation_path <- "~/MOCCS_paper_public/data/Fig5/"
Fig4C <- Fig4C_plot(target_TF)
plot(Fig4C)
ggsave(paste0("~/MOCCS_paper_public/plot/Fig4/Fig4C_", target_TF,".pdf"), Fig4C, width = 7, height = 7, dpi = 300)


## Fig. 4E --------
### scatter plot and bar plot
source("~/MOCCS_paper_public/function/Fig4E_plot.R")
path <- "~/MOCCS_paper_public/data/Fig5/"
target_TF <- "FOXA1"
Fig4E <- Fig4E_plot(target_TF, path)
plot(Fig4E[[1]])
plot(Fig4E[[2]])
ggsave(paste0("~/MOCCS_paper_public/plot/Fig4/Fig4E_scatter_", target_TF,".pdf"), Fig4E[[1]], width = 7, height = 7, dpi = 300)
ggsave(paste0("~/MOCCS_paper_public/plot/Fig4/Fig4E_bar_plot.pdf"), Fig4E[[2]], width = 7, height = 7)



##################
### Figure5 ###
#################
# Fig5BC, label---------
# devide by dMOCCS2score>0 and dMOCCS2score < 0 (label)
source("~/MOCCS_paper_public/function/Fig5BC_seihu_label.R")
figure <- "C" #C or D 
labels <- make_label(figure)

# dMOCCS2score > 0
plot(labels[[1]])
plot(labels[[2]])
plot(labels[[3]])
plot(labels[[4]])

# dMOCCS2score < 0
plot(labels[[5]])
plot(labels[[6]])
plot(labels[[7]])
plot(labels[[8]])

if(figure == "C"){
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_ID_plus.pdf", labels[[1]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_rs_plus.pdf", labels[[2]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_ref_plus.pdf", labels[[3]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_alt_plus.pdf", labels[[4]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_ID_minus.pdf", labels[[5]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_rs_minus.pdf", labels[[6]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_ref_minus.pdf", labels[[7]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_alt_minus.pdf", labels[[8]])
}else if(figure == "D"){
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_ID_plus.pdf", labels[[1]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_rs_plus.pdf", labels[[2]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_ref_plus.pdf", labels[[3]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_alt_plus.pdf", labels[[4]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_ID_minus.pdf", labels[[5]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_rs_minus.pdf", labels[[6]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_ref_minus.pdf", labels[[7]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_alt_minus.pdf", labels[[8]])
}


## Fig. 5B --------
source("~/MOCCS_paper_public/function/Fig5B_plot.R")
target_phenotype <- "SLE"
threshold <- 150
Fig5B <- Fig5B_plot(target_phenotype, threshold)
plot(Fig5B[[1]])
plot(Fig5B[[2]])
ggsave(paste0("~/MOCCS_paper_public/plot/Fig5/Fig5B_", target_phenotype,"_CL.pdf"), Fig5B[[1]], width = 7, height = 7)
#ggsave(paste0("~/MOCCS_paper_public/plot/Fig5/Fig5B_", target_phenotype,"_TF.pdf"), Fig5B[[2]], width = 7, height = 7)

## Fig. 5C --------
source("~/MOCCS_paper_public/function/Fig5C_plot.R")
target_phenotype <- "CD"
threshold <- 150
Fig5C <- Fig5C_plot(target_phenotype, threshold)
plot(Fig5C[[1]])
plot(Fig5C[[2]])
#ggsave(paste0("~/MOCCS_paper_public/plot/Fig6/Fig6D_", target_phenotype,"_CL.pdf"), Fig6D[[1]], width = 7, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/Fig5/Fig5C_", target_phenotype,"_TF.pdf"), Fig5C[[2]], width = 7, height = 7)


## Fig5BC devide by dMOCCS2score > 0 or dMOCCS2score < 0 -----
source("~/MOCCS_paper_public/function/Fig5BC_plot_seihu.R")
target_phenotype <- "SLE"
threshold <- 150
Fig5BC <- Fig5BC_plot_seihu(target_phenotype, threshold)
plot(Fig5BC[[1]])
plot(Fig5BC[[2]])
plot(Fig5BC[[3]])
plot(Fig5BC[[4]])

# Fig5BC label ---
source("~/MOCCS_paper_public/function/Fig5BC_seihu_label.R")
figure <- "D" #C or D 
labels <- make_label(figure)
plot(labels[[1]])
plot(labels[[2]])
plot(labels[[3]])
plot(labels[[4]])
if(figure == "C"){
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_ID.pdf", labels[[1]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_rs.pdf", labels[[2]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_ref.pdf", labels[[3]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/SLE_alt.pdf", labels[[4]])
}else if(figure == "D"){
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_ID.pdf", labels[[1]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_rs.pdf", labels[[2]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_ref.pdf", labels[[3]])
  ggsave("~/MOCCS_paper_public/plot/Fig5/CD_SLE_label/CD_alt.pdf", labels[[4]])
}

## Fig. 5D --------
source("~/MOCCS_paper_public/function/Fig5D_plot.R")
target_phenotype <- "CD"
target_rs <- "rs17293632"
target_position <- "chr15_67150258"
threshold <- 100
Fig5D <- Fig5D_plot(target_phenotype, target_rs, target_position ,threshold)
plot(Fig5D[[1]])
plot(Fig5D[[2]])
ggsave(paste0("~/MOCCS_paper_public/plot/Fig5/Fig5D_", target_phenotype,"_CL.pdf"), Fig5D[[1]], width = 7, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/Fig5/Fig5D_", target_phenotype,"_TF.pdf"), Fig5D[[2]], width = 7, height = 7)

### IBDの場合
# source("~/MOCCS_paper_public/function/Fig5D_plot_2.R")
# target_phenotype <- "IBD"
# target_rs <- "rs17293632"
# target_position <- "chr15_67150258"
# threshold <- 100
# Fig5D <- Fig5D_plot_2(target_phenotype, target_rs, target_position , threshold)
# plot(Fig5D[[1]])
# plot(Fig5D[[2]])
# ggsave(paste0("~/MOCCS_paper_public/plot/Fig5/Fig5D_", target_phenotype,"_CL.pdf"), Fig5D[[1]], width = 14, height = 14)
# ggsave(paste0("~/MOCCS_paper_public/plot/Fig5/Fig5D_", target_phenotype,"_TF.pdf"), Fig5D[[2]], width = 14, height = 14)


# devide by dMOCCS2score>0 and dMOCCS2score < 0
source("~/MOCCS_paper_public/function/Fig5D_plot_seihu.R")
target_phenotype <- "CD"
target_rs <- "rs17293632"
target_position <- "chr15_67150258"
threshold <- 100
Fig5D <- Fig5D_plot(target_phenotype, target_rs, target_position , threshold)
plot(Fig5D[[1]])
plot(Fig5D[[2]])



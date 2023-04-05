##################
### Supplement ###
##################
## Fig. S1  --------
### Fig. S1B --------
source("~/MOCCS_paper_public/function/FigS1B.R")
target_TF <- "CTCF"
annotation_path <- "~/MOCCS_paper_public/data/Fig1/"
FigS1B <- FigS1B_plot(target_TF, annotation_path)
plot(FigS1B)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS1/FigS1B_", target_TF,".pdf"), FigS1B, width = 7, height = 7)


### Fig. S1C  --------
source("~/MOCCS_paper_public/function/FigS1_plot.R")
path <- "~/MOCCS_paper_public/data/supplement/"
FigS1 <- FigS1_plot(path)
FigS1
ggsave(paste0("~/MOCCS_paper_public/plot/FigS1/FigS1C.pdf"), FigS1, width = 14, height = 14)


### Fig. S1D --------
source("~/MOCCS_paper_public/function/FigS1D.R")
annotation_path <- "~/MOCCS_paper_public/data/Fig1/"
FigS1D <- FigS1D_plot(annotation_path)
plot(FigS1D[[1]])
plot(FigS1D[[2]])
ggsave(paste0("~/MOCCS_paper_public/plot/FigS1/FigS1D_CL.pdf"), FigS1D[[1]], width = 7, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS1/FigS1D_TF.pdf"), FigS1D[[2]], width = 7, height = 7)


## Fig. S2 (FigS1 in new version (2022/02/21)) --------
source("~/MOCCS_paper_public/function/FigS2_plot.R")
simulation_path <- "~/MOCCS_paper_public/data/supplement/"
FigS2 <- FigS2_plot(simulation_path)
FigS2_plot <- FigS2[[1]] / FigS2[[2]] /FigS2[[3]]
ggsave(paste0("~/MOCCS_paper_public/plot/FigS2/FigS2_all.pdf"), FigS2_plot, width = 14, height = 14)


## Fig. S3 --------
source("~/MOCCS_paper_public/function/FigS3_plot.R")
annotation_path <- "~/MOCCS_paper_public/data/Fig1/"
path <- "~/MOCCS_paper_public/data/supplement/"
FigS3 <- FigS3_plot(annotation_path, path)
plot(FigS3[[1]])
plot(FigS3[[2]])
plot(FigS3[[3]])
ggsave(paste0("~/MOCCS_paper_public/plot/FigS3/FigS3_sample_hist.pdf"), FigS3[[1]], width = 7, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS3/FigS3_q005_density.pdf"), FigS3[[2]], width =7, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS3/FigS3_q001_density.pdf"), FigS3[[3]], width = 7, height = 7)S


## Fig. S9 --------
library(patchwork)
source("~/MOCCS_paper_public/function/FigS9_plot.R")
path <- "~/MOCCS_paper_public/data/supplement/FigS9/"
target_TF <- "GATA3" #FOXA1 or GATA3
#annotation_path <- "~/MOCCS_paper_public/data/Fig5/"
FigS9 <- FigS9_plot(path, target_TF)
FigS9_plot <- (FigS9[[1]][[1]] + FigS9[[1]][[2]] + FigS9[[1]][[3]]) / (FigS9[[1]][[4]] + FigS9[[1]][[5]] + FigS9[[1]][[6]])
plot(FigS9_plot)
plot(FigS9[[2]])
plot(FigS9[[3]])
plot(FigS9[[4]])
ggsave(paste0("~/MOCCS_paper_public/plot/FigS9/FigS9_position_", target_TF ,".pdf"), FigS9_plot, width = 21, height = 14)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS9/FigS9_", target_TF ,".pdf"), FigS9[[2]], width = 7, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS9/FigS9_allTF_PBS_Spe.pdf"), FigS9[[3]], width = 7, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS9/FigS9_allTF_PBS_Pear.pdf"), FigS9[[4]], width = 7, height = 7)


## Fig. S10 --------
source("~/MOCCS_paper_public/function/FigS10_plot_v5.R")
FigS10 <-FigS10_plot_v5()
plot(FigS10[[1]])
plot(FigS10[[2]])
plot(FigS10[[3]])
plot(FigS10[[4]])


## Fig. S11 --------
source("~/MOCCS_paper_public/function/Fig6C_plot.R")
target_phenotype <- "SLE"
annotation_path <- "~/MOCCS_paper_public/data/Fig6/"
threshold <- 150
Fig6C <- Fig6C_plot(target_phenotype, annotation_path, threshold)
plot(Fig6C[[1]])
plot(Fig6C[[2]])
ggsave(paste0("~/MOCCS_paper_public/plot/FigS11/FigS11_", target_phenotype,"_CL.pdf"), Fig6C[[1]], width = 7, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS11/FigS11_", target_phenotype,"_TF.pdf"), Fig6C[[2]], width = 7, height = 7)


### IBD
source("~/MOCCS_paper_public/function/Fig6C_plot_2.R")
target_phenotype <- "IBD"
annotation_path <- "~/MOCCS_paper_public/data/Fig6/"
threshold <- 150
Fig6C <- Fig6C_plot(target_phenotype, annotation_path, threshold)
plot(Fig6C[[1]])
plot(Fig6C[[2]])
ggsave(paste0("~/MOCCS_paper_public/plot/FigS11/FigS11_", target_phenotype,"_CL.pdf"), Fig6C[[1]], width = 7, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS11/FigS11_", target_phenotype,"_TF.pdf"), Fig6C[[2]], width = 7, height = 7)


## Fig. S12 ---------
source("~/MOCCS_paper_public/function/FigS12_plot.R")
target_phenotype <- "IBD"
path <- "~/MOCCS_paper_public/data/supplement/FigS12/"
FigS12 <- FigS12_plot(target_phenotype, path)
plot(FigS12[[1]])
plot(FigS12[[2]])
ggsave(paste0("~/MOCCS_paper_public/plot/FigS12/FigS12A_", target_phenotype,".pdf"), FigS12[[1]], width = 7, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS12/FigS12B_", target_phenotype,".pdf"), FigS12[[2]], width = 7, height = 7)

## Fig. S13 ---------
source("~/MOCCS_paper_public/function/FigS13_plot.R")
TF_list <- c("CTCF", "SPI1", "FOXA1")
FigS13 <- FigS13_plot(TF_list)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS13/FigS13_", names(FigS13)[1], ".pdf"),  FigS13[[1]], width = 14, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS13/FigS13_", names(FigS13)[2], ".pdf"),  FigS13[[2]], width = 14, height = 7)
ggsave(paste0("~/MOCCS_paper_public/plot/FigS13/FigS13_", names(FigS13)[3], ".pdf"),  FigS13[[3]], width = 14, height = 7)


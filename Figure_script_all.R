##################
### Figure1 ###
#################
## Fig. 1B
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1B_plot.R")
target_ID_Fig1B <- "SRX160844"
# /Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1にhg38 target_ID_Fig1BのMOCCS outputを移してから以下を実行
Fig1B <- Fig1B_plot(target_ID_Fig1B)
plot(Fig1B)

## Fig. 1C
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1C_plot.R")
target_TF <- "JUN"
Fig1C_plot(target_TF)

## Fig. 1D
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1D_plot.R")
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/"
Fig1D <- Fig1D_plot(annotation_path)
plot(Fig1D[[1]])
plot(Fig1D[[2]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig1/Fig1D_CL.pdf"), Fig1D[[1]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig1/Fig1D_TF.pdf"), Fig1D[[2]])

## Fig. 1E
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1E_plot.R")
target_TF <- "CTCF"
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/"
Fig1E <- Fig1E_plot(target_TF, annotation_path)
plot(Fig1E)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig1/Fig1E_", target_TF,".pdf"), Fig1E)


## Fig. 1F
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1F_plot.R")
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/"
Fig1F <- Fig1F_plot(annotation_path)
plot(Fig1F)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/Fig1/Fig1F.pdf"), Fig1F)


## Fig. 1G
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1G_plot.R")
target_TF <- "CTCF"
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/"
Fig1G_plot(target_TF, annotation_path)


##################
### Figure2 ###
#################

##################
### Figure3 ###
#################

##################
### Figure4 ###
#################

##################
### Figure5 ###
#################

##################
### Figure6 ###
#################

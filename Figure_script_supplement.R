##################
### Supplement ###
#################
## Fig. S1 (FigS2 in new version (2022/02/21))
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/FigS1_plot.R")
simulation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/supplement/"
FigS1 <- FigS1_plot(simulation_path)
FigS1_plot <- FigS1[[1]] / FigS1[[2]] 
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS1_all.pdf"), FigS1_plot, width = 21, height = 14)


## Fig. S2 (FigS1 in new version (2022/02/21))
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/FigS2_plot.R")
path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/supplement/"
FigS2 <- FigS2_plot(path)
FigS2[[1]]
FigS2[[2]]
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS2_before_filter.pdf"), FigS2[[1]], width = 14, height = 14)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS2_after_hardfilter.pdf"), FigS2[[2]], width = 14, height = 14)


## Fig. S3
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/FigS3_plot.R")
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig1/"
path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/supplement/"
FigS3 <- FigS3_plot(annotation_path, path)
plot(FigS3[[1]])
plot(FigS3[[2]])
plot(FigS3[[3]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS3_sample_hist.pdf"), FigS3[[1]], width = 14, height = 14)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS3_q005_density.pdf"), FigS3[[2]], width = 14, height = 14)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS3_q001_density.pdf"), FigS3[[3]], width = 14, height = 14)


## Fig. S6
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/FigS6_plot.R")
path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/supplement/simulation/"
simu_kind <- "l100m50"
FigS6 <- FigS6_plot(path, simu_kind)
plot(FigS6)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS6_", simu_kind,".pdf"), FigS6, width = 14, height = 14)



## Fig. S7
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/FigS7_plot.R")
path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/supplement/FigS7/"
target_tf <- "FOXA1" #FOXA1 or GATA3
annotation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig5/"
FigS7 <- FigS7_plot(path, target_tf, annotation_path)
plot(FigS7[[1]])
plot(FigS7[[2]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS7_", target_tf ,".pdf"), FigS7[[1]], width = 7, height = 7)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS7_position_", target_tf ,".pdf"), FigS7[[2]], width = 42, height = 7)



## Fig. S8
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/FigS8_plot.R")
path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig6/"
target_phenotype <- "CD"
target_tf <- "GATA3"
FigS8 <- FigS8_plot(path, target_phenotype, target_tf)
plot(FigS8[[1]])
plot(FigS8[[2]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS8_", target_tf, ".pdf"), FigS8[[1]], width = 7, height = 7)
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS8_all_barplot.pdf"), FigS8[[2]], width = 14, height = 7)

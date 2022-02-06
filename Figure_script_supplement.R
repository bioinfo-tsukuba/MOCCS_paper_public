##################
### Supplement ###
#################
## Fig. S1
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/FigS1_plot.R")
simulation_path <- "/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/supplement/"
FigS1 <- FigS1_plot(simulation_path)
plot(FigS1[[1]])
plot(FigS1[[2]])
plot(FigS1[[3]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS1_alpha.pdf"), FigS1[[1]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS1_N.pdf"), FigS1[[2]])
ggsave(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/plot/supplement/FigS1_sigma.pdf"), FigS1[[3]])


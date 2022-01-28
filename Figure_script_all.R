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

## Fig. 1E
source("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/function/Fig1E_plot.R")
target_TF <- "CTCF"
filter_name <- c("all", "hard", "soft") ## "all", "hard", "soft"
filter_auc <- list()
auc_df_all <- c()
for (i in 1:length(filter_name)) {
  target_filter <- filter_name[i]
  Fig1E_plot(target_TF, filter)
  filter_auc[[target_filter]] <- Fig1E_plot(target_TF, filter)
  filter_auc_df <- tibble(filter = rep(target_filter, nrow(filter_auc[[target_filter]])), auc = filter_auc[[target_filter]])
  if(i == i){
    auc_df_all <- filter_auc_df
  }else{
    auc_df_all <- rbind(auc_df_all, filter_auc_df)
  }
}

auc_df_all %>% ggplot(aes(x = filter, y = auc)) +
  geom_boxplot() +
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold")
  )




## Fig. 1F

## Fig. 1G

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

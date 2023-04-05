Fig4E_plot <- function(target_TF, path){
  
  library(tidyverse)
  options(timeout=1000)
  url_list <- read_csv("~/MOCCS_paper_public/data/Fig5/Fig5E_url_list.csv", col_names = FALSE)
  colnames(url_list) <- c("TF", "url")
  target_url <- url_list %>% filter(TF == target_TF) %>% .$url %>% as.character()
  if(target_TF != "CTCF"){
    df_all <- readRDS(url(target_url, "rb"))
    df_all2 <- df_all %>% mutate(position = end - 6 + posi)
    df3 <- df_all2 %>% mutate(fdrp_bh_ref_log = -log10(fdrp_bh_ref),
                              fdrp_bh_alt_log = -log10(fdrp_bh_alt), 
                              ASB = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, fdrp_bh_alt_log, fdrp_bh_ref_log), 
                              x_axis = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, "positive", "negative"))
    df4 <- df3 %>% filter(x_axis == "positive") %>% mutate(ASB_plot = ASB)
    df5 <- df3 %>% filter(x_axis == "negative") %>% mutate(ASB_plot = -ASB)
    df6 <- rbind(df4, df5)
    threshold <- -log10(0.05)
    df11 <- df6 %>% mutate(color1 = ifelse((ASB_plot < -threshold  & q_value < 0.05 & dMOCCS2score > 0) | (ASB_plot > threshold  & q_value < 0.05 & dMOCCS2score < 0), "concordant", "disconcordant")) 
    df12 <- df11 %>% filter(color1 == "concordant") %>% mutate(color2 = color1) %>%
      select(-fdrp_bh_ref, -fdrp_bh_alt, -MOCCS2score_before, -MOCCS2score_after, -fdrp_bh_ref_log, -fdrp_bh_alt_log)
    df13 <- df11 %>%  filter(color1 == "disconcordant") %>% mutate(color2 = ifelse((q_value < 0.05 & abs(ASB_plot) > threshold), "disconcordant", "nondifferent")) %>%
      select(-fdrp_bh_ref, -fdrp_bh_alt, -MOCCS2score_before, -MOCCS2score_after, -fdrp_bh_ref_log, -fdrp_bh_alt_log)
    rm(df11)
    df14 <- rbind(df12, df13)
    rm(df12, df13)
    
    p1 <- df14 %>% ggplot(aes(x = ASB_plot, y = dMOCCS2score, color = color2)) +
      geom_point(size = 0.5, alpha = 0.7) +
      xlab("ASB significance") +
      ylab("dMOCCS2score") +
      scale_color_manual(
        values = c(
          concordant = "#DC143C",
          disconcordant = "#1E90FF", 
          nondifferent = "gray"
        )
      )+
      ggtitle(paste0(target_TF, " ", "dMOCCS2score q<0.05, threshold q<0.05"  ))+
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
    
    
  }else{
    #CTCF
    #df_all1 <- readRDS(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/dMOCCS2score_output_all/ADASTRA_dMOCCS2score_", target_TF, "_all_part1.rds")) 
    df_all1 <- readRDS(url("https://figshare.com/ndownloader/files/34065725", "rb"))
    df_all2 <- readRDS(url("https://figshare.com/ndownloader/files/34065728", "rb"))
    df_all3 <- readRDS(url("https://figshare.com/ndownloader/files/34065737", "rb"))
    df_all4 <- readRDS(url("https://figshare.com/ndownloader/files/34065740", "rb"))
    df_all <- df_all1 %>% add_row(df_all2, df_all3, df_all4)
    
    df_all2 <- df_all %>% mutate(position = end - 6 + posi)
    rm(df_all)
    df3 <- df_all2 %>% mutate(fdrp_bh_ref_log = -log10(fdrp_bh_ref),
                              fdrp_bh_alt_log = -log10(fdrp_bh_alt), 
                              ASB = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, fdrp_bh_alt_log, fdrp_bh_ref_log), 
                              x_axis = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, "positive", "negative"))
    rm(df_all2)
    df4 <- df3 %>% filter(x_axis == "positive") %>% mutate(ASB_plot = ASB)
    df5 <- df3 %>% filter(x_axis == "negative") %>% mutate(ASB_plot = -ASB)
    rm(df3)
    df6 <- rbind(df4, df5)
    rm(df4, df5)
    threshold <- -log10(0.05)
    df11 <- df6 %>% mutate(color1 = ifelse((ASB_plot < -threshold  & q_value < 0.05 & dMOCCS2score > 0) | (ASB_plot > threshold  & q_value < 0.05 & dMOCCS2score < 0), "concordant", "disconcordant")) 
    rm(df6)
    df12 <- df11 %>% filter(color1 == "concordant") %>% mutate(color2 = color1) %>%
      select(-fdrp_bh_ref, -fdrp_bh_alt, -MOCCS2score_before, -MOCCS2score_after, -fdrp_bh_ref_log, -fdrp_bh_alt_log)
    df13 <- df11 %>%  filter(color1 == "disconcordant") %>% mutate(color2 = ifelse((q_value < 0.05 & abs(ASB_plot) > threshold), "disconcordant", "nondifferent")) %>%
      select(-fdrp_bh_ref, -fdrp_bh_alt, -MOCCS2score_before, -MOCCS2score_after, -fdrp_bh_ref_log, -fdrp_bh_alt_log)
    rm(df11)
    df14 <- rbind(df12, df13)
    rm(df12, df13)
    
    p1 <- df14 %>% ggplot(aes(x = ASB_plot, y = dMOCCS2score, color = color2)) +
      geom_point(size = 0.5, alpha = 0.7) +
      xlab("ASB significance") +
      ylab("dMOCCS2score") +
      scale_color_manual(
        values = c(
          concordant = "#DC143C",
          disconcordant = "#1E90FF", 
          nondifferent = "gray"
        )
      )+
      ggtitle(paste0(target_TF, " ", "dMOCCS2score q<0.05, threshold q<0.05"  ))+
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
    
    
  }
  
  # bar plot
    
    bar_df <- read_csv("~/MOCCS_paper_public/data/Fig5/Fig5E_barplot.csv")
    bar_df2 <- bar_df %>% pivot_longer(-c(Antigen, number_of_concordant_SNPs_q005, number_of_disconcordant_SNPs_q005, 
                                          nega_number_of_concordant_SNPs_q005, nega_number_of_disconcordant_SNPs_q005, 
                                          nega_concordant_ratio_q005, nega_disconcordant_ratio_q005 ,
                                          q_value),
                                       names_to = "concordant_or_disconcordant", values_to = "ratio")
    bar_df3 <- bar_df2 %>% mutate(color = gsub("_ratio_q005", "", bar_df2$concordant_or_disconcordant))
    
    bar_df3$color2 <- factor(bar_df3$color, levels = c("disconcordant", "concordant"))
    tf_factor_list <- bar_df3 %>% filter(color2 == "concordant") %>% arrange(desc(ratio)) %>% .$Antigen %>% as.character() %>% unique()
    bar_df3$Antigen <- factor(bar_df3$Antigen, levels = tf_factor_list)
    p2 <- bar_df3 %>% ggplot(aes(x = Antigen, y = ratio, fill = color2)) +
      geom_bar(stat = "identity", position = "fill") +
      #geom_text(aes(label = ratio), size = 2, hjust = 0.5, vjust = 3, position = "stack") +
      scale_fill_manual(values = c("#1E90FF", "#DC143C")) +
      theme(plot.title = element_text(face="bold",hjust = 0.5), 
            panel.grid.major = element_line(colour = "gray"),
            panel.grid.minor = element_line(colour="gray"),
            panel.background = element_blank(), 
            axis.line = element_line(colour="black"),
            axis.text=element_text(size=12,face="bold"),
            axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
            axis.text.y =element_text(size=10,face="bold"),
            axis.title=element_text(size=14,face="bold"),
            aspect.ratio = 1
      )+
      labs(fill = "")
  
  
  return(list(p1, p2))
  
}
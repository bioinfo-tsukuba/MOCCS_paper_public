FigS6_plot <- function(path, simu_kind){
  
  library(tidyverse)
  simu_result <- readRDS(paste0(path, "result_a01-02_N12000_varW5_stranded_sameA_func5_", simu_kind, ".rds"))
  simu_name <- paste0("a01-02_N12000_varW5_stranded_sameA_func5_", simu_kind)
  
  # remove 1bp shifted k-mer
  B1_and_B2 <- simu_result %>% filter(true_kmer_anotation == "B1" | true_kmer_anotation == "B2" | true_kmer_anotation == "A") %>% select(kmer, simu_num, true_kmer_anotation)
  B1_and_B2 %>% group_by(simu_num, true_kmer_anotation) %>% summarise(n=n())
  
  B1_and_B2_shifted <- list()
  for (j in 1:50) {
    target_B1_and_B2 <- B1_and_B2 %>% filter(simu_num == j) %>%.$kmer %>% as.character() %>% unique()
    target_B1_and_B2_shifted <- c()
    for (i in 1:length(target_B1_and_B2)) {
      target_kmer <- target_B1_and_B2[i]
      tmp <- substring(target_kmer, 2, 6)
      target_kmer_shifted <- c()
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp, "A"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp, "T"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp, "G"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp, "C"))
      
      tmp2 <- substring(target_kmer, 1, 5)
      target_kmer_shifted <- c(target_kmer_shifted,paste0("A", tmp2))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("T", tmp2))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("G", tmp2))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("C", tmp2))
      target_kmer_shifted <- unique(target_kmer_shifted)
      
      target_B1_and_B2_shifted <- c(target_B1_and_B2_shifted, target_kmer_shifted) 
      B1_and_B2_shifted[j] <- list(unique(target_B1_and_B2_shifted))
    }
  }
  
  for (j in 1:50) {
    target_rm_kmer <- B1_and_B2_shifted[[j]]
    target_df <-simu_result %>% filter(simu_num == j) %>% filter((true_kmer_anotation == "C") & (!kmer %in% target_rm_kmer))
    if(j == 1){
      simu_result_new <- target_df 
    }else{
      simu_result_new <- rbind(simu_result_new, target_df)
    }
  }
  
  only_A_B <- simu_result %>% filter(true_kmer_anotation == "A" | true_kmer_anotation == "B1" | true_kmer_anotation == "B2")
  simu_result_new <- rbind(simu_result_new, only_A_B)
  
  
  # remove 2bp shifted k-mer
  B1_and_B2_shifted <- list()
  for (j in 1:50) {
    target_B1_and_B2 <- B1_and_B2 %>% filter(simu_num == j) %>%.$kmer %>% as.character() %>% unique()
    target_B1_and_B2_shifted <- c()
    for (i in 1:length(target_B1_and_B2)) {
      target_kmer <- target_B1_and_B2[i]
      
      tmp <- substring(target_kmer, 2, 6)
      target_kmer_shifted <- c()
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp, "A"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp, "T"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp, "G"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp, "C"))
      
      tmp2 <- substring(target_kmer, 1, 5)
      target_kmer_shifted <- c(target_kmer_shifted,paste0("A", tmp2))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("T", tmp2))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("G", tmp2))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("C", tmp2))
      
      tmp3 <- substring(target_kmer, 3, 6)
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "AA"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "AT"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "AG"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "AC"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "TA"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "TT"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "TG"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "TC"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "GA"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "GT"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "GG"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "GC"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "CA"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "CT"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "CG"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0(tmp3, "CC"))
      
      tmp4 <- substring(target_kmer, 1, 4)
      target_kmer_shifted <- c(target_kmer_shifted,paste0("AA", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("AT", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("AG", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("AC", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("TA", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("TT", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("TG", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("TC", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("GA", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("GT", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("GG", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("GC", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("CA", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("CT", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("CG", tmp4))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("CC", tmp4))
      
      tmp5 <- substring(target_kmer, 2, 5)
      target_kmer_shifted <- c(target_kmer_shifted,paste0("A", tmp5, "A"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("A", tmp5, "T"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("A", tmp5, "G"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("A", tmp5, "C"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("T", tmp5, "A"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("T", tmp5, "T"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("T", tmp5, "G"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("T", tmp5, "C"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("G", tmp5, "A"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("G", tmp5, "T"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("G", tmp5, "G"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("G", tmp5, "C"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("C", tmp5, "A"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("C", tmp5, "T"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("C", tmp5, "G"))
      target_kmer_shifted <- c(target_kmer_shifted,paste0("C", tmp5, "C"))
      
      target_kmer_shifted <- unique(target_kmer_shifted)
      
      target_B1_and_B2_shifted <- c(target_B1_and_B2_shifted, target_kmer_shifted) 
      B1_and_B2_shifted[j] <- list(unique(target_B1_and_B2_shifted))
    }
  }
  
  for (j in 1:50) {
    target_rm_kmer <- B1_and_B2_shifted[[j]]
    target_df <- simu_result %>% filter(simu_num == j) %>% filter((true_kmer_anotation == "C") & (!kmer %in% target_rm_kmer))
    if(j == 1){
      simu_result_new2 <- target_df 
    }else{
      simu_result_new2 <- rbind(simu_result_new2, target_df)
    }
  }
  
  simu_result_new2 <- rbind(simu_result_new2, only_A_B)
  
  p <- simu_result_new2 %>% mutate(color = ifelse(q_value < 0.05, "differential", "nondifferential")) %>%
    ggplot(aes(x = MOCCS2score, y = MOCCS2score2, color = color)) +
    geom_point(size = 1) +
    scale_colour_manual(
      values = c(
        differential = "red",
        nondifferential = "gray"
      )
    )+
    theme(plot.title = element_text(face="bold",hjust = 0.5), 
          panel.grid.major = element_line(colour = "gray"),
          panel.grid.minor = element_line(colour="gray"),
          panel.background = element_blank(), 
          axis.line = element_line(colour="black"),
          axis.text=element_text(size=12,face="bold"),
          axis.text.x =element_text(size=10,face="bold", angle = 45, hjust = 1),
          axis.text.y =element_text(size=10,face="bold"),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_blank()
    ) +
    xlab("MOCCS2score1") +
    ggtitle(paste0(simu_name, "_q005"))
  
  return(p)
  
  
}
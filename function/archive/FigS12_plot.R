FigS12_plot <- function(target_phenotype, path){
  library(tidyverse)
  if(target_phenotype == "CD"){
    df_joined <- readRDS(url("https://figshare.com/ndownloader/files/34669810", "rb"))
  }else if(target_phenotype == "MS"){
    df_joined <- readRDS(url("https://figshare.com/ndownloader/files/34669816", "rb"))
  }else if(target_phenotype == "SLE"){
    df_joined <- readRDS(url("https://figshare.com/ndownloader/files/34669831", "rb"))
  }else if(target_phenotype == "IBD"){
    df_joined <- readRDS(url("https://figshare.com/ndownloader/files/34669801", "rb"))
  }
  df_joined %>%
    mutate(
      category_allele_frequency = cut_interval(risk_allele_frequency, n=5),
      diff_k_mer = (q_value < 0.05)
    ) -> df_tmp
  
  plot1 <- df_tmp %>%
    mutate(dMOCCS2score_ = case_when(
      q_value<0.05 ~ dMOCCS2score,
      TRUE ~ NaN
    )) %>%
    group_by(snp_for_join, category_allele_frequency) %>%
    summarise(max_asb_dMOCCS2score = max(abs(dMOCCS2score_),na.rm = TRUE)) %>%
    ggplot(aes(category_allele_frequency, max_asb_dMOCCS2score)) +
    geom_boxplot() +
    theme_bw()   +
    ggtitle(target_phenotype)+
    xlab("category of allele frequency") +
    ylab("absolute values of dMOCCS2score")
  
  plot2 <- df_tmp %>%
    group_by(snp_for_join, category_allele_frequency) %>%
    summarise(any_diff_kmer = sum((q_value < 0.05),na.rm = TRUE) >0 ) %>%
    group_by(category_allele_frequency) %>%
    summarise(fraction_any_diff_kmer = sum(any_diff_kmer)/n()) %>%
    ggplot(aes(category_allele_frequency, fraction_any_diff_kmer)) +
    geom_point(size = 3) +
    theme_bw()  +
    ylim(c(0,1)) +
    ggtitle(target_phenotype)+
    xlab("category of allele frequency") +
    ylab("ratio of SNPs with significant dMOCCS2score")
  
  return(list(plot1, plot2))
}
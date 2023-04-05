Fig1C_plot <- function(target_ID_Fig1C){
  
  library(tidyverse)
  MOCCS_output_path <- paste0("~/MOCCS_paper_public/data/Fig1/", target_ID_Fig1B, "_6mer_v2.auc_count.txt")
  MOCCS_output_target <- read_tsv(MOCCS_output_path)
  
  ## urlで読み込む場合
  #MOCCS_output_path <- "https://raw.githubusercontent.com/bioinfo-tsukuba/MOCCS-DB_paper/main/data/Fig1/SRX1156473_6mer_v2.auc_count.txt?token=GHSAT0AAAAAABN3UUI2U3BJTJOIVQF3JIBMYQLF4QA"
  #githubならrowのurlを使う
  #MOCCS_output_target <- read_tsv(url("https://figshare.com/ndownloader/files/34066733", "rb"))
  #figshareならdownloadのurlを使う
  #MOCCS_output_target <- read_tsv(url(MOCCS_output_path, "rb"))
  
  
  # calculate pvalue
  W <- 350
  if(nrow(MOCCS_output_target) != 0){
    # kmerごと(行ごと)に p valueを計算する
    p_list <- lapply(1:nrow(MOCCS_output_target), function(y){
      target_AUC <- as.numeric(MOCCS_output_target[y,2])
      target_kmer_count <- as.numeric(MOCCS_output_target[y,3])
      target_p <- 1-pnorm(target_AUC, mean = 0, sd = sqrt(W^2/12/target_kmer_count))
      return(target_p)
    })
        
    # sampleごとに多重検定補正
    p_list <- unlist(p_list)
    q_list <- p.adjust(p_list)
        
    # IDとpvalueとqvalueを足したtableにして返す 
    MOCCS_output_target_2 <- MOCCS_output_target %>% mutate(p_value = p_list, q_value = q_list, ID = rep(target_ID_Fig1B, nrow(MOCCS_output_target)))
  } #if
  saveRDS(MOCCS_output_target_2, paste0("~/MOCCS_paper_public/data/Fig1/", target_ID_Fig1B, "_MOCCSout_qval.rds"))
  
  
  # plot
  kmer_all <- MOCCS_output_target_2 %>% arrange(desc(MOCCS2score)) %>% .$kmer
  selected_kmer <- kmer_all[1]
  kmer_all[!kmer_all %in% selected_kmer] <- NA
  MOCCS_output_target_2 <- MOCCS_output_target_2 %>% arrange(desc(MOCCS2score)) %>% mutate(kmer_label = kmer_all) %>% arrange(desc(MOCCS2score))
  
  p <- MOCCS_output_target_2 %>% ggplot(aes(x = reorder(kmer, desc(MOCCS2score)), y = MOCCS2score)) +
    geom_col() +
    geom_text(aes(y = Inf,label=kmer_label),size=9,hjust=0, vjust = 1) +
    xlab("k-mer") +
    ylab("MOCCS2score") +
    ggtitle(paste0(target_ID_Fig1B ," example of MOCCS2score distribution"))+
    theme(axis.text.x = element_blank(),
          axis.line.x.bottom  = element_blank(),
          plot.title = element_text(hjust = 0.5),
          title=element_text(size=12,face="bold"),
          aspect.ratio = 1) +
    theme(aspect.ratio = 1)
  return(p)
  
}
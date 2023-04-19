make_label <- function(figure){
    
    
    if(figure == "C"){
      target_phenotype <- "SLE"
      GWAS_df <- suppressMessages(read_tsv(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig6/GWAS_catalog_phenotype/", target_phenotype, ".tsv")))
      GWAS_df2 <- GWAS_df  %>% drop_na(`SNPS`, `STRONGEST SNP-RISK ALLELE`) %>% 
        select(`CHR_ID`, `CHR_POS`, `SNPS`, `STRONGEST SNP-RISK ALLELE`) %>% 
        unite("position", c(`CHR_ID`, `CHR_POS`), sep = "_") %>%
        separate(`STRONGEST SNP-RISK ALLELE`, into = c("rs", "alt_bp"), sep = "-") %>%
        select(-rs) %>%
        unite("snp_for_join", c(position, alt_bp))
      GWAS_df2$snp_for_join <- paste0(rep("chr", nrow(GWAS_df2)), GWAS_df2$snp_for_join)
      GWAS_df2 <- distinct(GWAS_df2)
      colnames(GWAS_df2) <- c("snp_for_join", "SNP_rs")
      
      # join rs number -----
      annotation_path <- "/Users/saeko/MOCCS_paper_public/data/Fig6/"
      threshold <- 150
      target_df <- readRDS(paste0(annotation_path, target_phenotype, "_peak_rand_binded_all_qval.rds"))
      
      totalization <- readRDS("~/MOCCS_paper_public/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
      annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
      df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID") %>% filter(q_value < 0.05)
      df1 <- df_phenotype_binded_all_selected_annotated %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% select(ID, p_value, Cell_type_class, Cell_type, Antigen, abs_dMOCCS2score) %>% distinct()
      
      snp_for_join <- df_phenotype_binded_all_selected_annotated %>% unite("snp_for_join", c(position, alt_bp)) %>% .$snp_for_join %>% as.character()
      df_phenotype_binded_all_selected_annotated2 <- df_phenotype_binded_all_selected_annotated %>% select(-snp_for_join) %>% mutate(snp_for_join = snp_for_join)
      df2 <- df_phenotype_binded_all_selected_annotated2 %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% group_by(ID, snp_for_join, Antigen, Cell_type_class, Cell_type) %>% summarise(max_abs_dMOCCS2score = max(abs_dMOCCS2score))
      df3 <- df2 %>% left_join(GWAS_df2, by = "snp_for_join") 
      
      # ID -----
      x_lab_ID <- df3  %>% unite("ID_snp_for_join", c(ID, snp_for_join)) %>% filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% select(ID_snp_for_join)
      x_lab_ID2 <- x_lab_ID %>% separate(ID_snp_for_join, into = c("ID", "chr", "position", "alt_bp"), sep = "_") %>% .$ID %>% as.character() %>% rev()
      p1 <- df3  %>% unite("ID_snp_for_join", c(ID, snp_for_join)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
        ggplot(aes(x = reorder(ID_snp_for_join, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Cell_type_class))+
        geom_col() +
        theme(plot.title = element_text(face="bold",hjust = 0.5), 
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=12,face="bold"),
              axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
              axis.text.y =element_text(size=10,face="bold"),
              axis.title=element_text(size=14,face="bold"),
              #legend.position = 'none',
              legend.title = element_blank()
        )+
        scale_x_discrete("ID", labels = x_lab_ID2) +
        ylab("max(|dMOCCS2score|)") +
        ggtitle(target_phenotype) +
        coord_flip()
      
      # SNP_rs label ----- ##ここから！！
      #x_lab_ID <- df3  %>% #unite("ID_snp_for_join", c(ID, snp_for_join)) %>% 
        #filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% select(SNP_rs)
      #saveRDS(x_lab_ID, paste0("/Users/saeko/MOCCS_paper_public/data/Fig6/label/", target_phenotype, "_label_TF.rds"))
      x_lab_ID <- readRDS("/Users/saeko/MOCCS_paper_public/data/Fig6/label/SLE_label_TF_rs_added.rds")
      x_lab_ID2 <- x_lab_ID %>% .$SNP_rs %>% as.character() %>% rev()
      
      
      p2 <- df3  %>% unite("ID_snp_for_join", c(ID, snp_for_join)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
        ggplot(aes(x = reorder(ID_snp_for_join, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Cell_type_class))+
        geom_col() +
        theme(plot.title = element_text(face="bold",hjust = 0.5), 
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=12,face="bold"),
              axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
              axis.text.y =element_text(size=10,face="bold"),
              axis.title=element_text(size=14,face="bold"),
              #legend.position = 'none',
              legend.title = element_blank()
        )+
        scale_x_discrete("", labels = x_lab_ID2) +
        ylab("") +
        ggtitle(target_phenotype) +
        coord_flip()
      
      # 6C ref k-mer label ----
      df4 <- df_phenotype_binded_all_selected_annotated2 %>% 
        mutate(abs_dMOCCS2score = abs(dMOCCS2score))%>% 
        select(ID, snp_for_join, Antigen, Cell_type_class, Cell_type, kmer_before_tou, kmer_after_tou, 
               dMOCCS2score, q_value, abs_dMOCCS2score) %>%
        group_by(ID, snp_for_join, Antigen, Cell_type_class, Cell_type) %>% 
        filter(abs_dMOCCS2score == max(abs_dMOCCS2score)) %>% distinct()
      colnames(df4) <-   c("ID", "snp_for_join", "Antigen", "Cell_type_class", "Cell_type", "ref_kmer", "alt_kmer", "dMOCCS2score",
                           "q_value", "max_abs_dMOCCS2score")
      
      x_lab_ID <- df4  %>% 
        filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% 
        select(ref_kmer) 
      x_lab_ID2 <- x_lab_ID %>% .$ref_kmer %>% as.character() %>% rev()
      p3 <- df4  %>% unite("ID_snp_for_join", c(ID, snp_for_join)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
        ggplot(aes(x = reorder(ID_snp_for_join, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Cell_type_class))+
        geom_col() +
        theme(plot.title = element_text(face="bold",hjust = 0.5), 
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=12,face="bold"),
              axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
              axis.text.y =element_text(size=10,face="bold"),
              axis.title=element_text(size=14,face="bold"),
              #legend.position = 'none',
              legend.title = element_blank()
        )+
        scale_x_discrete("", labels = x_lab_ID2) +
        ylab("") +
        ggtitle(target_phenotype) +
        coord_flip()
      
      # 6C alt k-mer label ----
      x_lab_ID <- df4  %>% 
        filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% 
        select(alt_kmer) 
      x_lab_ID2 <- x_lab_ID %>% .$alt_kmer %>% as.character() %>% rev()
      p4 <- df4  %>% unite("ID_snp_for_join", c(ID, snp_for_join)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
        ggplot(aes(x = reorder(ID_snp_for_join, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Cell_type_class))+
        geom_col() +
        theme(plot.title = element_text(face="bold",hjust = 0.5), 
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=12,face="bold"),
              axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
              axis.text.y =element_text(size=10,face="bold"),
              axis.title=element_text(size=14,face="bold"),
              #legend.position = 'none',
              legend.title = element_blank()
        )+
        scale_x_discrete("", labels = x_lab_ID2) +
        ylab("") +
        ggtitle(target_phenotype) +
        coord_flip()
      
      
      
    }else if(figure == "D"){
      
      target_phenotype <- "CD"
      GWAS_df <- read_tsv(paste0("/Users/saeko/Documents/MOCCS/paper_figure/MOCCS-DB_paper/data/Fig6/GWAS_catalog_phenotype/", target_phenotype, ".tsv"))
      GWAS_df2 <- GWAS_df  %>% drop_na(`SNPS`, `STRONGEST SNP-RISK ALLELE`) %>% 
        select(`CHR_ID`, `CHR_POS`, `SNPS`, `STRONGEST SNP-RISK ALLELE`) %>% 
        unite("position", c(`CHR_ID`, `CHR_POS`), sep = "_") %>%
        separate(`STRONGEST SNP-RISK ALLELE`, into = c("rs", "alt_bp"), sep = "-") %>%
        select(-rs) %>%
        unite("snp_for_join", c(position, alt_bp))
      GWAS_df2$snp_for_join <- paste0(rep("chr", nrow(GWAS_df2)), GWAS_df2$snp_for_join)
      GWAS_df2 <- distinct(GWAS_df2)
      colnames(GWAS_df2) <- c("snp_for_join", "SNP_rs")
      
      # join rs number -----
      annotation_path <- "/Users/saeko/MOCCS_paper_public/data/Fig6/"
      threshold <- 150
      target_df <- readRDS(paste0(annotation_path, target_phenotype, "_peak_rand_binded_all_qval.rds"))
      
      totalization <- readRDS("~/MOCCS_paper_public/data/Fig1/MOCCSout_hg38_all_qval_annotated.rds")
      annotation <- totalization %>% select(ID, Antigen, Cell_type_class, Cell_type) %>% distinct()
      df_phenotype_binded_all_selected_annotated <- target_df %>% left_join(annotation, by = "ID") %>% filter(q_value < 0.05)
      df1 <- df_phenotype_binded_all_selected_annotated %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% select(ID, p_value, Cell_type_class, Cell_type, Antigen, abs_dMOCCS2score) %>% distinct()
      
      snp_for_join <- df_phenotype_binded_all_selected_annotated %>% unite("snp_for_join", c(position, alt_bp)) %>% .$snp_for_join %>% as.character()
      df_phenotype_binded_all_selected_annotated2 <- df_phenotype_binded_all_selected_annotated %>% select(-snp_for_join) %>% mutate(snp_for_join = snp_for_join)
      df2 <- df_phenotype_binded_all_selected_annotated2 %>% mutate(abs_dMOCCS2score = abs(dMOCCS2score)) %>% group_by(ID, snp_for_join, Antigen, Cell_type_class, Cell_type) %>% summarise(max_abs_dMOCCS2score = max(abs_dMOCCS2score))
      df3 <- df2 %>% left_join(GWAS_df2, by = "snp_for_join") 
      
      # ID -----
      x_lab_ID <- df3  %>% unite("ID_snp_for_join", c(ID, snp_for_join)) %>% filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% select(ID_snp_for_join)
      x_lab_ID2 <- x_lab_ID %>% separate(ID_snp_for_join, into = c("ID", "chr", "position", "alt_bp"), sep = "_") %>% .$ID %>% as.character() %>% rev()
      
      p1 <- df3  %>% unite("ID_snp_for_join", c(ID, snp_for_join)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
        ggplot(aes(x = reorder(ID_snp_for_join, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Antigen))+
        geom_col() +
        theme(plot.title = element_text(face="bold",hjust = 0.5), 
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=12,face="bold"),
              axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
              axis.text.y =element_text(size=10,face="bold"),
              axis.title=element_text(size=14,face="bold"),
              #legend.position = 'none',
              legend.title = element_blank()
        )+
        scale_x_discrete("ID", labels = x_lab_ID2) +
        ylab("max(|dMOCCS2score|)") +
        ggtitle(target_phenotype) +
        coord_flip()
      
      # SNP_rs label ----- 
      x_lab_ID <- readRDS("/Users/saeko/MOCCS_paper_public/data/Fig6/label/CD_label_TF_rs_added.rds")
      x_lab_ID2 <- x_lab_ID %>% .$SNP_rs %>% as.character() %>% rev()
      #x_lab_ID <- df3  %>% #unite("ID_snp_for_join", c(ID, snp_for_join)) %>% 
        #filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% select(SNP_rs)
      #x_lab_ID2 <- x_lab_ID %>% .$SNP_rs %>% as.character() %>% rev()
      #saveRDS(x_lab_ID, paste0("/Users/saeko/MOCCS_paper_public/data/Fig6/label/", target_phenotype, "_label_TF.rds"))
      p2 <- df3  %>% unite("ID_snp_for_join", c(ID, snp_for_join)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
        ggplot(aes(x = reorder(ID_snp_for_join, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Antigen))+
        geom_col() +
        theme(plot.title = element_text(face="bold",hjust = 0.5), 
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=12,face="bold"),
              axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
              axis.text.y =element_text(size=10,face="bold"),
              axis.title=element_text(size=14,face="bold"),
              #legend.position = 'none',
              legend.title = element_blank()
        )+
        scale_x_discrete("", labels = x_lab_ID2) +
        ylab("") +
        ggtitle(target_phenotype) +
        coord_flip()
      
      # ref k-mer label ----
      df4 <- df_phenotype_binded_all_selected_annotated2 %>% 
        mutate(abs_dMOCCS2score = abs(dMOCCS2score))%>% 
        select(ID, snp_for_join, Antigen, Cell_type_class, Cell_type, kmer_before_tou, kmer_after_tou, 
               dMOCCS2score, q_value, abs_dMOCCS2score) %>%
        group_by(ID, snp_for_join, Antigen, Cell_type_class, Cell_type) %>% 
        filter(abs_dMOCCS2score == max(abs_dMOCCS2score)) %>% distinct()
      colnames(df4) <-   c("ID", "snp_for_join", "Antigen", "Cell_type_class", "Cell_type", "ref_kmer", "alt_kmer", "dMOCCS2score",
                           "q_value", "max_abs_dMOCCS2score")
      
      x_lab_ID <- df4  %>% 
        filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% 
        select(ref_kmer) 
      x_lab_ID2 <- x_lab_ID %>% .$ref_kmer %>% as.character() %>% rev()
      
      p3 <- df4  %>% unite("ID_snp_for_join", c(ID, snp_for_join)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
        ggplot(aes(x = reorder(ID_snp_for_join, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Antigen))+
        geom_col() +
        theme(plot.title = element_text(face="bold",hjust = 0.5), 
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=12,face="bold"),
              axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
              axis.text.y =element_text(size=10,face="bold"),
              axis.title=element_text(size=14,face="bold"),
              #legend.position = 'none',
              legend.title = element_blank()
        )+
        scale_x_discrete("", labels = x_lab_ID2) +
        ylab("") +
        ggtitle(target_phenotype) +
        coord_flip()
      
      
      # 6C alt k-mer label ----
      x_lab_ID <- df4  %>% 
        filter(max_abs_dMOCCS2score > threshold) %>% arrange(desc(max_abs_dMOCCS2score)) %>% 
        select(alt_kmer) 
      x_lab_ID2 <- x_lab_ID %>% .$alt_kmer %>% as.character() %>% rev()
      p4 <- df4  %>% unite("ID_snp_for_join", c(ID, snp_for_join)) %>% filter(max_abs_dMOCCS2score > threshold) %>%
        ggplot(aes(x = reorder(ID_snp_for_join, max_abs_dMOCCS2score), y = max_abs_dMOCCS2score, fill = Antigen))+
        geom_col() +
        theme(plot.title = element_text(face="bold",hjust = 0.5), 
              panel.grid.major = element_line(colour = "gray"),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour="black"),
              axis.text=element_text(size=12,face="bold"),
              axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
              axis.text.y =element_text(size=10,face="bold"),
              axis.title=element_text(size=14,face="bold"),
              #legend.position = 'none',
              legend.title = element_blank()
        )+
        scale_x_discrete("", labels = x_lab_ID2) +
        ylab("") +
        ggtitle(target_phenotype) +
        coord_flip()
      
      
    }
    
  
  return(list(p1, p2, p3, p4))
}
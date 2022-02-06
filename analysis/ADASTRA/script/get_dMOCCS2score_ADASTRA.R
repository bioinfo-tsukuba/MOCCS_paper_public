get_dMOCCS2score <- function(target_ID = "SRX1041802",
                             target_tf = "GATA3",
                             input_bed_file_path = "/home/s-tahara/DROMPA/code_hg38_1/ftp.biosciencedbc.jp/archive/chip-atlas/data/hg38/eachData/bed05/",
                             pre_input_allele_binding_SNP_file_path = "/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/release_Susan/release_dump/TF_for_function/",
                             bed_output_dir = "/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/bed_output" ,
                             k_mer_num = 6,
                             ref_genome = "HG38",
                             input_MOCCS_file_dir = "/home/s-tahara/MOCCS-DB/WORK/MOCCS/",
                             fasta_file_name = "/home/s-tahara/MOCCS-DB/RESOURCE/UCSC/hg38/hg38.fa"
){
  
  library(tidyverse)
  
  # prepare target tf snp table
  input_allele_binding_SNP_file_path = paste0(pre_input_allele_binding_SNP_file_path, target_tf, "_for_function.tsv")
  df1 <- read_tsv(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/release_Susan/release_dump/TF/",target_tf,"_HUMAN.tsv"))
  df2 <- df1 %>% unite("snp", c(`#chr`, pos, ref, alt)) %>% mutate(tf = target_tf) %>% select(tf, snp, fdrp_bh_ref, fdrp_bh_alt)
  write_tsv(df2, paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/release_Susan/release_dump/TF_for_function/", target_tf, "_for_function.tsv"))
  
  kmer_all_bed_file <- paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/release_Susan/release_dump/TF_for_function/ADASTRA_SNP_all_for_dMOCCS2score_",target_tf, ".bed")
  
  if(file.exists(kmer_all_bed_file) == FALSE){
    source("/home/s-tahara/allele_binding_SNP/get_dMOCCS2score_ver5.R")
    SNP_table_5 <- get_dMOCCS2score(input_bed_file_path,input_allele_binding_SNP_file_path,bed_output_dir ,k_mer_num ,ref_genome,input_MOCCS_file_dir,fasta_file_name)
    saveRDS(SNP_table_5, paste0("/home/s-tahara/allele_binding_SNP/ADASTRA/release_Susan/release_dump/TF_for_function/ASTRADA_SNP_all_for_dMOCCS2score",target_tf,".rds"))
    write_tsv(SNP_table_5, paste0("/home/s-tahara/allele_binding_SNP/ADASTRA/release_Susan/release_dump/TF_for_function/ASTRADA_SNP_all_for_dMOCCS2score",target_tf,".tsv"), col_names = FALSE)
    system(paste0("cp /home/s-tahara/allele_binding_SNP/ADASTRA/release_Susan/release_dump/TF_for_function/ASTRADA_SNP_all_for_dMOCCS2score",target_tf,".tsv" , " ", kmer_all_bed_file))
  }
  
  # calculate dMOCCS2score
  W <- 350
  
  input_SRX_ID <- target_ID
  print(input_SRX_ID)
  if(file.exists(paste0(input_bed_file_path,input_SRX_ID, ".05.bed")) == FALSE){
    print(paste0("skip ", input_SRX_ID))
  }else{
    bed <- read_tsv(paste0(input_bed_file_path,input_SRX_ID, ".05.bed"), col_names = FALSE)
    
    if(nrow(bed)==0){
      print("no bed file")
    }else{
      ### 1. SRX ID
      bed_peak <- bed %>%  select(X1, X2, X3) # X1: chr1 etc., X2: 9900 etc., X3: 10436 etc.
      write_tsv(bed_peak, paste0(bed_output_dir, "/", input_SRX_ID, "_selected.bed"), col_names = FALSE) 
      
      selected_bed_file <- paste0(bed_output_dir, "/", input_SRX_ID, "_selected.bed")
      intersect_bed_file <- paste0(bed_output_dir, "/",input_SRX_ID, "_intersect.bed")
      output_fasta_file <- paste0(bed_output_dir, "/",input_SRX_ID, "_intersect.fasta")
      
      bedtools_exe_intersect_command <- paste0("~/miniconda3/bin/bedtools intersect -a ", kmer_all_bed_file, " -b ", selected_bed_file, " > ", intersect_bed_file)
      bedtools_exe_getfasta_command <- paste0("~/miniconda3/bin/bedtools getfasta -fi ", fasta_file_name, " -bed ", intersect_bed_file, " > ", output_fasta_file)
      
      system(bedtools_exe_intersect_command)
      system(bedtools_exe_getfasta_command)
      
      ### k-mer
      # snp table & k-mer
      intersect_bed <- read_tsv(paste0(bed_output_dir, "/",input_SRX_ID, "_intersect.bed"), col_names = FALSE)
      intersect_fasta <- read_tsv(paste0(bed_output_dir, "/",input_SRX_ID, "_intersect.fasta"), col_names = FALSE)
      if(nrow(intersect_bed) == 0){
        print("no common region")
      }else{
        
        # fasta
        intersect_fasta %>% filter(str_detect(X1, ">")) %>% .$X1 -> fasta_chr
        intersect_fasta %>% filter(!str_detect(X1, ">")) %>% .$X1 -> fasta_posi
        intersect_annotated <- intersect_bed %>% mutate(chr_fasta = fasta_chr, kmer_before = fasta_posi) %>%  filter(str_detect(kmer_before, "......"))
        write_tsv(intersect_annotated, paste0(bed_output_dir, "/", input_SRX_ID, "_intersect_annotated.tsv"))
        
        # MOCCS output (k-mer & MOCCS2score)
        genome <- ref_genome
        k <- k_mer_num
        SRX <- input_SRX_ID
        MOCCS_output <- read_tsv(paste0(input_MOCCS_file_dir, "MOCCS-OUTPUT-",genome,"-",k,"mer-ANALYSIS_v2/", SRX, "_",k ,"mer_v2.auc_count.txt"))
        
        # k-mer
        colnames(intersect_annotated) <- c("chr", "start","end", "ref_bp", "alt_bp", "fdrp_bh_ref", "fdrp_bh_alt", "posi", "chr_fasta", "kmer_before")
        kmer_before <- intersect_annotated$kmer_before %>% as.character()
        kmer_after <- c()
        for (x in 1:length(kmer_before)) {
          target_before_kmer <- kmer_before[x]
          snp_posi <- intersect_annotated[x,"posi"] %>% as.numeric()
          alt_bp <- intersect_annotated[x,"alt_bp"] %>% as.character()
          substring(target_before_kmer, snp_posi) <- alt_bp
          kmer_after <- c(kmer_after, target_before_kmer)
        }
        
        intersect_annotated3 <- intersect_annotated %>% mutate(kmer_after = kmer_after) %>% filter(!str_detect(kmer_after, ";"))
        
        # dMOCCS2score
        library(Biostrings)
        MOCCS2score_before <- c()
        MOCCS2score_after <- c()
        p_value_list <- c()
        
        for (y in 1:nrow(intersect_annotated3)) {
          target_kmer_before <- intersect_annotated3[y,"kmer_before"] %>% as.character() %>% toupper()
          target_kmer_after <- intersect_annotated3[y,"kmer_after"] %>% as.character() %>% toupper()
          target_MOCCS2score_before <- MOCCS_output %>% filter(kmer == target_kmer_before) %>% select(MOCCS2score) %>% as.numeric()
          target_MOCCS2score_after <- MOCCS_output %>% filter(kmer == target_kmer_after) %>% select(MOCCS2score) %>% as.numeric()
          target_auc_before <- MOCCS_output %>% filter(kmer == target_kmer_before) %>% select(auc) %>% as.numeric()
          target_auc_after <- MOCCS_output %>% filter(kmer == target_kmer_after) %>% select(auc) %>% as.numeric()
          target_count_before <- MOCCS_output %>% filter(kmer == target_kmer_before) %>% select(count) %>% as.numeric()
          target_count_after <- MOCCS_output %>% filter(kmer == target_kmer_after) %>% select(count) %>% as.numeric()
          
          if(is.na(target_MOCCS2score_before) == TRUE){
            target_kmer_before <- reverseComplement(DNAStringSet(target_kmer_before)) %>% as.character()
            target_MOCCS2score_before <- MOCCS_output %>% filter(kmer == target_kmer_before) %>% select(MOCCS2score) %>% as.numeric()
            target_auc_before <- MOCCS_output %>% filter(kmer == target_kmer_before) %>% select(auc) %>% as.numeric()
            target_count_before <- MOCCS_output %>% filter(kmer == target_kmer_before) %>% select(count) %>% as.numeric()
          }
          
          if(is.na(target_MOCCS2score_after) == TRUE){
            target_kmer_after <- reverseComplement(DNAStringSet(target_kmer_after)) %>% as.character()
            target_MOCCS2score_after <- MOCCS_output %>% filter(kmer == target_kmer_after) %>% select(MOCCS2score) %>% as.numeric()
            target_auc_after <- MOCCS_output %>% filter(kmer == target_kmer_after) %>% select(auc) %>% as.numeric()
            target_count_after <- MOCCS_output %>% filter(kmer == target_kmer_after) %>% select(count) %>% as.numeric()
          }
          
          var_i <- W^2/12/target_count_before
          var_j <- W^2/12/target_count_after
          
          # 2sample
          if(is.na(target_auc_before) == TRUE | is.na(target_auc_after) == TRUE){
            p_value_list <- c(p_value_list, -1)
            print("auc: na, skip!")
          }else if (target_auc_before >= target_auc_after){
            p <- 1 - pnorm((target_auc_before - target_auc_after)/sqrt(var_i + var_j) , 0, 1) 
            p_value_list <- c(p_value_list, p)
          }else{
            p <- 1 - pnorm((target_auc_after - target_auc_before)/sqrt(var_j + var_i) , 0, 1) 
            p_value_list <- c(p_value_list, p)
          }
          MOCCS2score_before <- c(MOCCS2score_before, target_MOCCS2score_before)
          MOCCS2score_after <- c(MOCCS2score_after, target_MOCCS2score_after)
        }
        q_value_list <- p.adjust(p_value_list)
        
        intersect_annotated4 <- intersect_annotated3 %>% mutate(kmer_before_tou = toupper(kmer_before),kmer_after_tou = toupper(kmer_after),MOCCS2score_before = MOCCS2score_before, MOCCS2score_after = MOCCS2score_after, p_value = p_value_list, q_value = q_value_list) %>% select(-kmer_before, -kmer_after) 
        intersect_annotated5 <- intersect_annotated4 %>% mutate(dMOCCS2score = MOCCS2score_before - MOCCS2score_after) %>% distinct()
        dMOCCS2score_df <- intersect_annotated5
        
      }# intersect_bed if
    } #
  } # bed,if
  
  return(dMOCCS2score_df)
  
}


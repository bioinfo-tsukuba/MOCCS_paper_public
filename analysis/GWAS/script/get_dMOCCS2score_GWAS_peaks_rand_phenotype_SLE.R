get_dMOCCS2score_GWAS_peaks_rand_phenotype <- function(target_phenotype, 
                                                       target_ID, 
                                                       input_bed_file_path, 
                                                       bed_output_dir,
                                                       fasta_file_name, 
                                                       input_MOCCS_file_dir, 
                                                       result_output_dir,
                                                       k_mer_num, 
                                                       simu_N, 
                                                       ref_genome){
  
  library(tidyverse)
  result <- list()
  
  # load data
  if(file.exists(paste0(input_bed_file_path,target_ID, ".05.bed")) == FALSE){
    print("no bed file")
  }else{
    bed <- suppressMessages(read_tsv(paste0(input_bed_file_path,target_ID, ".05.bed"), col_names = FALSE))
    kmer_all_bed_file_name <- paste0("/home/s-tahara/allele_binding_SNP_hg38/GWAS/data/GWAS_catalog_phenotype/SNP_all_for_function_", target_phenotype, ".rds")
    kmer_all_bed_file <- suppressMessages(readRDS(kmer_all_bed_file_name))
    kmer_all_bed_file_bed_name <- paste0("/home/s-tahara/allele_binding_SNP_hg38/GWAS/data/GWAS_catalog_phenotype/SNP_all_for_function_", target_phenotype, ".bed")
  
    # modulate chr name
    tmp <- kmer_all_bed_file 
    chr_new <- str_replace(tmp$chr, pattern="chr", replacement="")
    #chr_new <- paste0(rep("chr", nrow(tmp)), chr_new)
    tmp$chr <- suppressWarnings(as.numeric(chr_new))
    tmp2 <- tmp %>% drop_na(chr,start, end)
    #tmp2 <- tmp %>% filter(chr != "chrNA") %>% drop_na(chr)
    write_tsv(tmp2, paste0("/home/s-tahara/allele_binding_SNP_hg38/GWAS/data/GWAS_catalog_phenotype/SNP_all_for_function_", target_phenotype, ".bed"), col_names = FALSE)
  
    # get dMOCCS2score
    if(nrow(bed)==0){
      print("no bed file")
    }else{
      ### 1. get peak call bed file
      bed_peak <- bed %>% select(X1, X2, X3) # X1: chr1 etc., X2: 9900 etc., X3: 10436 etc.
      X1_new <- str_replace(bed_peak$X1, pattern="chr", replacement="") 
      X1_new <- suppressWarnings(as.numeric(X1_new))
      #X1_new <- paste0(rep("chr", length(X1_new)), X1_new)
      bed_peak$X1 <- X1_new
      bed_peak <- bed_peak %>% drop_na(X1)
      write_tsv(bed_peak, paste0(bed_output_dir, "/", target_ID, "_selected.bed"), col_names = FALSE) 
    
      ### 2. get fasta
      selected_bed_file <- paste0(bed_output_dir, "/", target_ID, "_selected.bed")
      intersect_bed_file <- paste0(bed_output_dir, "/",target_ID, "_intersect.bed")
      intersect_bed_file_name <- paste0(bed_output_dir, "/",target_ID, "_intersect.bed")
      output_fasta_file <- paste0(bed_output_dir, "/",target_ID, "_intersect.fasta")
    
      bedtools_exe_intersect_command <- paste0("/home/s-tahara/miniconda3/bin/bedtools intersect -a ",  selected_bed_file, " -b ",  kmer_all_bed_file_bed_name," > ", intersect_bed_file_name)
      bedtools_exe_getfasta_command <- paste0("/home/s-tahara/miniconda3/bin/bedtools getfasta -fi ", fasta_file_name, " -bed ", intersect_bed_file_name, " > ", output_fasta_file)
    
      # command preparation
      ## within peaks
      system(bedtools_exe_intersect_command)
      intersect_bed_file <- suppressMessages(read_tsv(paste0(bed_output_dir, "/",target_ID, "_intersect.bed"), col_names = FALSE))
      
      if(nrow(intersect_bed_file)==0){
        print("no common region in peaks")
      }else{
        intersect_bed_file$X1 <- paste0(rep("chr", nrow(intersect_bed_file)),intersect_bed_file$X1)
        write_tsv(intersect_bed_file, paste0(bed_output_dir, "/", target_ID, "_intersect.bed"), col_names = FALSE)
        intersect_bed <- suppressMessages(read_tsv(paste0(bed_output_dir, "/", target_ID, "_intersect.bed"), col_names = FALSE))
        system(bedtools_exe_getfasta_command)
        intersect_fasta <- suppressMessages(read_tsv(paste0(bed_output_dir, "/", target_ID, "_intersect.fasta"), col_names = FALSE))
      
        ## out of peaks
        # peak, random
        for (j in 1:simu_N) {
          non_intersect_bed_file <- paste0(bed_output_dir, "/", target_ID, "_non_intersect.bed")
          output_non_fasta_file <- paste0(bed_output_dir, "/", target_ID, "_" , j ,"_non_intersect.fasta")
          bedtools_exe_non_intersect_command <- paste0("/home/s-tahara/miniconda3/bin/bedtools intersect -a ",  kmer_all_bed_file_bed_name, " -b ", selected_bed_file, " -v > ", non_intersect_bed_file)
        
          system(bedtools_exe_non_intersect_command) 
          non_intersect_bed <- suppressMessages(read_tsv(paste0(bed_output_dir, "/", target_ID, "_non_intersect.bed"), col_names = FALSE))
        
          if (nrow(intersect_bed) < nrow(non_intersect_bed)){
            rand <- sample(nrow(intersect_bed))
            non_intersect_bed <- non_intersect_bed[rand,] 
          }
        
          non_intersect_bed$X1 <- paste0(rep("chr", nrow(non_intersect_bed)),non_intersect_bed$X1)
          non_intersect_bed <- non_intersect_bed %>% select(X1, X2, X3)  
        
          write_tsv(non_intersect_bed,paste0(bed_output_dir, "/", target_ID, "_" , j ,"_non_intersect2.bed"), col_names = FALSE)
          non_intersect_bed_file <- paste0(bed_output_dir, "/", target_ID, "_" , j ,"_non_intersect2.bed")
        
          bedtools_exe_non_getfasta_command <- paste0("/home/s-tahara/miniconda3/bin/bedtools getfasta -fi ", fasta_file_name, " -bed ", non_intersect_bed_file, " > ", output_non_fasta_file)
          system(bedtools_exe_non_getfasta_command)
          non_intersect_fasta <- suppressMessages(read_tsv(output_non_fasta_file, col_names = FALSE))
        
          ### 4. get dMOCCS2score 
          # 4-1 snp table & k-mer
          if(nrow(intersect_bed) == 0 | nrow(non_intersect_bed) == 0){
            print("no common region")
          }else{
          
            intersect_fasta %>% filter(str_detect(X1, ">")) %>% .$X1 -> fasta_chr
            intersect_fasta %>% filter(!str_detect(X1, ">")) %>% .$X1 -> fasta_posi
            non_intersect_fasta %>% filter(str_detect(X1, ">")) %>% .$X1 -> non_fasta_chr
            non_intersect_fasta %>% filter(!str_detect(X1, ">")) %>% .$X1 -> non_fasta_posi
          
            intersect_annotated <- intersect_bed %>% mutate(chr_fasta = fasta_chr, kmer_before = fasta_posi) %>%  
              #filter(str_detect(kmer_before, "......")) %>% 
              mutate(peak = "within peaks")
            non_intersect_annotated <- non_intersect_bed %>% mutate(chr_fasta = non_fasta_chr, kmer_before = non_fasta_posi) %>%  
              #filter(str_detect(kmer_before, "......")) %>% 
              mutate(peak = "without peaks")
            intersect_annotated <- rbind(intersect_annotated, non_intersect_annotated)
            snp_for_join <- intersect_annotated %>% unite("snp_for_join", c(X1, X2, X3), sep = "_") %>% .$snp_for_join %>% as.character()
            intersect_annotated2 <- intersect_annotated %>% mutate(snp_for_join = snp_for_join)
            GWAS_annotation <- suppressMessages(read_tsv(kmer_all_bed_file_bed_name, col_names = FALSE))
            GWAS_annotation$X1 <- paste0(rep("chr", nrow(GWAS_annotation)), GWAS_annotation$X1)
            snp_for_join <- GWAS_annotation %>% unite("snp_for_join", c(X1, X2, X3), sep = "_") %>% .$snp_for_join %>% as.character()
            GWAS_annotation2 <- GWAS_annotation %>%  mutate(snp_for_join = snp_for_join) %>% select(X4, X5, X6, X7, X8, X9, X10,snp_for_join)
            intersect_annotated3 <- intersect_annotated2 %>% left_join(GWAS_annotation2, by = "snp_for_join") %>% select(-chr_fasta, snp_for_join)
            write_tsv(intersect_annotated3, paste0(bed_output_dir, "/", target_ID, "_" , j , "_intersect_annotated.tsv"))
          
            # 4-2 MOCCS output (k-mer & MOCCS2score)
            genome <- ref_genome
            k <- k_mer_num
            SRX <- target_ID
            MOCCS_output <- suppressMessages(read_tsv(paste0(input_MOCCS_file_dir, "MOCCS-OUTPUT-",genome,"-",k,"mer-ANALYSIS_v2/", SRX, "_",k ,"mer_v2.auc_count.txt")))      #hg38, low, output
          
            # 4-3 k-mer MOCCS2score
            colnames(intersect_annotated3) <- c("chr", "start","end", "kmer_before", "peak", "snp_for_join","p_value", "OR" ,"Mapped_gene", "Trait", "snp", "position" ,"snp_posi")
          
            ## add alt-bp annotation 
            alt_bp_annotation <- suppressMessages(read_tsv("/home/s-tahara/allele_binding_SNP_hg38/GWAS/data/SNP_all_for_function.bed", col_names = FALSE))
            alt_bp_snp_for_join <- alt_bp_annotation %>% unite("snp_for_join", c(X1, X2, X3)) %>% .$snp_for_join %>% as.character()
            alt_bp_snp_for_join <- paste0(rep("chr", length(alt_bp_snp_for_join)), alt_bp_snp_for_join)
            alt_bp_annotation2 <- alt_bp_annotation %>% mutate(snp_for_join = alt_bp_snp_for_join) %>% select(snp_for_join, X4)
            colnames(alt_bp_annotation2) <- c("snp_for_join", "alt_bp")
            alt_bp_annotation3 <- alt_bp_annotation2 %>% filter(alt_bp == "?" | alt_bp == "A" | alt_bp == "C" | alt_bp == "T" | alt_bp == "G")
          
            intersect_annotated3_2 <- intersect_annotated3 %>% left_join(alt_bp_annotation3, by = "snp_for_join") %>% mutate(snp_posi_k = 6-(end-snp_posi)) %>% distinct()
          
            kmer_before <- intersect_annotated3_2$kmer_before %>% as.character()
            kmer_after <- c()
            for (z in 1:length(kmer_before)) {
              target_before_kmer <- kmer_before[z]
              snp_posi <- intersect_annotated3_2[z,"snp_posi_k"] %>% as.numeric()
              alt_bp <- intersect_annotated3_2[z,"alt_bp"] %>% as.character()
              substring(target_before_kmer, snp_posi) <- alt_bp
              kmer_after <- c(kmer_after, target_before_kmer)
            }
          
            intersect_annotated4 <- intersect_annotated3_2 %>% mutate(kmer_after = kmer_after) %>% filter(!str_detect(kmer_after, ";"))
            write_tsv(intersect_annotated4, paste0(bed_output_dir, "/", target_ID, "_" , j ,"_intersect_annotated3.tsv")) 
          
          
            # 4-4 MOCCS2score
            library(Biostrings)
            MOCCS2score_before <- c()
            MOCCS2score_after <- c()
            for (y in 1:nrow(intersect_annotated4)) {
              target_kmer_before <- intersect_annotated4[y,"kmer_before"] %>% as.character() %>% toupper()
              target_kmer_after <- intersect_annotated4[y,"kmer_after"] %>% as.character() %>% toupper()
              target_MOCCS2score_before <- MOCCS_output %>% filter(kmer == target_kmer_before) %>% select(MOCCS2score) %>% as.numeric()
              target_MOCCS2score_after <- MOCCS_output %>% filter(kmer == target_kmer_after) %>% select(MOCCS2score) %>% as.numeric()
            
              if(is.na(target_MOCCS2score_before) == TRUE){
                target_kmer_before <- reverseComplement(DNAStringSet(target_kmer_before)) %>% as.character()
                target_MOCCS2score_before <- MOCCS_output %>% filter(kmer == target_kmer_before) %>% select(MOCCS2score) %>% as.numeric()
              }
            
              if(is.na(target_MOCCS2score_after) == TRUE){
                if(str_detect(target_kmer_after, pattern = "\\?") == TRUE){
                  target_MOCCS2score_after <- NA
                }else if(str_detect(target_kmer_after, pattern = "I") == TRUE | str_detect(target_kmer_after, pattern = "<") == TRUE | str_detect(target_kmer_after, pattern = "_") == TRUE){ 
                  target_MOCCS2score_after <- NA
                }else{
                  target_kmer_after <- reverseComplement(DNAStringSet(target_kmer_after)) %>% as.character()
                  target_MOCCS2score_after <- MOCCS_output %>% filter(kmer == target_kmer_after) %>% select(MOCCS2score) %>% as.numeric()
                }
              }
              MOCCS2score_before <- c(MOCCS2score_before, target_MOCCS2score_before)
              MOCCS2score_after <- c(MOCCS2score_after, target_MOCCS2score_after)
            }
            intersect_annotated5 <- intersect_annotated4 %>% mutate(kmer_before_tou = toupper(kmer_before),kmer_after_tou = toupper(kmer_after),MOCCS2score_before = MOCCS2score_before, MOCCS2score_after = MOCCS2score_after) %>% select(-kmer_before, -kmer_after) 
            intersect_annotated6 <- intersect_annotated5 %>% mutate(dMOCCS2score = MOCCS2score_before - MOCCS2score_after) %>% distinct()
            result[[paste0(target_ID, "_" , j )]] <- intersect_annotated6
            #saveRDS(result, paste0(result_output_dir, target_ID , "_peak_rand.rds"))
            #dMOCCS2score_df[[paste0(target_ID, "_" , j )]] <- intersect_annotated6
          }# if, intersect_bed
        } #for rand (j)
      } #if
    } #if
  }
  return(result)
} #funcion
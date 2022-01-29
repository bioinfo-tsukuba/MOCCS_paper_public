# SNP-SELEX PBS (Fig.4B)
get_dMOCCS2score <- function(input_bed_file_path = "~/ftp.biosciencedbc.jp/archive/chip-atlas/data/hg19/eachData/bed05/",
                             input_allele_binding_SNP_file_path = "~/dMOCCS_files/file/GSE118725_pbs.novel_batch.tsv",
                             bed_output_dir = "~/dMOCCS_files",
                             k_mer_num = 6,
                             ref_genome = "HG19",
                             target_tf = "EBF",
                             input_SRX_ID = "SRX100432",
                             input_MOCCS_file_dir = "~/dMOCCS_files",
                             fasta_file_name = "~/allele_binding_SNP/hg19.fa"
){
  
  library(tidyverse)
  
  # 1. SRX ID
  print(input_SRX_ID)
  bed <- read_tsv(paste0(input_bed_file_path,input_SRX_ID, ".05.bed"), col_names = FALSE)
  if(nrow(bed)==0){
    print("no bed file")
  }else{
    
  bed_peak <- bed %>%  select(X1, X2, X3) # X1: chr1 etc., X2: 9900 etc., X3: 10436 etc.
  write_tsv(bed_peak, paste0(bed_output_dir, "/", input_SRX_ID, "_selected.bed"), col_names = FALSE)
  
  # 2. snp list
    
    SNP_table <- read_tsv(input_allele_binding_SNP_file_path)
    # annotation
    snp <- SNP_table$snp
    SNP_table2 <- SNP_table %>% separate(snp, into = c("chr", "number", "before", "after"), sep = "_") %>% mutate(snp = snp)
    SNP_table2$number <- SNP_table2$number %>% as.numeric() 
    SNP_table2 <- SNP_table2 %>% filter(after == "A" | after == "T" | after == "G" | after == "C")
    SNP_table3 <- SNP_table2 %>% select(chr, start, end, snp, number, tf, before, after, pbs, pval)
    SNP_table3 <- SNP_table3 %>% filter(tf == target_tf)
    
    if(nrow( SNP_table3) == 0){
      print("no target TF")
    }else{
     
    
    k <- k_mer_num
    k_1 <- k - 1
    SNP_table3 <- SNP_table3 %>% filter(pbs != "-" & pbs != "+")
    N <- nrow(SNP_table3)
    
    # snp
    library(pbapply)
    SNP_table4 <- bind_rows(
      pblapply(1:N, function(x){
        target_row <- SNP_table3[x,] %>% as_tibble()
        table <- as_tibble(target_row)
        if(nrow(table) == 1){
          table$pbs <- as.numeric(table$pbs)
          table$pval <- as.numeric(table$pval)
        }
        
        for (y in 1:k_1) {
          chr <- target_row$chr %>% as.character()
          start <- target_row$start %>% as.numeric() + y
          end <- target_row$end %>% as.numeric() + y
          target_snp <- target_row$snp %>% as.character()
          number <- target_row$number %>% as.numeric()
          tf <- target_row$tf %>% as.character()
          before <- target_row$before %>% as.character()
          after <- target_row$after %>% as.character()
          pbs <- target_row$pbs  %>% as.numeric()
          pval <- target_row$pval %>% as.numeric()
          add_tib <- tibble(chr = chr, start = start, end = end, snp = target_snp, number = number, tf = tf, before = before, after = after, pbs = pbs, pval = pval)
          table <- table %>% add_row(add_tib)
        }#y
        return(table)
      })#lapply
    )#bind_rows

    posi <- rep(k:1, N)
    SNP_table5 <- SNP_table4 %>% mutate(snp_posi = posi) %>% select(-snp, -number) %>% mutate(number = end + snp_posi - 6)
    write_tsv(SNP_table5, paste0(bed_output_dir, "/",input_SRX_ID,"_kmer_all.bed"), col_names = FALSE) # TO BE MODIFIED about file preset
  #}
  
    
    # 3.fasta (shell)
    kmer_all_bed_file <- paste0(bed_output_dir, "/",input_SRX_ID,"_kmer_all.bed")
    selected_bed_file <- paste0(bed_output_dir, "/", input_SRX_ID, "_selected.bed")
    intersect_bed_file <- paste0(bed_output_dir, "/",input_SRX_ID, "_intersect.bed")
    output_fasta_file <- paste0(bed_output_dir, "/",input_SRX_ID, "_intersect.fasta")
    
    bedtools_exe_intersect_command <- paste0("~/miniconda3/bin/bedtools intersect -a ", kmer_all_bed_file, " -b ", selected_bed_file, " > ", intersect_bed_file)
    bedtools_exe_getfasta_command <- paste0("~/miniconda3/bin/bedtools getfasta -fi ", fasta_file_name, " -bed ", intersect_bed_file, " > ", output_fasta_file)
    
    #system("~/miniconda3/bin/bedtools")
    system(bedtools_exe_intersect_command)
    system(bedtools_exe_getfasta_command)
    
    # 4. k-mer, MOCCS2score
    
    # 4-1 snp table & k-mer
    intersect_bed <- read_tsv(paste0(bed_output_dir, "/",input_SRX_ID, "_intersect.bed"), col_names = FALSE)
    intersect_fasta <- read_tsv(paste0(bed_output_dir, "/",input_SRX_ID, "_intersect.fasta"), col_names = FALSE)
    if(nrow(intersect_bed) == 0){
      print("no common region")
    }else{
    intersect_bed <- intersect_bed %>% filter(X6 != "TRUE")
    if(nrow(intersect_bed) != 0){
    
    intersect_fasta %>% filter(str_detect(X1, ">")) %>% .$X1 -> fasta_chr
    intersect_fasta %>% filter(!str_detect(X1, ">")) %>% .$X1 -> fasta_posi
    intersect_annotated <- intersect_bed %>% mutate(chr_fasta = fasta_chr, kmer_before = fasta_posi) %>%  filter(str_detect(kmer_before, "......"))
    write_tsv(intersect_annotated, paste0(bed_output_dir, "/", input_SRX_ID, "_intersect_annotated.tsv"))
    
    # 4-2 MOCCS output (k-mer & MOCCS2score)
    genome <- ref_genome
    k <- k_mer_num
    SRX <- input_SRX_ID
    MOCCS_output <- read_tsv(paste0("/home/s-tahara/MOCCS-DB/WORK/MOCCS/MOCCS-OUTPUT-", genome, "-6mer-ANALYSIS_v2/",SRX, "_6mer_v2.auc_count.txt")) # --lowthreshold option
        
    # 4-3 k-mer, MOCCS2score
    # kmer_after
    intersect_annotated2 <- intersect_annotated
    colnames(intersect_annotated2) <- c("chr", "start", "end", "tf", "before_bp", "after_bp","pbs","pval","posi","snp_position" ,"chr_fasta", "kmer_before")
    kmer_before <- intersect_annotated2$kmer_before %>% as.character()
    kmer_after <- c()
    for (x in 1:length(kmer_before)) {
      target_before_kmer <- kmer_before[x]
      snp_posi <- intersect_annotated2[x,"posi"] %>% as.numeric()
      after_bp <- intersect_annotated2[x,"after_bp"] %>% as.character()
      if(after_bp == "A" | after_bp == "T" | after_bp == "G" | after_bp == "C"){
        substring(target_before_kmer, snp_posi) <- after_bp
        kmer_after <- c(kmer_after, target_before_kmer)
      }
    }
    
    intersect_annotated3 <- intersect_annotated2 %>% mutate(kmer_after = kmer_after)
    write_tsv(intersect_annotated3, paste0(bed_output_dir, "/", input_SRX_ID, "_intersect_annotated3.tsv"))
    
    # 4-4 MOCCS2score
    library(Biostrings)
    MOCCS2score_before <- c()
    MOCCS2score_after <- c()
    for (y in 1:nrow(intersect_annotated3)) {
      target_kmer_before <- intersect_annotated3[y,"kmer_before"] %>% as.character() %>% toupper()
      target_kmer_after <- intersect_annotated3[y,"kmer_after"] %>% as.character() %>% toupper()
      target_MOCCS2score_before <- MOCCS_output %>% filter(kmer == target_kmer_before) %>% select(MOCCS2score) %>% as.numeric()
      target_MOCCS2score_after <- MOCCS_output %>% filter(kmer == target_kmer_after) %>% select(MOCCS2score) %>% as.numeric()
      
      if(is.na(target_MOCCS2score_before) == TRUE){
        target_kmer_before <- reverseComplement(DNAStringSet(target_kmer_before)) %>% as.character()
        target_MOCCS2score_before <- MOCCS_output %>% filter(kmer == target_kmer_before) %>% select(MOCCS2score) %>% as.numeric()
      }
      if(is.na(target_MOCCS2score_after) == TRUE){
        target_kmer_after <- reverseComplement(DNAStringSet(target_kmer_after)) %>% as.character()
        target_MOCCS2score_after <- MOCCS_output %>% filter(kmer == target_kmer_after) %>% select(MOCCS2score) %>% as.numeric()
      }
      MOCCS2score_before <- c(MOCCS2score_before, target_MOCCS2score_before)
      MOCCS2score_after <- c(MOCCS2score_after, target_MOCCS2score_after)
    }
    intersect_annotated4 <- intersect_annotated3 %>% mutate(kmer_before_tou = toupper(kmer_before),kmer_after_tou = toupper(kmer_after),MOCCS2score_before = MOCCS2score_before, MOCCS2score_after = MOCCS2score_after) %>% select(-kmer_before, -kmer_after) 
    intersect_annotated5 <- intersect_annotated4 %>% mutate(dMOCCS2score = MOCCS2score_before - MOCCS2score_after) %>% distinct()
    write_tsv(intersect_annotated5, paste0(bed_output_dir, "/", input_SRX_ID, "_intersect_annotated5.tsv"))
    
    return(intersect_annotated5)
    }#内側のif
    }# intersect_bed, if
    
    }# target tf , if
  }#bed , if
}

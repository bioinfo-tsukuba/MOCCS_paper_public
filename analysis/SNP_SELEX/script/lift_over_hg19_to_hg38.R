# lift over of SNP_SELEX data
input_allele_binding_SNP_file_path = "/home/s-tahara/allele_binding_SNP/file/GSE118725_pbs.novel_batch.tsv"
SNP_table <- read_tsv(input_allele_binding_SNP_file_path)
# annotation
snp <- SNP_table$snp
SNP_table2 <- SNP_table %>% separate(snp, into = c("chr", "number", "before", "after"), sep = "_") %>% mutate(snp = snp)
SNP_table2$number <- SNP_table2$number %>% as.numeric() 
SNP_table_3 <- SNP_table2 %>% mutate(number_start = number-6) %>% select(chr, number_start, number, snp, tf, pbs, pval) 
colnames(SNP_table_3) <- c("chr", "start", "end", "snp", "tf", "pbs", "pval")
write_tsv(SNP_table_3, "/home/s-tahara/allele_binding_SNP_hg38/SNP_SELEX/data/GSE118725_pbs_novel_batch_hg19.bed" , col_names = FALSE)

#tmp_bed <- read_tsv("/home/s-tahara/allele_binding_SNP/file/GSE118725_pbs.novel_batch.bed", col_names = FALSE)
#tmp_bed2 <- tmp_bed %>% separate(X2, c("chr", "start", ""))
#k_mer_all_bed_file_path = "/home/s-tahara/allele_binding_SNP/file/GSE118725_kmer_all.bed"

# after lift over
input_allele_binding_SNP_file_path2 = "/home/s-tahara/allele_binding_SNP_hg38/SNP_SELEX/data/GSE118725_pbs_novel_batch_hg38.bed"
SNP_table <- read_tsv(input_allele_binding_SNP_file_path2, col_names = FALSE)
colnames(SNP_table) <- c("chr", "start", "end", "snp", "tf", "pbs", "pval")
write_tsv(SNP_table, "/home/s-tahara/allele_binding_SNP_hg38/SNP_SELEX/data/GSE118725_pbs_novel_batch_hg38.tsv")

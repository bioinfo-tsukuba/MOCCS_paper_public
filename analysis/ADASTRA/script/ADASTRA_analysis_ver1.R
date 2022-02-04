target_tf <- "GATA3"
ID_list <- read_tsv(paste0("/home/s-tahara/singularity_ADASTRA/ID_list/", target_tf, ".txt"))
df_all <- readRDS(paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/dMOCCS2score_output_all/ADASTRA_dMOCCS2score_", target_tf, "_all.rds"))

# 1. number of IDs 
length(ID_list)
length(unique(df_all$ID))

# 2. number of snps
df_all2 <- df_all %>% mutate(position = end - 6 + posi) 
snp <- df_all2 %>% unite("snp", c(chr, position)) %>% .$snp %>% unique() %>% as.character()
length(unique(snp))

# 3. ASB significance summary
df3 <- df_all2 %>% mutate(fdrp_bh_ref_log = -log10(fdrp_bh_ref),
                          fdrp_bh_alt_log = -log10(fdrp_bh_alt), 
                          ASB = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, fdrp_bh_alt_log, fdrp_bh_ref_log), 
                          x_axis = ifelse(fdrp_bh_ref_log <= fdrp_bh_alt_log, "positive", "negative"))
df4 <- df3 %>% filter(x_axis == "positive") %>% mutate(ASB_plot = ASB)
df5 <- df3 %>% filter(x_axis == "negative") %>% mutate(ASB_plot = -ASB)
df6 <- rbind(df4, df5)
summary(df6$ASB_plot)

# 4. dMOCCS2score summary
summary(df$dMOCCS2score)

# 5. number of differential k-mer

# 6. number of concordant and disconcordant snps

# 7. negative control summary

# 8. p value from negative control's experimental distribution
# define soft filter  IDs
annotation <- readRDS("/home/s-tahara/DROMPA/code_hg38_2/experimentList_tab4.rds")
hg38_annotation <- annotation %>% filter(Genome == "hg38" & Antigen_class == "TFs and others")
Antigen_list <- read_tsv("/home/s-tahara/DROMPA/code_hg38_2/Antigen_list_hg19.txt", col_names = FALSE)
Antigen_list <- Antigen_list$X1 %>% as.character()
hg38_annotation <- hg38_annotation %>% filter(ID %in% Antigen_list)
length(unique(hg38_annotation$ID))

# 1.read
summary(hg38_annotation$reads)
hg38_annotation %>%
  ggplot() +
  geom_histogram(aes(x = reads), bins = 100) + 
  #xlim(0, 100000)
  labs(title = "hg38_#_of_reads_1E-05") 

hg38_annotation %>%
  ggplot() +
  geom_histogram(aes(x = log(reads)), bins = 100) + 
  #xlim(0, 100000)
  labs(title = "hg38_#_of_reads_1E-05") 

# 2.peak
summary(hg38_annotation$peaks)
hg38_annotation %>%
  ggplot() +
  geom_histogram(aes(x = peaks), bins = 100) + 
  labs(title = "hg38_#_of_peaks_1E-05")

hg38_annotation %>%
  ggplot() +
  geom_histogram(aes(x = log(peaks)), bins = 100) + 
  labs(title = "hg38_#_of_peaks_1E-05")

# 3. add %mapped from bowtie2 output
mapped_df1 <- read_tsv("/home/s-tahara/DROMPA/code_hg38_2/mapped_hg38.txt", col_names = FALSE)
colnames(mapped_df1) <- c("log", "mapped_bowtie2")
mapped_df1$mapped_bowtie2 <- gsub('%', '', mapped_df1$mapped_bowtie2) %>% as.numeric()

array_SRX_df <- read_tsv("/home/s-tahara/DROMPA/code_hg38_2/array_log_SRX.txt", col_names = FALSE)
colnames(array_SRX_df) <- c("log", "ID")
array_SRX_df$log <- gsub('drompa_ssp_hg38_ddbj_v2.sh.o', 'drompa_ssp_hg38_ddbj_v2.sh.e', array_SRX_df$log)
mapped_df2 <- mapped_df1 %>% left_join(array_SRX_df, by = "log")
mapped_df2$ID <- gsub("/home/s-tahara/DROMPA/FASTQ/", "", mapped_df2$ID)
mapped_df3 <- mapped_df2 %>% select(ID, mapped_bowtie2)
hg38_annotation2 <- hg38_annotation %>% left_join(mapped_df3, by = "ID")
hg38_annotation2$mapped_bowtie2[is.na(hg38_annotation2$mapped_bowtie2)] <- 0
hg38_annotation2 <- hg38_annotation2 %>% filter(mapped_bowtie2 != 0)


# 4.define mapping ratio threshold
threshold1 <- mean(hg38_annotation2$mapped_bowtie2) -2*sd(hg38_annotation2$mapped_bowtie2)

summary(hg38_annotation2$mapped_bowtie2)
hg38_annotation2 %>%
  ggplot() +
  geom_histogram(aes(x = mapped_bowtie2), bins = 100) + 
  labs(title = "hg38_%_of_mapped_1E-05_threshold1") +
  geom_vline(xintercept = threshold1)+
  theme(plot.title = element_text(face="bold",hjust = 0.5), 
        panel.grid.major = element_line(colour = "gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(), 
        axis.line = element_line(colour="black"),
        axis.text=element_text(size=12,face="bold"),
        axis.text.x =element_text(size=6,face="bold", angle = 45, hjust = 1),
        axis.text.y =element_text(size=10,face="bold"),
        axis.title=element_text(size=14,face="bold"),
        #legend.position = 'none',
        legend.title = element_blank()
  )



# 5. get soft filter ID
nrow(hg38_annotation2)
hg38_annotation_soft <- hg38_annotation2 %>% filter(mapped_bowtie2 > threshold1 & reads > 10000000 & peaks > 100) %>% distinct()
nrow(hg38_annotation_soft)

summary(hg38_annotation_soft$reads)
hg38_annotation_soft %>% ggplot() + geom_histogram(aes(x = reads))
summary(hg38_annotation_soft$peaks)
hg38_annotation_soft %>% ggplot() + geom_histogram(aes(x = peaks))
summary(hg38_annotation_soft$mapped_bowtie2)
hg38_annotation_soft %>% ggplot() + geom_histogram(aes(x = mapped_bowtie2))

ID_soft_filter_hg38 <- hg38_annotation_soft$ID %>% as.character()
length(ID_soft_filter_hg38)
write(ID_soft_filter_hg38, "/home/s-tahara/soft_filter/ID_soft_filter_hg38.txt")
saveRDS(ID_soft_filter_hg38, "/home/s-tahara/soft_filter/ID_soft_filter_hg38.rds")


# see relation of variable
g_pair <- hg38_annotation2 %>% select(reads, peaks, mapped_bowtie2) %>% ggpairs(lower=list(continuous=wrap("points",size=0.05)))+
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
g_pair


g_pair_soft <- hg38_annotation_soft %>% select(reads, peaks, mapped_bowtie2) %>% ggpairs(lower=list(continuous=wrap("points",size=0.05)))+
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
g_pair_soft
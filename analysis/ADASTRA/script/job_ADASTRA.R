# シェルスクリプトからSGE_TASK_IDを引き継ぐ
args <- commandArgs(TRUE)
target_ID <- args[1] # 引き継いだSGE_TASK_ID
target_tf <- args[2] # 引き継いだAntigen

# 関数を実行する
source("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/script/get_dMOCCS2score_ADASTRA.R")

input_bed_file_path = "/home/s-tahara/DROMPA/code_hg38_1/ftp.biosciencedbc.jp/archive/chip-atlas/data/hg38/eachData/bed05/"
pre_input_allele_binding_SNP_file_path = "/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/release_Susan/release_dump/TF_for_function/"
bed_output_dir = "/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/bed_output" 
k_mer_num = 6
ref_genome = "HG38"
input_MOCCS_file_dir = "/home/s-tahara/MOCCS-DB/WORK/MOCCS/"
fasta_file_name = "/home/s-tahara/MOCCS-DB/RESOURCE/UCSC/hg38/hg38.fa"

df <- get_dMOCCS2score(target_ID, 
                       target_tf, 
                       input_bed_file_path, 
                       pre_input_allele_binding_SNP_file_path, 
                       bed_output_dir, 
                       k_mer_num, 
                       ref_genome,
                       input_MOCCS_file_dir,
                       fasta_file_name)

saveRDS(df, paste0("/home/s-tahara/allele_binding_SNP_hg38/ADASTRA/data/dMOCCS2score_output/ADASTRA_dMOCCS2score_",target_tf, "_", target_ID,"_qval.rds"))






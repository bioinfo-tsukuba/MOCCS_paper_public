## FOR NIG
args <- commandArgs(TRUE)

if (length(args) == 1){
  eval(parse(text = args)) # $SGE_TASK_ID is substituted to x
} else {
  stop ()
}

library(dplyr)

# Read data
df_p_1 <- readRDS("/home/takaho-tsuchiya/MOCCS-DB/QSUB/2022_02_05/INPUT/RDS/df_p_1.rds")
# df_p_1 <- readRDS("~/MOCCS-DB_paper/results/df_p_1.rds")

bed_dir <- "/home/takaho-tsuchiya/MOCCS-DB/QSUB/2022_02_05/INPUT/BED"
# bed_dir <- "~/MOCCS-DB_paper/function/tmp_input/BED/hg19/eachData/bed05"

# Function for calculating POi between given pairs
Get_POi_each <- function(ID1, ID2, bed_dir){
  
  # Prepare get n1 / n2 / n1all / n2all commands as character variables
  bed_command_1 <- "bedtools intersect -u -a "
  bed_command_2 <- " -b "
  bed_command_3 <- " | wc -l"
  
  cat_command_1 <- "cat "
  cat_command_2 <- " | wc -l"
  
  ID1_full_path <- paste0(bed_dir, "/", ID1, ".05.bed")
  ID2_full_path <- paste0(bed_dir, "/", ID2, ".05.bed")
  
  n1_command <- paste0(bed_command_1, ID1_full_path, bed_command_2, ID2_full_path, bed_command_3)
  n2_command <- paste0(bed_command_1, ID2_full_path, bed_command_2, ID1_full_path, bed_command_3)
  
  n1all_command <- paste0(cat_command_1, ID1_full_path, cat_command_2)
  n2all_command <- paste0(cat_command_1, ID2_full_path, cat_command_2)
  
  # Do get commands
  n1 <- as.double(system(n1_command, intern = TRUE))
  n2 <- as.double(system(n2_command, intern = TRUE))
  
  n1all <- as.double(system(n1all_command, intern = TRUE))
  n2all <- as.double(system(n2all_command, intern = TRUE))
  
  n1_by_n1all <- n1 / n1all
  n2_by_n2all <- n2 / n2all
  n1_by_n1all_plus_n2_by_n2all <- n1_by_n1all + n2_by_n2all
  
  return(list(n1 = n1,
              n2 = n2,
              n1all = n1all,
              n2all = n2all,
              n1_by_n1all = n1_by_n1all,
              n2_by_n2all = n2_by_n2all,
              n1_by_n1all_plus_n2_by_n2all = n1_by_n1all_plus_n2_by_n2all))
  
}

# Function for calculating POi of specified pairs
Get_POi_pairs <- function(df_p_1, bed_dir, x){
  
  library(dplyr)
  
  pair_num <- nrow(df_p_1)
  
  job_num <- 200
  range_num <-  pair_num / job_num
  
  start_ind_vec <- range_num * (seq(1:job_num) - 1) + 1
  end_ind_vec <- range_num * seq(1:job_num)
  
  colnames_df_poi <- c("ID1", "ID2",
                       "n1", "n2", "n1all", "n2all",
                       "n1_by_n1all", "n2_by_n2all",
                       "n1_by_n1all_plus_n2_by_n2all")
  df_poi <- matrix(0, nrow = pair_num, ncol = length(colnames_df_poi))
  colnames(df_poi) <- colnames_df_poi
  
  for (target_pair in start_ind_vec[x]:end_ind_vec[x]){
    # for (target_pair in 1:1000){ # for debug
    
    id1 <- df_p_1[target_pair, "ID1"]
    id2 <- df_p_1[target_pair, "ID2"]
    
    # Get POi
    poi_each <- Get_POi_each(id1, id2, bed_dir)
    
    # Collect
    df_poi[target_pair, "ID1"] <- as.character(id1)
    df_poi[target_pair, "ID2"] <- as.character(id2)
    
    df_poi[target_pair, "n1"] <- poi_each$n1 %>% as.double()
    df_poi[target_pair, "n2"] <- poi_each$n2
    df_poi[target_pair, "n1all"] <- poi_each$n1all
    df_poi[target_pair, "n2all"] <- poi_each$n2all
    df_poi[target_pair, "n1_by_n1all"] <- poi_each$n1_by_n1all
    df_poi[target_pair, "n2_by_n2all"] <- poi_each$n2_by_n2all
    df_poi[target_pair, "n1_by_n1all_plus_n2_by_n2all"] <- poi_each$n1_by_n1all_plus_n2_by_n2all
  
    # if (target_pair %% 100 == 0){
    #  print(target_pair)
    # }

  }
  
  df_poi <- as_tibble(df_poi)
  df_poi[,"n1"] <- as.double(df_poi$n1)
  df_poi[,"n2"] <- as.double(df_poi$n2)
  df_poi[,"n1all"] <- as.double(df_poi$n1all)
  df_poi[,"n2all"] <- as.double(df_poi$n2all)
  df_poi[,"n1_by_n1all"] <- as.double(df_poi$n1_by_n1all)
  df_poi[,"n2_by_n2all"] <- as.double(df_poi$n2_by_n2all)
  df_poi[,"n1_by_n1all_plus_n2_by_n2all"] <- as.double(df_poi$n1_by_n1all_plus_n2_by_n2all)
  
  return(df_poi)
  
}

# x <- 1 # tmp
df_poi <- Get_POi_pairs(df_p_1, bed_dir, x)

saveRDS(df_poi, paste0("/home/takaho-tsuchiya/MOCCS-DB/QSUB/2022_02_05/OUTPUT/2022_02_05_df_poi_", x, ".rds"))
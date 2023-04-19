Collect_poi <- function(){
  
  pair_num <- 4426800
  job_num <- 200
  range_num <-  pair_num / job_num
  
  start_ind_vec <- range_num * (seq(1:job_num) - 1) + 1
  end_ind_vec <- range_num * seq(1:job_num)
  
  data_pre <- "~/MOCCS-DB_paper/data/Fig2/poi/df_poi_"
  df_poi <- c()

  for (i in 1:job_num){
  # for (i in 1:20){
    print(i)
    data_path <- paste0(paste0(data_pre, i), ".rds")
    data_each <- readRDS(data_path)
    df_poi <- rbind(df_poi, data_each[start_ind_vec[i]:end_ind_vec[i], ])
  }
  
  return(df_poi)
  
}
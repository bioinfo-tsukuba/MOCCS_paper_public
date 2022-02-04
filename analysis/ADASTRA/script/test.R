# 実行するR script側(仮にここではexe.Rとする)の1行目に下記のコードを挿入する
# xに$SGE_TASK_IDが格納される
library(tidyverse)
print("start now!!")
args <- commandArgs(TRUE)

#if (length(args) == 1){
  #eval(parse(text = args))
#} else {
  #stop ()
#}
print(paste0("args : ", args))
#print(paste0("x : ", x))

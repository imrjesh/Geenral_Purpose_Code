library(readxl)
library(tidyr)
setwd("/Users/kumarr9/Desktop/rajesh_projects/fusion_gene/arriba_detected/arriba_high_confidence/arriba_sorted_discordant")
getwd()
my_files <- list.files(pattern = ".xls")
my_files
rajesh = lapply(my_files, function(i){
  x = read_excel(i, sheet = 1)
  x$my_file = i
  x
})
rajesh[[1]]
rajesh = do.call("rbind.data.frame", rajesh)
write.table(rajesh, file ="/Users/kumarr9/Desktop/all_fusion_rajesh.csv", row.names = FALSE, col.names = TRUE, sep = ',', append = FALSE);

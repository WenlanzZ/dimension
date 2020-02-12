## code to prepare `lung` dataset goes here
lung <- read.table("/Users/wz262/Projects/dimension/data-raw/lung.txt", sep="\t", header=TRUE)
lung <- as.matrix(lung)
usethis::use_data(lung, overwrite = TRUE)

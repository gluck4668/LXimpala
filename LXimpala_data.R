
library(openxlsx)
meta_data_example <- read.csv("meta_ORA_results.csv")
gene_meta_example <- read.csv("gene_meta_ORA_results.csv")
protein_meta_example <- read.csv("pro_meta_ORA_results.csv")

usethis::use_data(meta_data_example,overwrite = T)
usethis::use_data(gene_meta_example,overwrite = T)
usethis::use_data(protein_meta_example,overwrite = T)

rm(list=ls())

data(meta_data_example)
data(gene_meta_example)
data(protein_meta_example)

require(tidyverse)
require(DESeq2)

metadata <- read.csv("../code/sample-metadata2.tsv",sep = "\t")

Function_table <- read.csv("../results/8.DEG/q2-picrust2_output/Functional_profiling_ko/pred_metagenome_unstrat_descrip.tsv",sep = "\t",check.names = F) %>% 
  dplyr::select(-`function`) %>% mutate(Function = str_extract(description,pattern = "(?<=;)(.*)(?=)"))


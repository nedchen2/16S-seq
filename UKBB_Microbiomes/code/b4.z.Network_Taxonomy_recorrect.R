require(tidyverse)

taxonomy_original <- read.csv(file ="../results/4.Diversity_ana/taxonomy-corrected.tsv",sep = "\t")

taxonomy_blast <- read.csv(file = "../results/blast/0.97-5hit/taxonomy.tsv",sep = "\t")

taxonomy_blast_2 <- read.csv(file = "../results/blast/0.99-1hit/taxonomy.tsv",sep = "\t")

combine_taxonomy <- left_join(taxonomy_original,taxonomy_blast,by="Feature.ID")

write_tsv(combine_taxonomy,"../results/6.Network_analysis/test.tsv")

uninomial_to_bionomial <- function(test_taxa){
  species <- str_extract(test_taxa,"(?<=s__)(.*)")
  
  print("Species assigned")
  
  genus <- str_extract(test_taxa,"(?<=g__)(.*)(?=; s__?)")
  
  bionomial_nomen <- paste0(genus,"_",species)
  
  taxonomy_before <- str_extract(test_taxa,"(.*)(?=; s__?)")
  
  new_taxonomy <- paste0(taxonomy_before,"; s__",bionomial_nomen)
  
  return(new_taxonomy)
}


ICare <- c("Snodgrassella","Gilliamella","Arsenophonus")

Combine_unassigned_0.97 <- function(){
  taxonomy_blast <- read.delim ("../results/blast/0.97-5hit/taxonomy.tsv")
  colnames(taxonomy_blast) <- c("Feature.ID","Taxon_blast","Consensus")
  idx_need_binomial <- str_detect(taxonomy_blast$Taxon_blast,"s__")
  taxonomy_blast[idx_need_binomial,"Taxon_blast"] = uninomial_to_bionomial(taxonomy_blast$Taxon_blast[idx_need_binomial])
  
  df_test2 <- left_join(taxonomy_original,taxonomy_blast,by="Feature.ID") 
  
  df_final2 <- df_test2 %>% mutate(Taxon_blast_revised =str_replace(str_replace(str_replace(str_replace(Taxon_blast,"k__","d__"),"c__Betaproteobacteria; o__Burkholderiales","c__Gammaproteobacteria; o__Burkholderiales"),"g__Escherichia","g__Escherichia-Shigella"),"g__Shigella","g__Escherichia-Shigella")) %>% 
    mutate(Taxon_blast_revised=str_replace(Taxon_blast_revised,"s__Snodgrassella_alvi wkB2","s__Snodgrassella_alvi"))%>%
    mutate(Taxon_blast_revised=str_replace(Taxon_blast_revised,"s__Phytobacter_massiliensis JC163","s__Phytobacter_massiliensis")) %>%
    mutate(Taxon_blast_revised=str_replace(Taxon_blast_revised,"s__Rhizorhapis_suberifaciens NBRC 15211","s__Rhizorhapis_suberifaciens"))%>%
    mutate(Taxon_blast_revised=str_replace(Taxon_blast_revised,"d__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Neisseriales; f__Neisseriaceae; g__Snodgrassella","d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Neisseriaceae; g__Snodgrassella")) %>%
    dplyr::select("Feature.ID","Taxon","Confidence","Taxon_blast_revised","Consensus")
  
  
  
  count = 0 
  for (i in seq(nrow(df_final2))){
    if (grepl("g__", df_final2[i,"Taxon"] ,ignore.case = T)){
      if (grepl("s__", df_final2[i,"Taxon"] ,ignore.case = T )){
        genus_original = str_extract(df_final2[i,"Taxon"],pattern = "(.*)(?=; s__?)")
        genus_blast = str_extract(df_final2[i,"Taxon_blast_revised"],pattern = "(.*)(?=; s__?)")
      }else{
        genus_original = df_final2[i,"Taxon"]
        genus_blast = str_extract(df_final2[i,"Taxon_blast_revised"],pattern = "(.*)(?=; s__?)")
      }
      
      idx = !is.na(genus_blast == genus_original)
      idx0 =  genus_blast== genus_original
      
      idx1= idx0 && idx
      idx2 = grepl("s__", df_final2[i,"Taxon_blast_revised"] ,ignore.case = T )
      if (idx1 && idx2){
        print("we could remove this")
        count = count + 1
        print (df_final2[i,"Taxon_blast_revised"])
        df_final2[i,"Taxon"] = df_final2[i,"Taxon_blast_revised"]
      }
    }
  }

  write_tsv(df_final2,file="../results/6.Network_analysis//Unassigned_blastn-result-noblank-corrected_blast_97.tsv")
  df_final3 = df_final2[,1:3] %>% rename(`Feature ID` = Feature.ID)
  write_tsv(df_final3,file="../results/6.Network_analysis//taxonomy-corrected-again.tsv")
  return(df_final2)
}

# qiime tools import --input-path ../results/6.Network_analysis/taxonomy-corrected-again.tsv --type 'FeatureData[Taxonomy]' --output-path ../results/6.Network_analysis/taxonomy-corrected



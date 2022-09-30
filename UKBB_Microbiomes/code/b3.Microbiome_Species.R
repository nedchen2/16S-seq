require(tidyverse)

table_genus = read.csv("../results/5.Enterotype/genus_abundance_data.csv",row.names = 1, check.names = F)
metadata <- read.csv("../code/metadata_Microbiome.csv") %>% mutate(SampleID = as.character(SampleID)) %>% dplyr::select(SampleID, Species)

relative.feature.table <- apply(table_genus, 2, function(x) x/sum(x))

#convert relative_abundance into binary data


relative.feature.table[relative.feature.table > 0] = 1
data_prevalence <- as.data.frame(relative.feature.table) %>% rownames_to_column("Genus") %>%
  pivot_longer(cols = colnames(relative.feature.table),names_to = "SampleID",values_to = "relative_abundance") %>% left_join(metadata,by = "SampleID") %>% 
  group_by(Genus,Species) %>% summarise(prevalence = mean(relative_abundance))


relative.feature.table <- apply(table_genus, 2, function(x) x/sum(x))

substr_relative_genus_abundances_heatmap <- as.data.frame(relative.feature.table) %>% rownames_to_column("Genus") %>%
  pivot_longer(cols = colnames(relative.feature.table),names_to = "SampleID",values_to = "relative_abundance") %>% left_join(metadata,by = "SampleID") %>% 
  group_by(Genus,Species) %>% summarise(mean_abundance = mean(relative_abundance)) %>% left_join(data_prevalence) 

#substr_relative_genus_abundances_heatmap <- as.data.frame(relative.feature.table) %>% rownames_to_column("Genus") %>%
#  pivot_longer(cols = colnames(relative.feature.table),names_to = "SampleID",values_to = "relative_abundance") %>% left_join(metadata,by = "SampleID") %>% 
#  group_by(Genus,Species) %>% summarise(mean_abundance = mean(relative_abundance)) %>% pivot_wider(names_from = Species,values_from = mean_abundance) %>% column_to_rownames("Genus")

write.csv(substr_relative_genus_abundances_heatmap,"../results/7.Final_graph/Genus_level_Species_mean_relative_abundance.csv")

data <- substr_relative_genus_abundances_heatmap %>% 
  mutate(Percent_Abundance = round(mean_abundance * 100,2),
         Percent_Prevalence = round(prevalence * 100,2))  %>% mutate(NormAbundance=log10(Percent_Abundance+0.01),
                                                                     NormPrevalence=log10(Percent_Prevalence+0.01))

dataTerr <- data %>% filter(Species == "B.terrestris") %>% arrange(mean_abundance) %>% tail(10)

dataTerr$Genus <- factor(dataTerr$Genus, levels = unique(dataTerr$Genus))

dataTerr <- dataTerr %>% pivot_longer(cols = c("Percent_Abundance", "Percent_Prevalence"),names_to = "Stat",values_to = "abundance")



dataRud <- data %>% filter(Species == "B.ruderatus") %>% arrange(mean_abundance) %>% tail(10)

dataRud$Genus <- factor(dataRud$Genus, levels = unique(dataRud$Genus))

dataRud <- dataRud %>% pivot_longer(cols = c("Percent_Abundance", "Percent_Prevalence"),names_to = "Stat",values_to = "abundance")


dataHor <- data %>% filter(Species == "B.hortorum") %>% arrange(mean_abundance) %>% tail(10)

dataHor$Genus <- factor(dataHor$Genus, levels = unique(dataHor$Genus))

dataHor<- dataHor %>% pivot_longer(cols = c("Percent_Abundance", "Percent_Prevalence"),names_to = "Stat",values_to = "abundance")


plot_heatmap_2 <- function(dataAbundance){
  p <- ggplot(data = dataAbundance, aes(x= Stat, y=Genus, fill= log10(abundance+0.01))) +
    geom_tile(color = "white",lwd = 1.5,linetype = 1) +
    coord_fixed() + 
    geom_text(aes(label =  abundance),color = "white",size = 3) + 
    theme(axis.text.x=element_text(size = 0,face = "bold",angle = 90,vjust = 0.2)) + 
    scale_fill_gradient(low = "navy",high = "firebrick",limits=c(-1,2.1)) +
    labs( y = "", x = "") +
    guides(fill= guide_colorbar(title="log10(% Abundance/Prevalence)")) + theme(legend.position = "left")
  
  return(p)
}


p = plot_heatmap_2(dataAbundance = dataTerr)
ggsave(p,filename = "../results/7.Final_graph/Ter_MeanAbundance_Heatmap.png",width = 2, height = 4,dpi = "print")


p = plot_heatmap_2(dataAbundance = dataHor)
ggsave(p,filename = "../results/7.Final_graph/Hor_MeanAbundance_Heatmap.png",width = 2, height = 4,dpi = "print")


p = plot_heatmap_2(dataAbundance = dataRud)
ggsave(p,filename = "../results/7.Final_graph/Rud_MeanAbundance_Heatmap.png",width = 3, height = 4,dpi = "print")


# descriptive alpha_diversity 

alpha_matrix <- read.csv("../results/4.Diversity_ana/Alpha_Diversity_Matrix.csv") %>% mutate(SampleID = as.character(sample_name)) %>% left_join(metadata) 

terr_matrix <- alpha_matrix %>% filter(Species == "B.terrestris") 

rud_matrix <- alpha_matrix %>% filter(Species == "B.ruderatus" )

Hor_matrix <- alpha_matrix %>% filter(Species == "B.hortorum" )

library(pastecs)
stat.desc(terr_matrix$faith_pd)
stat.desc(rud_matrix$faith_pd)
stat.desc(Hor_matrix$faith_pd)



# =================== SIMPER ANALYSIS

library("vegan")

test_simper = as.data.frame(t(table_genus))

metadata <- read.csv("../code/metadata_Microbiome.csv") %>% mutate(SampleID = as.character(SampleID)) 

metadata2 <- metadata[colnames(table_genus)%in%metadata$SampleID,] %>% column_to_rownames("SampleID")

metadata2 <- metadata2[row.names(test_simper),] %>% rownames_to_column("SampleID")

test_simper <- test_simper %>% rownames_to_column("SampleID") %>% dplyr::select(!"SampleID")

model = simper(test_simper,metadata2$Species,permutations = 999)

B.hortorum_B.terrestris <- model$B.hortorum_B.terrestris

table_genus 

list_of_genus = c("Snodgrassella" ,"Lactobacillus" ,"Arsenophonus" ,"Pseudomonas", "Apibacter", "Fructobacillus", "Apibacter" ,"Pantoea"  ,"Dubosiella" , "Gilliamella" )

# ========== rename genus 

DATA <- as.data.frame(table_genus) %>% rownames_to_column("Genus") %>% mutate(OTU = paste0("OTU",rownames(.)))

Taxon = DATA %>% dplyr::select(OTU,Genus)

OTU = DATA %>% dplyr::select(-Genus) %>% column_to_rownames("OTU") %>% t() %>% as.data.frame()

metadata3 <- metadata %>% dplyr::select("Species")

simper.pretty = function(x, metrics, interesting, perc_cutoff, low_cutoff, low_val, output_name){
  library(vegan)
  interesting = c("Species")
  
  metrics = metadata3
  
  x = test_simper
  
  output_name = "Test"
  perc_cutoff = 0.8
  
  for(variables in interesting){
    test_1=with(metrics, simper(x, metrics[[variables]]))
    #parsing through simper output, saving desired info to table
    for(name in names(test_1)){
      testmx=matrix(ncol=length(interesting))
      testmx=cbind(test_1[[name]]$ord,test_1[[name]]$cusum)
      sorted=testmx[order(testmx[,1]),]
      sorted=cbind(sorted,test_1[[name]]$species)
      sorted=sorted[order(sorted[,2]),]
      t=matrix(sorted[sorted[,2]<=perc_cutoff,],ncol=3)
      i=nrow(t)
      #converting percents to percent of whole
      while(i>1){
        t[i,2]=as.character(as.numeric(t[i,2])-as.numeric(t[i-1,2]))
        i=i-1
      }
      t[,1]=name
      write.table(t,file=paste(output_name,'_simper.csv',sep=""), append=TRUE, sep=",", col.names = FALSE)
       }
    }
  y=read.table(paste(output_name,'_simper.csv',sep=""), header=FALSE,sep=",",fill = TRUE,row.names = NULL)
  file.remove(paste(output_name,'_simper.csv',sep = ""))
  y=y[-c(1)]
  colnames(y) = c("Comparison", "SIMPER", "OTU")

  #prevents changing of colnames if OTU table

  write.csv(y,file=paste(output_name,'_clean_simper.csv', sep=''))
}

#Using the function
simper.pretty(test_simper, metadata3, c('Species'), perc_cutoff=0.5, low_cutoff = 'y', low_val=0.01, 'name')

y=read.table('./Test_clean_simper.csv', header=TRUE,sep=",",fill = TRUE,row.names = NULL)

list <- unique(y$OTU)

relative_genus_table <- t(relative.feature.table) %>% as.data.frame()

metadata2 <- metadata2 %>% dplyr::select(SampleID, Species)

select_taxa <- relative_genus_table[colnames(relative_genus_table) %in% list] %>% rownames_to_column("SampleID") %>% left_join(metadata2,by="SampleID")

krusk = c()
d.f = c()
stat = c()
Bh_relative_mean = c()
Bh_sd = c()
Br_relative_mean = c()
Br_sd = c()
Bt_relative_mean = c()
Bt_sd = c()
tax=c()

select_taxa = subset(select_taxa, select_taxa$Species != "B.terrestris")
select_taxa = subset(select_taxa, select_taxa$Species != "B.hortorum")
select_taxa = subset(select_taxa, select_taxa$Species != "B.ruderatus")

for (otus in unique(y$OTU)){
  print (otus)
  result = kruskal.test(select_taxa[[otus]], select_taxa$Species)
  
  
    
  krusk=append(krusk, result$p.value)
  d.f = append(d.f, result$parameter)
  stat =  append(stat, result$statistic)
  fdr=p.adjust(krusk, method='fdr')
  df_Bh = subset(select_taxa, select_taxa$Species == "B.hortorum")
  Bh_relative_mean = append(Bh_relative_mean, mean(df_Bh[[otus]]))
  Bh_sd = append(Bh_sd, sd(df_Bh[[otus]]))
  
  df_Br =  subset(select_taxa, select_taxa$Species == "B.ruderatus")
  Br_relative_mean = append(Br_relative_mean, mean(df_Br[[otus]]))
  Br_sd = append(Br_sd, sd(df_Br[[otus]]))
  
  df_Bt =  subset(select_taxa, select_taxa$Species == "B.terrestris")
  Bt_relative_mean = append(Bt_relative_mean, mean(df_Bt[[otus]]))
  Bt_sd = append(Bt_sd, sd(df_Bt[[otus]]))
  
  tax =  append(tax, otus)
}



o_csv = matrix(nrow = 10,ncol = 10)

o_csv[,1]=round(Bh_relative_mean*100,2)
o_csv[,2]=round(Bh_sd*100,2)
o_csv[,3]=round(Br_relative_mean*100,2)
o_csv[,4]=round(Br_sd*100,2)
o_csv[,5]=round(Bt_relative_mean*100,2)
o_csv[,6]=round(Bt_sd*100,2)
o_csv[,7]=round(fdr,3)
o_csv[,8]=d.f
o_csv[,9]=round(stat,2)
o_csv[,10]=tax

colnames(o_csv) = c("Bh_relative_mean_abundance","Bh_sd","Br_relative_mean_abundance","Br_sd","Bt_relative_mean_abundance","Bt_sd","fdr_krusk_p.val","d.f","stat","genus")

write.csv(o_csv,file = "../results/7.Final_graph/SIMPER_tax_Bh_vs_Br.csv",row.names = F)
write.csv(o_csv,file = "../results/7.Final_graph/SIMPER_tax_Bt_vs_Br.csv",row.names = F)
write.csv(o_csv,file = "../results/7.Final_graph/SIMPER_tax_Bh_vs_Bt.csv",row.names = F)


# ====== Apibacter box plot
Crithidia_binomial <- metadata %>% select("SampleID","Crithidia_binomial")
Nosema_binomial <-  metadata %>% select("SampleID","Nosema_binomial")
  
df = select_taxa %>% left_join(Crithidia_binomial) %>% left_join(Nosema_binomial) %>% mutate(CrithidiaInfection = ifelse(Crithidia_binomial == 1, "infected", "not infected"),
                                                                                             NosemaInfection = ifelse(Nosema_binomial == 1, "infected", "not infected"), 
                                                                                             ApibacterPrevalence = ifelse(Apibacter == 0, 0, 1))

p = ggplot(df, aes(x=CrithidiaInfection, y=log10(Apibacter*100), color=CrithidiaInfection)) +
  geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent")  + theme_classic() + scale_fill_gdocs() +
  #geom_text(data=df, aes(x=CrithidiaInfection, y=y, color=group, label=stat)) + 
  labs(x="", y="Apibacter rel.abu (log 10%)", color="Crithidia") +
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+
  theme(text=element_text(family="sans", size=20)) + 
  scale_color_manual(values = c("navy","firebrick","gold3","grey")) +
  theme(legend.position = c(0.2,0.8))

df %>% group_by(CrithidiaInfection) %>% summarise(apibacter_count = mean(ApibacterPrevalence))
df %>% group_by(NosemaInfection) %>% summarise(apibacter_count = mean(ApibacterPrevalence))


ggsave(p, filename= "../results/7.Final_graph/Boxplot_Apibacter_Crithidia.png", height = 5, dpi = "print")

p = ggplot(df, aes(x=NosemaInfection, y=log10(Apibacter*100), color=NosemaInfection)) +
  geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent")  + theme_classic() + scale_fill_gdocs() +
  #geom_text(data=df, aes(x=CrithidiaInfection, y=y, color=group, label=stat)) + 
  labs(x="", y="Apibacter rel.abu(log10 %)", color="Nosema") +
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+
  theme(text=element_text(family="sans", size=20)) + 
  scale_color_manual(values = c("navy","firebrick","gold3","grey")) +
  theme(legend.position = c(0.2,0.8))

ggsave(p, filename= "../results/7.Final_graph/Boxplot_Apibacter_Nosema.png", height = 5, dpi = "print")

















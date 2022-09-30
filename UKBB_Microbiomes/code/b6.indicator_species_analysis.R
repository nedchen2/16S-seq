library(indicspecies)
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}

genus_abundances <- read.csv("../results/5.Enterotype/genus_abundance_data.csv",row.names = 1,check.names = F) 

relative_genus_abundances <- apply(genus_abundances, 2, function(x) x/sum(x))

data.denoized=noise.removal(relative_genus_abundances, percent=0.01) * 100

relative_genus_abundances_filtered <- data.denoized %>% t() %>% as.data.frame() %>% 
  rownames_to_column("SampleID")

metadata <- read.csv("./metadata_Microbiome.csv") 
groups_species <- metadata %>% dplyr::select(SampleID,Species,Enterotype,Crithidia_binomial,Nosema_binomial) %>% mutate(SampleID = as.character(SampleID))

MergeF <- left_join(groups_species,relative_genus_abundances_filtered )

abund = MergeF[,6:ncol(MergeF)]
Enterotypes = MergeF$Enterotype
Species = MergeF$Species
Crithidia_infected <- MergeF$Crithidia_binomial
Nosema_infected  <-  MergeF$Nosema_binomial

inv = multipatt(abund, Enterotypes,func =  "r.g", control = how(nperm=9999))
summary(inv,indvalcomp=TRUE,alpha = 0.05)

inv2 = multipatt(abund, Species, func = "r.g",  control = how(nperm=9999))
summary(inv2,indvalcomp=TRUE,alpha = 0.05)

inv2 = multipatt(abund, Species,   control = how(nperm=9999),duleg = T)
summary(inv2,indvalcomp=TRUE,alpha = 0.05)


inv3 = multipatt(abund, Crithidia_infected,func = "r.g", control = how(nperm=9999))
summary(inv3,indvalcomp=TRUE,alpha = 0.05)

inv3 = multipatt(abund, Crithidia_infected,func = "r.g.", control = how(nperm=9999))
summary(inv3,indvalcomp=TRUE,alpha = 0.05)

inv4 = multipatt(abund, Nosema_infected, func = "r.g",  control = how(nperm=9999))
summary(inv4,indvalcomp=TRUE,alpha = 0.05)

# ============= Expand Apibacter

MergeF %>% ggplot(aes(x = as.factor(Crithidia_binomial), y = Apibacter)) + geom_boxplot()


MODEL <- aov(Apibacter~Crithidia_binomial,data=MergeF)
summary(MODEL)

MergeF %>% ggplot(aes(x = as.factor(Nosema_binomial), y = Apibacter)) + geom_boxplot()

MODEL <- aov(Apibacter~Nosema_binomial,data=MergeF)
summary(MODEL)

# richness of Apibacter

MergeF %>% mutate(Apibacter_binomial = ifelse(Apibacter == 0, 0 , 1)) %>% group_by(Nosema_binomial) %>% summarise(richness = sum(Apibacter_binomial),
                                                                                                                  Richness = n(),
                                                                                                                  Percentage= richness/Richness)

MergeF %>% mutate(Apibacter_binomial = ifelse(Apibacter == 0, 0 , 1)) %>% group_by(Crithidia_binomial) %>% summarise(richness = sum(Apibacter_binomial),
                                                                                                                  Richness = n(),
                                                                                                                  Percentage= richness/Richness)
data <- MergeF %>% mutate(Apibacter_binomial = ifelse(Apibacter == 0, 0 , 1))

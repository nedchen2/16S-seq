require(tidyverse)
require(qiime2R)

# load the parasite count data 
# load the infection and Environmental factor

metadata <- read.csv("../../UKBB_Microbiomes//code/UKBB_NED_filtered.csv") %>% dplyr::rename(SampleID =`Bee_ID` ,
                                                                                             Species = `finalised_species`,
                                                                                             CollectionSite = `Site`)  %>% mutate(interaction = paste0(Species, "_",CollectionSite)) %>%
  filter(parasite_counts_available == "YES" | RADseq_available == "YES")


TABLE <- t(table_genus) %>% as.data.frame() %>% rownames_to_column("SampleID")

newmetadata <- metadata%>% left_join(TABLE)  %>% left_join(metadata_alpha_diversity)

newmetadata %>% write_csv('../data/UKBB_Metadata.csv')

# ------- prevalence and inbreeding coefficient and alpha microbial diversity

metadata2 <- read.csv("../../UKBB_Microbiomes//code/UKBB_NED_filtered.csv") %>% dplyr::rename(SampleID =`Bee_ID` ,
                                                                                              Species = `finalised_species`,
                                                                                              CollectionSite = `Site`)  %>% 
  mutate(interaction = paste0(Species, "_",CollectionSite)) %>% filter(!is.na(Crithidia_binomial))


id <- read.table("../../UKBB_RADseq/results/d.population/distance/plink.mdist.id") %>%  mutate(sample_name =  str_extract(V1,"(?<=PG_)(.*)(?=_1.sorted?)")) %>% pull(sample_name)



metadata2 <- metadata2 %>% filter(SampleID %in% id)

cal_mean_non_zero <- function(data){
  idx = data != 0
  return(mean(data[idx],na.rm=T))
  
}



mean_parasite_load <- metadata2 %>% group_by(interaction) %>% summarize(meanCrithidia = mean(Crithidia_total_per_bee,na.rm = T),
                                                                      meanNosema = mean(Nosema_total_per_bee,na.rm = T),
                                                                      meanApicystis = mean(Apicystis_total_per_bee,na.rm = T),
                                                                      meanCrithidiaNon.Zero = cal_mean_non_zero(Crithidia_total_per_bee),
                                                                      meanNosemaNon.Zero =  cal_mean_non_zero(Nosema_total_per_bee),
                                                                      meanApicystisNon.Zero = cal_mean_non_zero(Apicystis_total_per_bee)) 





genetic_df <- read.delim("../../UKBB_RADseq/results/c.ref_map/test.txt") 

# The result of alpha diversity
metadata_alpha_diversity <- read.csv("../../UKBB_Microbiomes/results/4.Diversity_ana/Alpha_Diversity_Matrix.csv")  %>% mutate(SampleID = as.character(sample_name))

test_alpha_diversiy = metadata_alpha_diversity %>% left_join(metadata2) %>% mutate(Ruderatus = ifelse(Species == "B.ruderatus", "Yes","No"))

kruskal.test(test_alpha_diversiy$faith_pd,test_alpha_diversiy$Ruderatus)



metadata2 <- metadata2 %>% left_join(metadata_alpha_diversity)


# ==== if the inbreeding coefficient within different species differ


Crithidia_stack_plot <- metadata2 %>% 
  group_by(interaction,Species) %>% 
  summarise(Parasite.present = sum(Crithidia_binomial),frequency = n(), FaithAlphaDiversity = mean(faith_pd,na.rm =  T)) %>%
  mutate(CrithidiaIncidence =Parasite.present/frequency,
         no.Parasite.present = frequency - Parasite.present, 
         total = sum(frequency),
         parasite = "Crithidia"
        ) %>% left_join(genetic_df)




model = aov(formula = Fis ~ Species ,data = Crithidia_stack_plot)
summary(model)

TukeyHSD(model)


ggplot(data=Crithidia_stack_plot, aes(x=Fis,y = Percent)) + geom_point() 

hist(scale(genetic_df$Fis))
hist(scale(genetic_df$Pi))

summary(lm(CrithidiaIncidence~Pi,data = Crithidia_stack_plot))



Nosema_stack_plot <- metadata2 %>% 
  group_by(interaction,Species) %>% 
  summarise(Parasite.present = sum(Nosema_binomial),frequency = n()) %>%
  mutate(NosemaIncidence =Parasite.present/frequency,
         no.Parasite.present = frequency - Parasite.present, 
         total = sum(frequency),
         parasite = "Nosema") %>% left_join(genetic_df)

Apicystis_stack_plot <- metadata2 %>% 
  group_by(interaction,Species) %>% 
  summarise(Parasite.present = sum(Apicystis_binomial),frequency = n()) %>%
  mutate(ApicystisIncidence =Parasite.present/frequency,
         no.Parasite.present = frequency - Parasite.present, 
         total = sum(frequency),
         parasite = "Apicystis") %>% left_join(genetic_df) %>% dplyr::select(interaction,ApicystisIncidence)

info <- metadata2 %>% select(interaction,Species,CollectionSite)

# ======== important 
dataForPairs = Nosema_stack_plot %>% dplyr::select(interaction,NosemaIncidence) %>% left_join(Crithidia_stack_plot) %>% left_join(Apicystis_stack_plot)%>% 
  dplyr::select(interaction,CrithidiaIncidence,NosemaIncidence,ApicystisIncidence,Fis,Pi, FaithAlphaDiversity)  %>% mutate(Fis = Fis,Pi = Pi) %>% left_join(mean_parasite_load) %>% 
  dplyr::select(-meanCrithidia,-meanNosema,-meanApicystis) %>% 
  tidyr::replace_na(list(meanNosemaNon.Zero = 0,meanApicystisNon.Zero =0)) %>% dplyr::rename(NosemaCount = meanNosemaNon.Zero,
                                                                                      ApicystisCount = meanApicystisNon.Zero,
                                                                                      CrithidiaCount = meanCrithidiaNon.Zero) %>% left_join(info, by = "interaction") %>% distinct()


cor.test(dataForPairs$FaithAlphaDiversity,dataForPairs$Fis)
cor.test(dataForPairs$CrithidiaIncidence,dataForPairs$Fis,method = "kendall")
cor.test(dataForPairs$NosemaIncidence,dataForPairs$Fis,method = "kendall")
cor.test(dataForPairs$ApicystisIncidence,dataForPairs$Fis,method = "spearman")
cor.test(dataForPairs$ApicystisIncidence,dataForPairs$Fis,method = "spearman")
cor.test(dataForPairs$NosemaIncidence,dataForPairs$Fis,method = "spearman")

dataForPairs

model_Crithidia = lmer(formula = CrithidiaIncidence ~ Fis + (1|CollectionSite)  + (1|Species)  , data = dataForPairs)
summary(lmer(formula = ApicystisIncidence ~ Fis + (1|CollectionSite) + (1|Species) , data = dataForPairs))

model_Nosema = lmer(formula = NosemaIncidence ~ Fis + (1|CollectionSite) + (1|Species)  , data = dataForPairs)
summary(lmer(formula = CrithidiaIncidence ~ Fis + (1|CollectionSite) + (1|Species) , data = dataForPairs))

model_Apcystis = lmer(formula = ApicystisIncidence ~ Fis + (1|CollectionSite) , data = dataForPairs)
summary(lmer(formula = NosemaIncidence ~ Fis + (1|CollectionSite) + (1|Species) , data = dataForPairs))


summary(lmer(formula = FaithAlphaDiversity  ~ Fis + (1|CollectionSite) + (1|Species) , data = dataForPairs))

require(sjPlot)

tab_model(model_Crithidia,model_Nosema,model_Apcystis,collapse.ci = T,show.reflvl = F,show.stat = T,prefix.labels = "varname",show.icc = F,file = "../../UKBB_Microbiomes//results/7.Final_graph/Linearmixed_table_2.html")


cor(dataForPairs$ApicystisIncidence,dataForPairs$Fis)

hist(dataForPairs$Fis)
hist(dataForPairs$ApicystisIncidence)

model3 <- glmer( CrithidiaIncidence ~ Fis + (1|interaction),data = dataForPairs,family = poisson)
summary(model3)




write.csv(dataForPairs, file = "../../UKBB_Microbiomes/data/Population_Prevalence_inbreeding_alpha_diversity.csv" )
  
write.csv(dataForPairs, file = "../../UKBB_Microbiomes/results/7.Final_graph/inbreeding_coefficient_Table.csv" )


summary(lm( formula = Fis ~ scale(FaithAlphaDiversity),data = dataForPairs ))

summary(lm( formula = Fis ~ scale(CrithidiaIncidence),data = dataForPairs ))

summary(lm( formula = Fis ~ scale(NosemaIncidence),data = dataForPairs ))

summary(lm( formula = Fis ~ scale(ApicystisIncidence),data = dataForPairs ))

summary(lm( formula = Fis ~ ApicystisIncidence +  FaithAlphaDiversity + CrithidiaIncidence + NosemaIncidence,data = dataForPairs ))


chart.Correlation.linear <-function (R, histogram = TRUE, method=c("pearson", "kendall", "spearman"), ...)
{ # @author R Development Core Team
  # @author modified by Peter Carl & Marek Lahoda
  # Visualization of a Correlation Matrix. On top the (absolute) value of the correlation plus the result 
  # of the cor.test as stars. On botttom, the bivariate scatterplots, with a linear regression fit. 
  # On diagonal, the histograms with probability, density and normal density (gaussian) distribution.
  
  x = checkData(R, method="matrix")
  
  if(missing(method)) method=method[1] #only use one
  cormeth <- method
  
  # Published at http://addictedtor.free.fr/graphiques/sources/source_137.R
  panel.cor <- function(x, y, digits=2, prefix="", use="pairwise.complete.obs", method=cormeth, cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use=use, method=method) # MG: remove abs here
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 2.5
    
    test <- cor.test(as.numeric(x),as.numeric(y), method=method)
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " "))
    # MG: add abs here and also include a 30% buffer for small numbers
    text(0.5, 0.5, txt, cex = cex)
    text(.8, .8, Signif, cex=cex, col=2)
  }
  
  #remove method from dotargs
  dotargs <- list(...)
  dotargs$method <- NULL
  rm(method)
  
  hist.panel = function (x, ...=NULL ) {
    par(new = TRUE)
    hist(x,
         col = "light gray",
         probability = TRUE,
         axes = FALSE,
         main = "",
         breaks = "FD")
    lines(density(x, na.rm=TRUE),
          col = "red",
          lwd = 1)
    # adding line representing density of normal distribution with parameters correponding to estimates of mean and standard deviation from the data 
    ax.x = seq(min(x), max(x), 0.1)                                                  # ax.x containts points corresponding to data range on x axis
    density.est = dnorm(ax.x, mean = mean(x), sd = sd(x))   # density corresponding to points stored in vector ax.x 
    lines(ax.x, density.est, col = "blue", lwd = 1, lty = 1)                                # adding line representing density into histogram
    rug(x)
  }
  
  # Linear regression line fit over points
  reg <- function(x, y, ...) {
    points(x,y, cex=0.5, pch = 19,...)
    abline(lm(y~x), col = "red") 
  }
  
  # Draw the chart
  if(histogram)
    pairs(x, gap=0, lower.panel=reg, upper.panel=panel.cor, diag.panel=hist.panel)
  else
    pairs(x, gap=1, lower.panel=reg, upper.panel=panel.cor) 
}


Incidence <- dataForPairs %>% dplyr::select(CrithidiaIncidence,NosemaIncidence,ApicystisIncidence,Fis,Pi)

COUNT <-  dataForPairs %>% dplyr::select(CrithidiaCount,NosemaCount,ApicystisCount,Fis,Pi)

library(outliers)

grubbs.test(COUNT$ApicystisCount)

pdf(file = "../results/7.Final_graph/CorrelationPlot_Crithidia_Nosema_Prevalence.pdf",width = 9,pointsize = 20)
chart.Correlation.linear(Incidence,histogram = T)
dev.off()

pdf(file = "../results/7.Final_graph/CorrelationPlot_Crithidia_Nosema_Count.pdf",width = 9)
chart.Correlation.linear(COUNT ,histogram = T)
dev.off()


# ============================================================
fixed = metadata %>% 
  dplyr::select(SampleID,CollectionSite,Species,"Crithidia_total_per_bee","Apicystis_total_per_bee","Nosema_total_per_bee","gut_parasite_richness","Crithidia_binomial","Nosema_binomial","Apicystis_binomial")

# load the diversity matrix 
faith_pd_vector<-read_qza("../results/4.Diversity_ana/core-metrics-results/faith_pd_vector.qza")$data %>% rownames_to_column("SampleID") 
shannon_vector <- read_qza("../results/4.Diversity_ana/core-metrics-results/shannon_vector.qza")$data %>% rownames_to_column("SampleID") 
evenness_vactor <- read_qza("../results/4.Diversity_ana/core-metrics-results/evenness_vector.qza")$data %>% rownames_to_column("SampleID") 

alpha_matrix = as.data.frame(left_join(faith_pd_vector,shannon_vector))  %>% left_join(evenness_vactor)
str(alpha_matrix)



# bacterial loads
list_of_biomarker = c( "Snodgrassella","Bifidobacterium", "Apibacter", "Gilliamella","Lactobacillus")

#bacterial_relative_phylum <-  apply(merged_table2, 2, function(x) x *100/sum(x)) %>% t() %>% as.data.frame() %>% rownames_to_column("SampleID") 

bacterial_relative <-  read.csv("../results/5.Enterotype/genus_abundance_data.csv",check.names = F,row.names = 1)

bacterial_relative <- apply(bacterial_relative, 2, function(x) x *100/sum(x)) %>% t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>%
  dplyr::select(list_of_biomarker,SampleID)


bacterial <- read.csv("../results/5.Enterotype/genus_abundance_data.csv",check.names = F,row.names = 1) %>% t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>%
  dplyr::select(list_of_biomarker,SampleID)




# load the Hostgenetics 

F_list <- read.csv("../../UKBB_RADseq/code/metadata_RADseq.csv") %>% dplyr::select(index_list,F_list) %>% 
  mutate(SampleID = as.character(index_list))
str(F_list)



map_list =  read.csv(file = "../../UKBB_RADseq//genome/metadata_individuals_2.csv")  %>% 
  dplyr::select(name_list,site_list, field_ID_list,F_list,index_list) %>% 
  mutate(name_list = paste0(name_list,".sorted")) %>% 
  mutate(Pop.ID = paste0(field_ID_list,"_",site_list)) %>% mutate(SampleID = as.character(index_list)) %>% dplyr::select(-index_list)

Final_table <- left_join(map_list, genetic_df, by = "Pop.ID") 

Genetic_table <- Final_table # n = 206
# combine four 

DataForLM <- left_join(fixed,F_list) %>% left_join(alpha_matrix) %>% left_join(bacterial_relative, by = "SampleID")

DataForLM <- left_join(fixed,F_list) %>% left_join(alpha_matrix) %>% left_join(bacterial, by = "SampleID")

#DataForLM <- left_join(fixed,F_list) %>% left_join(alpha_matrix) %>% left_join(bacterial_relative_phylum, by = "SampleID")

require(lme4)
require(sjPlot)
require(sjmisc)
require(sjlabelled)
require(ggpubr)
require(pscl)
require(lmerTest)
require(performance)

# ========= explore parasite count data 

Crithidia_lmer <- DataForLM %>% filter(Crithidia_total_per_bee != 0, !is.na(DataForLM$faith_pd))
Nosema_lmer <- DataForLM %>% filter(Nosema_total_per_bee != 0, !is.na(DataForLM$faith_pd))
Apicystis <- DataForLM %>% filter(Apicystis_total_per_bee != 0, !is.na(DataForLM$faith_pd))




model1 <- lmer(log10(Crithidia_total_per_bee) ~ log10(Snodgrassella+1) + log10(Bifidobacterium+1) + log10(Apibacter+1) + log10(Gilliamella+1) + log10(Lactobacillus+1)  + (1|Species), data=Crithidia_lmer)
summary(model1)

model2 <- lmer(log10(Nosema_total_per_bee) ~ log10(Snodgrassella+1) + log10(Bifidobacterium+1) + log10(Apibacter+1) + log10(Gilliamella+1) + log10(Lactobacillus+1)  + (1|Species), data=Nosema_lmer)
summary(model2)

tab_model(model1,model2,collapse.ci = T,show.reflvl = F,show.stat = T,prefix.labels = "varname",show.icc = F,file = "../results/7.Final_graph/Linearmixed_table.html")

tab_model(model2,collapse.ci = T,show.reflvl = F,show.stat = T,prefix.labels = "varname",show.icc = F,file = "../results/7.Final_graph/Linearmixed_table2.html")
#model <- lmer(Apicystis_total_per_bee ~  Snodgrassella + Bifidobacterium + Apibacter + Gilliamella + Lactobacillus + (1|Species) , data=Apicystis)
#summary(model)


model3 <- glmer(Crithidia_binomial ~ log10(Snodgrassella+1) + log10(Bifidobacterium+1) + log10(Apibacter+1) + log10(Gilliamella+1) + log10(Lactobacillus+1)  + (1|CollectionSite),data = DataForLM,family = binomial)

summary(model3)

model3 <- glmer(Crithidia_binomial ~ Snodgrassella + Bifidobacterium + Apibacter + Gilliamella + Lactobacillus +  (1|SampleID) , data = DataForLM,family = binomial)

summary(model3)


# ================

p = ggplot(metadata, aes(Crithidia_total_per_bee )) + geom_histogram()

ggsave(p,filename = "../results/7.Final_graph/distribution/cithidia_no.png",dpi = "print")

p = ggplot(metadata, aes(log10(Crithidia_total_per_bee + 1 ))) + geom_histogram()

ggsave(p,filename = "../results/7.Final_graph/distribution/crithidia_yes.png",dpi = "print")

p = ggplot(metadata, aes(Nosema_total_per_bee)) + geom_histogram()

ggsave(p,filename = "../results/7.Final_graph/distribution/Nosema_no.png",dpi = "print")

p = ggplot(metadata, aes(log10(Nosema_total_per_bee + 1))) + geom_histogram()

ggsave(p,filename = "../results/7.Final_graph/distribution/Nosema_yes.png",dpi = "print")

p = ggplot(metadata, aes(Apicystis_total_per_bee)) + geom_histogram()

ggsave(p,filename = "../results/7.Final_graph/distribution/Api_yes.png",dpi = "print")

p = ggplot(metadata, aes(log10(Apicystis_total_per_bee + 1))) + geom_histogram()

ggsave(p,filename = "../results/7.Final_graph/distribution/Api_no.png",dpi = "print")



# ===========================BINOMIAL ===============================

Crithidia_binomial2 <- DataForLM %>% mutate(CrithidiaInfection = (Crithidia_total_per_bee == 0)) %>% filter(!is.na(shannon_entropy)) %>% mutate(Priority = ifelse(Species == "B.ruderatus","Yes","No"))

kruskal.test(Crithidia_binomial2$faith_pd,Crithidia_binomial2$Priority)

model = glm(formula = CrithidiaInfection ~ Firmicutes   ,data = Crithidia_binomial2,family =  "binomial")

summary(model)

model = glm(formula = CrithidiaInfection ~ Lactobacillus  + Species   ,data = Crithidia_binomial2,family =  "binomial")

summary(model)


model = glm(formula = CrithidiaInfection ~ Bifidobacterium + Species   ,data = Crithidia_binomial2,family =  "binomial")

summary(model)

model = glm(formula = CrithidiaInfection ~ Apibacter   ,data = Crithidia_binomial2,family =  "binomial")

summary(model)

model = glm(formula = CrithidiaInfection ~ Snodgrassella ,data = Crithidia_binomial2,,family =  "binomial")

summary(model)

model = glm(formula = CrithidiaInfection ~ Gilliamella ,data = Crithidia_binomial2,,family =  "binomial")

summary(model)


model = glm(formula = CrithidiaInfection ~ CollectionSite + Species, data = Crithidia_binomial2, family = "binomial")


logistic_reg2 = glm(formula = CrithidiaInfection ~ Gilliamella + Bifidobacterium + Snodgrassella  + Lactobacillus + CollectionSite + Species,data = Crithidia_binomial2, family = "binomial")

summary(logistic_reg2)


ggplot(Crithidia_binomial2, aes( x = CrithidiaInfection, y = Snodgrassella)) +
  geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x= "", color =  group1) + theme_classic() + scale_fill_gdocs() + facet_grid(~ Species, scales = "free_x",switch = "x") +
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+ 
  theme(text=element_text(family="sans", size=20)) + theme(axis.text = element_text(size = 20)) + theme(axis.text.y = element_text(angle = 45)) + theme(legend.position = "top") + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","yellow"))


pairwise.wilcox.test(Crithidia_binomial2$Apibacter,Crithidia_binomial2$CrithidiaInfection,p.adjust.method = "BH")

pairwise.wilcox.test(Nosema_binomial2$Snodgrassella,Nosema_binomial2$NosemaInfection,p.adjust.method = "BH")



Crithidia_binomial3 <- Crithidia_binomial2 %>% filter(CrithidiaInfection == T)

p = ggplot(Crithidia_binomial3, aes( x = CollectionSite, y = Snodgrassella,color = CollectionSite)) +
  geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x= "", color =  "CollectionSite") + theme_classic() + scale_fill_gdocs() +  facet_grid(~ Species, scales = "free_x",switch = "x") +
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+ 
  theme(text=element_text(family="sans", size=20)) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.text.y = element_text(angle = 45)) + theme(legend.position = "top") + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","grey"))






Nosema_binomial2 <- DataForLM %>% mutate(NosemaInfection = Nosema_total_per_bee == 0) %>% filter(!is.na(Snodgrassella))


model = glm(formula = NosemaInfection ~ Bifidobacterium    ,data =Nosema_binomial2,family = "binomial")

summary(model)

model = glm(formula = NosemaInfection ~ Lactobacillus   ,data =Nosema_binomial2,family = "binomial")

summary(model)


model = glm(formula = NosemaInfection ~ Apibacter  ,data = Nosema_binomial2,family = "binomial")
summary(model)

model = glm(formula = NosemaInfection ~ CollectionSite + Species  ,data = Nosema_binomial2,family = "binomial")

summary(model)

model = glm(formula = NosemaInfection ~ Snodgrassella ,data = Nosema_binomial2,family = "binomial")

summary(model)

model = glm(formula = NosemaInfection ~ Gilliamella ,data = Nosema_binomial2,family = "binomial")

summary(model)


model = glm(formula = NosemaInfection ~ Snodgrassella + CollectionSite + Species,data = Nosema_binomial2, family = "binomial")

summary(model)


logistic_reg = glm(formula = NosemaInfection ~ Gilliamella + Bifidobacterium + Snodgrassella + Lactobacillus + CollectionSite + Species,data = Nosema_binomial2, family = "binomial")

summary(logistic_reg)



#==================================visualize


require(sjPlot)

p = plot_model(logistic_reg,grid = T, show.values = TRUE, value.offset = .3, vline.color = "grey",title = "",show.intercept = T,) + theme_few() + theme( panel.border = element_rect(colour = "white"))

ggsave(p,filename = "../../UKBB_Microbiomes//results/7.Final_graph/forestplot_Nosema.png",dpi="print",scale = 0.8, width = 6,height = 8)

p = plot_model(logistic_reg2,grid = T, show.values = TRUE, value.offset = .3, vline.color = "grey",title = "",show.intercept = T) + theme_few() + theme( panel.border = element_rect(colour = "white"))

ggsave(p,filename = "../../UKBB_Microbiomes//results/7.Final_graph/forestplot_Crithidia.png",dpi="print",scale = 0.8, width = 6,height = 8)











#============ boxplot


Nosema_binomial3 <- Nosema_binomial2 %>% filter(NosemaInfection == 1) 
Nosema_binomial4 <- Nosema_binomial2 %>% filter(NosemaInfection == 0) 

p = ggplot(Nosema_binomial3, aes( x = CollectionSite, y = Snodgrassella,color = CollectionSite)) +
  geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x= "", color =  "CollectionSite") + theme_classic() + scale_fill_gdocs() +  facet_grid(~ Species, scales = "free_x",switch = "x") +
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+ 
  theme(text=element_text(family="sans", size=20)) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.text.y = element_text(angle = 45)) + theme(legend.position = "top") + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","grey"))

ggsave(p, filename = "../../UKBB_Microbiomes/results/7.Final_graph/Nosema_infected_CollectionSite_Species.png",width = 12, height = 6,dpi = "print")

# --- uninfected
p = ggplot(Nosema_binomial4, aes( x = CollectionSite, y = Snodgrassella,color = CollectionSite)) +
  geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x= "", color =  "CollectionSite") + theme_classic() + scale_fill_gdocs() +  facet_grid(~ Species, scales = "free_x",switch = "x") +
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+ 
  theme(text=element_text(family="sans", size=20)) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.text.y = element_text(angle = 45)) + theme(legend.position = "top") + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","grey"))

ggsave(p, filename = "../../UKBB_Microbiomes/results/7.Final_graph/Nosema_uninfected_CollectionSite_Species.png",width = 12, height = 6,dpi = "print")


# ===== infected apibacter

p = ggplot(Nosema_binomial3, aes( x = CollectionSite, y = Apibacter,color = CollectionSite)) +
  geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x= "", color =  "CollectionSite") + theme_classic() + scale_fill_gdocs() +  facet_grid(~ Species, scales = "free_x",switch = "x") +
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+ 
  theme(text=element_text(family="sans", size=20)) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.text.y = element_text(angle = 45)) + theme(legend.position = "top") + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","grey"))

ggsave(p, filename = "../../UKBB_Microbiomes/results/7.Final_graph/Nosema_uninfected_CollectionSite_Species.png",width = 12, height = 6,dpi = "print")

# --- uninfected
p = ggplot(Nosema_binomial4, aes( x = CollectionSite, y = Apibacter,color = CollectionSite)) +
  geom_boxplot(alpha=1, outlier.shape = NA, outlier.size=0, size=0.7, width=0.5, fill="transparent") +
  labs(x= "", color =  "CollectionSite") + theme_classic() + scale_fill_gdocs() +  facet_grid(~ Species, scales = "free_x",switch = "x") +
  geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)+ 
  theme(text=element_text(family="sans", size=20)) + theme(axis.text = element_text(size = 20)) + 
  theme(axis.text.y = element_text(angle = 45)) + theme(legend.position = "top") + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","grey"))

ggsave(p, filename = "../../UKBB_Microbiomes/results/7.Final_graph/Nosema_uninfected_CollectionSite_Species.png",width = 12, height = 6,dpi = "print")






# test for bacterial loads

plot_lm_scatter <- function(data = lm_genus, response = "Crithidia_total_per_bee",explantory = "shannon_entropy",Group = "Enterotype", confidence_interval = T,group = F,add = "Name") {
  library(ggpubr)
  library(ggplot2)
  library(ggthemes)
  
  df_subset <- subset(data, select=c(response,explantory,Group))
  
  
  #log10 transform to response variable
  df_subset[,response] <- log10(df_subset[,response] + 1) 
 # df_subset[,explantory] <- log10(df_subset[,explantory] + 1) 
  colnames(df_subset) <- c("response","explantory","group")
  
  
  #replace the "-/Inf"
  if( "Inf" %in% df_subset[,"response"] | "-Inf" %in% df_subset[,"response"]){
    tmp <- df_subset[which(df_subset[,"response"] != "Inf" & df_subset[,"response"] != "-Inf"),]
    max = max(tmp[,"response"])
    min = min(tmp[,"response"])
    df_subset[,"response"] <-  as.numeric(sub("-Inf", 0, df_subset[,"response"]))
    df_subset[,"response"] <- as.numeric(sub("Inf", max, df_subset[,"response"]))
  }
  
  
  if( "Inf" %in% df_subset[,"explantory"] | "-Inf" %in% df_subset[,"explantory"]){
    tmp <- df_subset[which(df_subset[,"explantory"] != "Inf" & df_subset[,"explantory"] != "-Inf"),]
    max = max(tmp[,"explantory"])
    min = min(tmp[,"explantory"])
    df_subset[,"explantory"] <- as.numeric(sub("-Inf", min,  df_subset[,"explantory"]))
    df_subset[,"explantory"] <- as.numeric(sub("Inf", max, df_subset[,"explantory"]))
  }
  # remove uninfected samples
  df_subset <- subset(df_subset, response !=0)
  
  fit_model <- try(lm(data = df_subset,response ~  explantory),silent = T)
  if (as.vector(summary(fit_model))[2] != "try-error"){
    
    report_tmp <- as.data.frame(coef(summary(fit_model)))
    report_tmp["Parameter"] <- rownames(report_tmp) #add parameter to one col
    # convert to longer sheet
    report_tmp2 <- report_tmp 
    report_tmp2["Explantory"] <- explantory
    report_tmp2["Response"]   <- response
    report_tmp2["SampleSize"] <- nrow(df_subset)
    report_tmp2["Adjusted_R"] <- summary(fit_model)$adj.r.squared
    report_tmp2["p_value"]    <- pf(summary(fit_model)$fstatistic[1],summary(fit_model)$fstatistic[2],summary(fit_model)$fstatistic[3],lower.tail = F)
  }else{
    cat("========================")
    stop("please check the table")
  }
  
  if (group){
    p = ggplot(df_subset,aes(x=explantory,y=response)) + 
      geom_point(aes(color=group),size = 3) + 
      geom_smooth(data = df_subset, aes(x=explantory,y=response), method = "lm",formula = y~x, se = confidence_interval) +
      stat_cor(data = df_subset, method = "pearson",label.sep = "\n",label.x.npc = 0.01,label.y.npc = 0.1) + 
      theme_clean() + 
      labs(y = bquote(paste(.(response),(Log[10]))),x =  bquote(paste(.(add),.("bacterial load"),(Log[10]))), color = "Species") + scale_color_manual(values = c("navy","firebrick","gold3")) + theme(legend.position = "top")
    
    
  }else{
    p = ggplot(df_subset,aes(x=explantory,y=response)) + 
      geom_point(size = 3) + 
      geom_smooth(method = "lm",formula = y~x, se = confidence_interval) +
      stat_cor(data = df_subset, method = "pearson",label.sep = " ",label.x.npc = 0.8) + 
      theme_clean() + 
      labs(y =  bquote(paste(.(response),(Log[10]))),x = explantory, title = "linear regression") 
  }
  
  return(list(picture = p, report = report_tmp2))
}

AmpiconData <- DataForLM %>% filter(!is.na(DataForLM$Gilliamella)) %>% mutate(bacterial_load4 = rowSums(.[12:15])) %>% mutate(bacterial_load3 = rowSums(.[c(12,13,15)]), 
                                                                                                                              bacterial_load2 = rowSums(.[c(12,15)]),
                                                                                                                              bacterial_load2_2 = rowSums(.[c(12,13)]),
                                                                                                                              bacterial_load2_3 = rowSums(.[c(15,13)]),
                                                                                                                              bacterial_load2_4 = rowSums(.[c(14,13)]))

p = plot_lm_scatter(data =  AmpiconData, response = "Crithidia_total_per_bee", explantory = "Bifidobacterium", Group = "Species", confidence_interval = F,group = T, add = "Bifidobacterium ")$picture

ggsave(p,filename = "../results/7.Final_graph/bacterial_load_bifi_reg1.png",dpi = "print")

p = plot_lm_scatter(data =  AmpiconData, response = "Crithidia_total_per_bee", explantory = "Gilliamella", Group = "Species", confidence_interval = F,group = T, add = "Gilliamella ")$picture

ggsave(p,filename = "../results/7.Final_graph/bacterial_load_Gillia_reg1.png",dpi = "print")

p = plot_lm_scatter(data =  AmpiconData, response = "Crithidia_total_per_bee", explantory = "Apibacter", Group = "Species", confidence_interval = F,group = T, add = "Apibacter ")$picture

ggsave(p,filename = "../results/7.Final_graph/bacterial_load_Apibacter_reg1.png",dpi = "print")

p = plot_lm_scatter(data =  AmpiconData, response = "Crithidia_total_per_bee", explantory = "Snodgrassella", Group = "Species", confidence_interval = F,group = T, add = "Snodgrassella ")$picture

ggsave(p,filename = "../results/7.Final_graph/bacterial_load_Snodgrassella_reg1.png",dpi = "print")

p = plot_lm_scatter(data =  AmpiconData, response = "Crithidia_total_per_bee", explantory = "Snodgrassella", Group = "Species", confidence_interval = F,group = T, add = "Snodgrassella ")$picture

ggsave(p,filename = "../results/7.Final_graph/bacterial_load_Snodgrassella_reg1.png",dpi = "print")


p = plot_lm_scatter(data =  AmpiconData, response = "Crithidia_total_per_bee", explantory = "faith_pd", Group = "Species", confidence_interval = F,group = F, add = "Snodgrassella ")$picture

ggsave(p,filename = "../results/7.Final_graph/faith_pd_reg1.png",dpi = "print")







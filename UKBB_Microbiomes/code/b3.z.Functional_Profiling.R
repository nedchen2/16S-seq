df_function = read.table("../results/3.Taxonomy_ana/q2-picrust2_output/Functional_profiling_ko/pred_metagenome_unstrat_descrip.tsv",sep = "\t",header = T,check.names = F) %>%
  column_to_rownames(var = "function") %>% select(!c(description)) %>%
  mutate(Frequency = rowSums(.)) %>% 
  arrange(desc(Frequency)) %>% select(!c(Frequency)) %>% drop_na()


Create_Table_topN <- function(data = df_function  ,topN=10) {

  other = colSums(data[(topN+1):dim(data)[1],])
  
  df_relative_abundance_top =data[1:(topN),]
  
  df_relative_abundance_top = rbind(df_relative_abundance_top, other)
  
  rownames(df_relative_abundance_top)[topN+1] = c("Other")
  
  df_relative_abundance_top$Function = rownames(df_relative_abundance_top)
  
  data_all = as.data.frame(melt(df_relative_abundance_top, id.vars = c("Function")))
  
  #To let the bar plot sorted with abundance
  data_all$Function  = factor(data_all$Function, levels = rownames(df_relative_abundance_top))
  
  return (data_all)
  
}

data_all <- Create_Table_topN(data = df_function ,topN=10)

colors = brewer.pal(11,"Paired")



plot_stack_barplot <- function(data = data_all, subgroup = FALSE){ 
  
  if (!subgroup){
    p=ggplot(data, aes(x = variable, y = value, fill = Function)) +
      geom_bar(stat = "identity",
               position = "fill",
               width = 0.7) +
      scale_y_continuous(labels = scales::percent) +
      #facet_grid(~ group, scales = "free_x", switch = "x") +
      theme(strip.background = element_blank()) +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
      xlab("Samples") + ylab("Percentage (%)") +
      theme_classic() + theme(axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )) +
      theme(text = element_text(family = "sans", size = 10)) +
      scale_fill_manual(values = colors)
    #scale_fill_brewer(palette="Set3")
  }else{
    p=ggplot(data, aes(x = variable, y = value, fill = Function)) +
      geom_bar(stat = "identity",
               position = "fill",
               width = 0.7) +
      scale_y_continuous(labels = scales::percent) +
      facet_grid(~ Group, scales = "free_x", switch = "x") +
      theme(strip.background = element_blank()) +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
      xlab("Samples") + ylab("Percentage (%)") +
      theme_classic() + theme(axis.text.x = element_text(
        angle = 45,
        vjust = 1,
        hjust = 1
      )) +
      theme(text = element_text(family = "sans", size = 10)) +
      scale_fill_manual(values = colors)
  }
  return(p)
}

plot_stack_barplot()# the result is wired

require(ggplot2)
require(tidyverse)


metadata <- read.csv("../../UKBB_Microbiomes//code/UKBB_NED_filtered.csv") %>% dplyr::rename(SampleID =`Bee_ID` ,
                                                                                             Species = `finalised_species`,
                                                                                             CollectionSite = `Site`)  %>% mutate(interaction = paste0(Species, "_",CollectionSite)) %>%
  filter(parasite_counts_available == "YES" | RADseq_available == "YES", SampleID != "5")

size = metadata %>% group_by(CollectionSite) %>% count() 



# transform nos and est
df = metadata %>% dplyr::select("SampleID","nos_list","eos_list","CollectionSite") %>% drop_na()

library(tidyverse)
library(sf)
#> Linking to GEOS 3.6.1, GDAL 2.2.3, proj.4 4.9.3
df_lat <- df %>%
  st_as_sf(coords = c("eos_list", "nos_list"), crs = 27700) %>%
  st_transform(4326) %>%
  dplyr::mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2]) %>%
  as_tibble() 
#> # A tibble: 6 x 2
#> 
#> 

library("rnaturalearth")
library("rnaturalearthdata")

UK <- ne_countries(scale = "medium", returnclass = "sf",country = "United Kingdom")

p = ggplot(data = UK) +
  geom_sf() +
  geom_point(data =df_lat, mapping = aes(x = lon, y = lat, color = CollectionSite),size = 6) + 
  coord_sf(ylim = c(50, 54),xlim = c(-6,2)) + 
  theme_minimal() + theme(text = element_text(size = 20)) + theme(legend.position = c(0.15,0.8))+
  labs(x = "",y = "")


ggsave(p, filename = "../results/7.Final_graph/GIS_graph.pdf",dpi = "print")

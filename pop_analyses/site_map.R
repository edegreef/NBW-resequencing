# making map with site locations

library(raster)
library(rworldmap)
library(rworldxtra)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsn)
library(rgeos)
#devtools::install_github("MikkoVihtakari/ggOceanMapsData")
#devtools::install_github("MikkoVihtakari/ggOceanMaps")
library(ggOceanMaps)
library(ggOceanMapsData)
library(marmap)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps/pca")

# set map boundary (xmin, xmax, ymin, ymax)
boundary = extent(-75, 16, 37, 73)
boundary

# get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

# download ocean depth map
ocean_map <- getNOAA.bathy(lon1 = -75, lon2 = 16, lat1 = 37, lat2 = 73, resolution = 4)

# load sample info for site locations
sample_info <- read.csv("C:/Users/Evelien de Greef/Dropbox/NBW-me/sample_info_updated_25Jan2021.csv", header=T)

# filter out individuals need to remove
sample_info_36 <- subset(sample_info, remove_indiv!="Y")

# set colors for each pop (assigning pop color to each sample, in order of sample_info file)
site_manual_fill <- c("#287D8EFF","#287D8EFF","#287D8EFF",
                      "#95D840FF",
                      "#481567FF","#481567FF","#481567FF","#481567FF",
                      "#95D840FF",
                      "#FDE725FF","#FDE725FF","#FDE725FF","#FDE725FF",
                      "#95D840FF","#95D840FF","#95D840FF","#95D840FF",
                      "#FDE725FF", "#FDE725FF", "#FDE725FF", "#FDE725FF", "#FDE725FF",
                      "#481567FF","#481567FF","#481567FF","#481567FF","#481567FF","#481567FF","#481567FF",
                      "#C51B7D","#C51B7D","#C51B7D","#C51B7D","#C51B7D","#C51B7D","#C51B7D")

# Plot with ggplot2
nbw_map_ocean <- autoplot(ocean_map, geom=c("raster"), colour="white", size=1) + 
  scale_fill_gradient2(low="#C6DBEF", high="#DEEBF7")+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), fill="grey80",
               colour="grey72", size=0.5)+
  geom_point(data=sample_info_36, aes(x=longitude, y=latitude), pch=21, size=4, colour="black", bg=site_manual_fill)+
  theme(legend.position="none")+
  ylab("latitude")+
  xlab("longitude")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  #ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
  ggsn::scalebar(data = map.outline, dist = 500, dist_unit = "km", height = 0.01,
                 transform = TRUE, model = "WGS84", 
                 location = "bottomleft", anchor = c(x = -70, y = 38),
                 st.bottom = FALSE, st.size = 2.5, st.dist = 0.015)+
  annotate("text", x=-28, y=45, label= "North Atlantic Ocean",fontface="italic", size=3, color="#6191c2")+
  annotate("text", x=-50, y=57, label="Labrador Sea", fontface="italic", size=3, color="#6191c2")+
  annotate("text", x=-2, y=67, label="Norweigan Sea", fontface="italic", size=3, color="#6191c2")
nbw_map_ocean
ggsave("site_map_n36_legend_color_latorder.png", width=9, height=6, dpi=300)

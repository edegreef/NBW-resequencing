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

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps/map")

# set map boundary (xmin, xmax, ymin, ymax)
boundary = extent(-75, 16, 37, 73)
boundary

# get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

# download ocean depth map
ocean_map <- getNOAA.bathy(lon1 = -75, lon2 = 16, lat1 = 37, lat2 = 73, resolution = 4)

# load site info
site_info <-  read.csv("site_locations.csv", header=T)

# set site point colors
site_manual_fill <- c("#4575B4", "#ABD9E9", "#FEE090", "#F46D43", "#A50026")

# for using times new roman font need to do this
library(extrafont)
font_import()
loadfonts(device="win")
fonts()

# plot with ggplot2
nbw_map_ocean <- autoplot(ocean_map, geom=c("raster")) + 
  scale_fill_stepsn(n.breaks=20,colors=c("#125ca1", "#70AED3","#BEDAEC","white", "gray"))+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), 
               fill="gray75",colour="gray40", size=0.5)+
  geom_point(data=site_info, aes(x=longitude, y=latitude), pch=21, size=4, stroke=1, colour="black", bg=site_manual_fill)+
  theme(legend.position="none")+
  ylab("Latitude (°N)")+
  xlab("Longitude (°W)")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  #ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
  ggsn::scalebar(data = map.outline, dist = 500, dist_unit = "km", height = 0.01,
                 transform = TRUE, model = "WGS84", family="Times New Roman",
                 location = "bottomleft", anchor = c(x = -70, y = 38),
                 st.bottom = FALSE, st.size = 2.5, st.dist = 0.015)+
  annotate("text", x=-28, y=45, label= "North Atlantic Ocean",fontface="italic", size=4, color="black",family="Times New Roman")+
  annotate("text", x=-46, y=57, label="Labrador Sea", fontface="italic", size=4, color="black",family="Times New Roman")+
  annotate("text", x=-2, y=67, label="Norweigan Sea", fontface="italic", size=4, color="black", family="Times New Roman")+
  theme(text=element_text(family="Times New Roman", size=15))
nbw_map_ocean

# can stop here, but want to calculate distances between sites outside continental shelf
# calculate distance with -500 depth minimum, making ppath impossible in waters shallower than -500 meters depth

# prepare site coords in necessary format
site_coords <- site_info[,c("longitude", "latitude")]
colnames(site_coords) <- c("x", "y")

trans <- trans.mat(ocean_map,min.depth=-500)
out <- lc.dist(trans,site_coords,res="path")

# 'out' has lists for all the points in distance lines
# extract the output for adding to ggplot map
A <- as.data.frame(out[1])
B <- as.data.frame(out[2])
C <- as.data.frame(out[3])
D <- as.data.frame(out[4])
E <- as.data.frame(out[5])
F <- as.data.frame(out[6])
G <- as.data.frame(out[7])
H <- as.data.frame(out[8])
I <- as.data.frame(out[9])
J <- as.data.frame(out[10])


# add the distance lines to map (adding geom_point again at the end to be on top of the lines)

nbw_map_ocean +
  geom_line(data=A, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_line(data=B, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_line(data=C, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_line(data=D, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_line(data=E, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_line(data=F, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_line(data=H, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_line(data=J, aes(x=x, y=y), color="white", lwd=0.7, linetype="solid")+
  geom_point(data=site_info, aes(x=longitude, y=latitude), pch=21, size=4, stroke=1, colour="black", bg=site_manual_fill)

ggsave("site_map_solo_redone_TNRfont.png", width=9, height=6, dpi=2000)

# actual distance to use in IBD analyses
dist <- lc.dist(trans,site_coords,res="dist")
dist

# 1=IC, 2=AR, 3=LB, 4=NF, 5=SS

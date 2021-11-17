# Making map with site locations and also calculating distances between points

library(raster)
library(rworldmap)
library(rworldxtra)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsn)
library(rgeos)
library(marmap)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/fst")

# Set map boundary (xmin, xmax, ymin, ymax)
boundary = extent(-75, 5, 37, 73)
boundary

# Get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

# Download ocean depth map. Lower resolution number is more fine-scale.
ocean_map <- getNOAA.bathy(lon1 = -75, lon2 = 5, lat1 = 37, lat2 = 73, resolution = 4)

# Load site info
site_info <-  read.csv("site_locations_up.csv", header=T)

# Set site point colors
site_manual_fill <- c("#4575B4", "#ABD9E9", "#FEE090", "#F46D43", "#A50026")

# Plot with ggplot2
nbw_map_ocean <- autoplot(ocean_map, geom=c("raster")) + 
  scale_fill_stepsn(n.breaks=20,colors=c("#125ca1", "#70AED3","#BEDAEC","white", "gray"))+
  geom_polygon(data=map.outline, aes(x=long, y=lat, group=group), 
               fill="gray75",colour="gray40", size=0.5)+
  geom_point(data=site_info, aes(x=longitude, y=latitude), pch=21, size=4, stroke=1, colour="black", bg=site_manual_fill)+
  theme(legend.position="none")+
  ylab("Latitude (°N)")+
  xlab("Longitude (°W)")+
  theme(panel.border = element_rect(colour = "ablack", fill=NA, size=1))+
  #ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft")+
  ggsn::scalebar(data = map.outline, dist = 500, dist_unit = "km", height = 0.01,
                 transform = TRUE, model = "WGS84", family="serif",
                 location = "bottomleft", anchor = c(x = -70, y = 38),
                 st.bottom = FALSE, st.size = 2.5, st.dist = 0.015)+
  annotate("text", x=-28, y=45, label= "North Atlantic Ocean",fontface="italic", size=4, color="black",family="serif")+
  annotate("text", x=-46, y=57, label="Labrador Sea", fontface="italic", size=4, color="black",family="serif")+
  #annotate("text", x=-2, y=67, label="Norweigan Sea", fontface="italic", size=4, color="black", family="serif")+
  theme(text=element_text(family="serif", size=15))
nbw_map_ocean

# Can stop here, but want to calculate distances between sites outside continental shelf
# Calculate distance with -500 depth minimum, making path impossible in waters shallower than -500 meters depth

# Prepare site coords in necessary format
site_coords <- site_info[,c("longitude", "latitude")]
colnames(site_coords) <- c("x", "y")

trans <- trans.mat(ocean_map,min.depth=-500)
out <- lc.dist(trans,site_coords,res="path")

# 'out' has lists for all the points in distance lines
# Extract the output for adding to ggplot map
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

# Add the distance lines to map (adding geom_point again at the end to be on top of the lines)

nbw_map_ocean +
  geom_path(data=A, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_path(data=B, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_path(data=C, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_path(data=D, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_path(data=E, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_path(data=F, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_path(data=H, aes(x=x, y=y), color="white", lwd=0.7,linetype="solid")+
  geom_path(data=J, aes(x=x, y=y), color="white", lwd=0.7, linetype="solid")+
  geom_point(data=site_info, aes(x=longitude, y=latitude), pch=21, size=4, stroke=1, colour="black", bg=site_manual_fill)

ggsave("site_map_up_redo.png", width=9, height=6, dpi=2000)

# Actual distance to use in IBD analyses
dist <- lc.dist(trans,site_coords,res="dist")
dist

# Convert to matrix
mat <- as.matrix(dist)

# Save distance matrix to csv
write.csv(mat, "distance_matrix.csv")

# 1=IC, 2=AR, 3=LB, 4=NF, 5=SS

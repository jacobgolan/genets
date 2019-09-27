################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
################################ ################################ ################################ ################################ 
#Pop maps

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)


setwd("~/Dropbox/Genet Size/Jacob Analyses")

cor<-read.csv("New_Genets_worldmap.csv", head=T, stringsAsFactors = F)
cor.EU<-cor[cor$Map=="EU",]
cor.EU$lab<-paste(cor.EU$Population.Code, cor.EU$year)
cor.US<-cor[cor$Map=="US",]
cor.US$lab<-paste(cor.US$Population.Code, cor.US$year)
cor.CA<-cor[cor$Zoom=="CA",]
cor.CA$lab<-paste(cor.CA$Population.Code, cor.CA$year)

cor.EA<-cor.US[1:3,]


cor.EU.AFLP<-cor.EU[cor.EU$Technology=="AFLP",]
cor.EU.SNP<-cor.EU[cor.EU$Technology=="SNP",]

map.EU <- get_map(location = 'Europe', source ='google', maptype = 'satellite', zoom = 4)
map.CAst <- get_map(location = 'California', source ='google', maptype = 'satellite', zoom = 6)
map.EA <- get_map(location = 'New Jersey', source ='google', maptype = 'satellite', zoom = 6)
map.US <- get_map(location = 'USA', source ='google', maptype = 'satellite', zoom = 4)
map.CA <- get_map(location = 'Point Reyes National Seashore', source ='google', maptype = 'satellite', zoom = 12)
map.Dr <- get_map(location = c(lon=	-122.834134, lat=38.05521201), source ='google', maptype = 'satellite', zoom = 17)

usa <- map_data("usa") # we already did this, but we can do it again

ggplot() + 
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(1.3)+theme_classic()+
  geom_point(data = cor.US, aes(x = Longitude, y = Latitude), color = "red", size = 4)

states <- map_data("state")

ggplot(data = states) + 
  geom_polygon(aes(x = long, y = lat, fill = region, group = group), color = "white") + 
  coord_fixed(1.3) +
  guides(fill=FALSE)  # do this to leave off the color legend

socal <- subset(states, region %in% c("california"))
ggplot(data = socal) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "palegreen", color = "black") + 
  coord_fixed(1.3)


# califonria
ca_df <- subset(states, region == "california")
counties <- map_data("county")
ca_county <- subset(counties, region == "california")
ca_base <- ggplot(data = ca_df, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "gray")
ca_base + theme_nothing()+geom_point(data = cor.CA, aes(x = cor.CA$Longitude, y = cor.CA$Latitude), color = "red", size = 4,inherit.aes = FALSE)
ca_base +  theme_nothing()+
  geom_polygon(data = ca_county, fill = NA, color = "white") +
  geom_polygon(color = "black", fill = NA) + 
  geom_point(data=cor.CA, aes(x = Longitude, y = Latitude),color = "red", size = 4,inherit.aes = FALSE)




# east coast
east <- subset(states, region %in% c("new york", "new jersey"))
counties <- map_data("county")
east_county <- subset(counties, region %in% c("new york", "new jersey"))
east_base <- ggplot(data = east, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "gray")
east_base + theme_nothing()+
  geom_polygon(data = east_county, fill = NA, color = "white") +
  geom_polygon(color = "black", fill = NA) +
  geom_point(data=cor.EA, aes(x = Longitude, y = Latitude),color = "red", size = 4,inherit.aes = FALSE)



####### Europe
europe<-map_data("world")

SNPEU<-read.csv("SNP_EU_latlong.csv", head=T, stringsAsFactors = F)
df1<-cor.EU[1:9,1:3]
df2<-SNPEU[,c(1,6,7)]
colnames(df1)<-c("pop", "Latitude", "Longitude")     
colnames(df2)<-c("pop", "Latitude", "Longitude")     

correct.cor.EU<-rbind(df1,df2)

DIYEU <- subset(europe, region %in% c("France", "Serbia", "Spain", "UK", "Portugal", "Andorra", "Belgium", "Netherlands", "Liechtenstein",
                                      "Luxembourg", "Germany", "Poland","Austria", "Czech Republic", "Ireland", "Slovakia", "Hungary",
                                      "Romania", "Italy", "Switzerland", "Slovenia", "Croatia", "Denmark", "Sweden", "Norway", "Bosnia and Herzegovina",
                                      "Montenegro","Monserrat", "Albania", "Greece", "Bulgaria", "Macedonia", "Kosovo"))

DIYEU.cl<-DIYEU[DIYEU$lat <= 72,]

ggplot(data = DIYEU.cl, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "gray")+
  geom_point(data=correct.cor.EU, aes(x = Longitude, y = Latitude),color = "red", size = 4,inherit.aes = FALSE)+ theme_nothing()


####### Portugal
portugal <- subset(europe, region %in% c("Portugal"))
cor.EU.port<-cor.EU[cor.EU$Population.Code ==c("Agraria", "Vilarinho", "Mira"),]

ggplot(data = portugal, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "gray")+
  geom_point(data=cor.EU.port, aes(x = Longitude, y = Latitude),color = "red", size = 4,inherit.aes = FALSE)+ theme_nothing()

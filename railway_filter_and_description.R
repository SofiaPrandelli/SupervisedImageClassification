########### Italian railways manipulation and descriptive statistics

######## setup #########
setwd("C:/Users/Sofia/PhD/Analisi/railway") # set working directory 

# packages installation
if (!require("pacman")) install.packages("pacman")
pacman::p_load( # p_load function to install.packages + library
  terra, # handle raster data
  sf, # vector data operations
  #  dplyr, # data wrangling
  tidyr, # data wrangling
  data.table, # data wrangling
  tmap, # for mapping
  tiff,
  magrittr,
  ggplot2,
  sp,
  bigreadr,
  tidyverse,
  ggspatial,
  RColorBrewer,
  rgeos,
  maptools,
  rnaturalearth,
  maps,
  devtools,
  viridis,
  rgdal,
  geodata,
  corrplot,
  ggspatial
)

######## Download of Railway data #########

# Data used for this analysis were downloaded from [Geofabrik](http://download.geofabrik.de/europe/italy.html), a consulting and software development firm based in Karlsruhe, Germany, specialized in OpenStreetMap services: data are accurated and up-to-date. After importing and aggregating the data in r, we obtained single sfdf objects.
# Data ARE UP TO DATE 9/04/2024.
CENTRO_rail<- st_read("C:/Users/Sofia/PhD/Analisi/dati_territ/infrastructures/centro/gis_osm_railways_free_1.shp")
ISOLE_rail <- st_read("C:/Users/Sofia/PhD/Analisi/dati_territ/infrastructures/isole/gis_osm_railways_free_1.shp")
NE_rail<- st_read("C:/Users/Sofia/PhD/Analisi/dati_territ/infrastructures/nord_est/gis_osm_railways_free_1.shp")
NO_rail <- st_read("C:/Users/Sofia/PhD/Analisi/dati_territ/infrastructures/nord_ovest/gis_osm_railways_free_1.shp")
SUD_rail <- st_read("C:/Users/Sofia/PhD/Analisi/dati_territ/infrastructures/sud/gis_osm_railways_free_1.shp")
# Dataframe union
railway <- rbind(NO_rail, NE_rail, CENTRO_rail, ISOLE_rail, SUD_rail)

######## Total railway line plot ####### 
ggplot() + 
  geom_sf(data=railway, aes(color=fclass)) +
  ggtitle("Italy railways") + 
  coord_sf()

railway %>% 
  mutate(len = st_length(.)) %>%
  st_drop_geometry() %>% 
  count('fclass')
# "fclass"     n
# 1   fclass 91280

railway %>% 
  mutate(len = st_length(.)) %>% 
  st_drop_geometry() %>% #elimina le geometrie: calcolo più veloce
  as.data.frame() %>% 
  ungroup() %>% 
  group_by(fclass) %>% 
  summarise(tot_len = sum(len, na.rm = T))
# fclass              tot_len
# <chr>                   [m]
# 1 funicular            50870.
# 2 light_rail           45562.
# 3 miniature_railway     6247.
# 4 monorail             16198.
# 5 narrow_gauge       1484042.
# 6 rack                  4313.
# 7 rail              35717408.
# 8 subway              559597.
# 9 tram                707251.

######## Filtering railway data #########
# removing funicular, subway, tram, miniature railway, rack, monorail, tunnel and bridge
rail <- subset(railway, !(fclass %in% c("subway", "funicular", "tram", "miniature_railway", "rack", "monorail")))
rail <- subset(rail, !(tunnel %in% "T"), !(bridge %in% "T"))

rail$osm_id <- as.numeric(rail$osm_id)
st_write(rail, "C:/Users/Sofia/PhD/Analisi/dati_territ/rail2024.shp", append=FALSE)

# -   Miniature railway: are found in parks, public/private gardens and can be tourist attractions in their own right
# -   Monorail: distinctive railway where trains run on one single rail, often suspended above streets
# -   Narrow gauge: passenger or freight trains on narrower tracks than the standard gauge for the country or state
# -   Rack: whatever railway has toothed rack rail that allows operation of trains on very steep gradients

rail %>% 
  mutate(len = st_length(.)) %>% 
  st_drop_geometry() %>% #elimina le geometrie: calcolo più veloce
  as.data.frame() %>% 
  ungroup() %>% 
  group_by(fclass) %>% 
  summarise(tot_len = sum(len, na.rm = T))
# fclass         tot_len
# <chr>              [m]
# 1 light_rail      44771.
# 2 narrow_gauge  1395132.
# 3 rail         32940304.

rail %>% 
  mutate(len=st_length(.)) %>% 
  st_drop_geometry() %>% 
  count('fclass')
# "fclass"     n
# 1   fclass 79032

rail %>% 
  mutate(len = st_length(.)) %>% 
  st_drop_geometry() %>% #elimina le geometrie: calcolo più veloce
  as.data.frame() %>% 
  ungroup() %>% 
  group_by(name) %>% 
  summarise(tot_len = sum(len, na.rm = T))

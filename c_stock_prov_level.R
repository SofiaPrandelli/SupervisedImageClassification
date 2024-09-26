###################### Carbon stock for land use classes - provinces level

######## setup ##########
setwd("/media/r_projects/sofia.prandelli2/wcmc_rail")

#setwd("C:/Users/Sofia/PhD/Analisi/c_stock/provinces_level/") 
#setwd("~/Desktop/analisi/c_stock/")

if (!require("pacman")) install.packages("pacman")
pacman::p_load( # p_load function to install.packages + library
  terra, # handle raster data
  sf, # vector data operations
  #  dplyr, # data wrangling
  tidyr, # data wrangling
  data.table, # data wrangling
  tmap, # for mapping
  stars,
  tiff,
  magrittr,
  ggplot2,
  sp,
  bigreadr,
  tidyverse,
  ggspatial,
  skimr,
  RColorBrewer,
  exactextractr,
  rgeos,
  maptools,
  rnaturalearth,
  maps,
  devtools,
  viridis,
  mapview,
  rgdal,
  geodata,
  Hmisc,
  corrplot,
  ggspatial,
  rio
)

# for ggploot old version: (ci sono problemi nel plottare, potrebbe essere la versione aggiornata)
packageurl <- "http://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

######## WCMC raster #########
# Filtering the dataset in gee: 
wcmc_rast <- rast("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_italy.tif")
# wcmc_df <- as.data.frame(wcmc_rast, xy = TRUE) #convert raster to spat points df
# coordinates(wcmc_df) <- c("x", "y") #convert spatial points df to spat object
# wcmc_sdf <- st_as_sf(wcmc_df) #convert to sf object
wcmc_rast

######## CORINE land cover level III raster #########
corine_rast <- rast("/Volumes/Hard_disk_Sofia/PhD/dati_territ/corinelc_IV/corineIII_italy.tif")

# # crop extent of Corine
# extent_wcmc <- ext(wcmc_100m)
# corine_rast_crop <- crop(corine_rast, extent_wcmc)
# #corine_rast_crop_poly <- as.polygons(corine_rast_crop)

# Disaggregate wcmc raster from 300 m to 100 m res (perché altrimenti ci sarebbero problemi nel determinare le classi di clc cambiando la risoluzione da 100 a 300 m )
wcmc_100m <- disagg(wcmc_rast, fact = 3, method ="bilinear")


######## COPERNICUS LC #########
# Filtrato su gee con extent Italia, tutte le bande di land cover, 100 m res
# https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_Landcover_100m_Proba-V-C3_Global
copern <- rast("/Volumes/Hard_disk_Sofia/PhD/dati_territ/copernicus_lc_ita.tif")

if (crs(wcmc_100m) == crs(corine_rast) && crs(corine_rast) == crs(copern)) {
  print("same CRS")} else {print("different CRS")} #same crs

######## Rasters Stack #########
# Resample rasters: transfers values between non matching Raster objects (in terms of origin and resolution)
wcmc_100_rast_resample <- resample(wcmc_100m, corine_rast, method="bilinear") # weighted average of 4 pixels in the original image nearest to the new pixel location
copern_resample <- resample(copern, corine_rast, method="bilinear")
# Rasters stack (Italy level)
corine_cop_wcmc_stack <- c(wcmc_100_rast_resample, corine_rast, copern_resample) # combines SpatRaster objects with same extent and res
plot(corine_cop_wcmc_stack)

######## Save raster stack ######

# writeRaster(corine_cop_wcmc_stack, "/Users/sofiaprandelli/PhD/c_stock/corine_cop_wcmc_stack.tif")
# analisi su aree_convertibili.R
writeRaster(corine_cop_wcmc_stack, "/media/r_projects/sofia.prandelli2/convertible_rail_areas/data/corine_cop_wcmc_stack_eco_value_classes.tif", overwrite=T)
writeRaster(corine_cop_wcmc_stack_subset, "/media/r_projects/sofia.prandelli2/convertible_rail_areas/data/corine_cop_wcmc_stack_subset_eco_value_classes.tif", overwrite=T)

######## Import raster stack corine_cop_wcmc_stack ########
corine_cop_wcmc_stack <- rast("/media/r_projects/sofia.prandelli2/convertible_rail_areas/data/corine_cop_wcmc_stack_eco_value_classes.tif")
corine_cop_wcmc_stack_subset <- rast("/media/r_projects/sofia.prandelli2/convertible_rail_areas/data/corine_cop_wcmc_stack_subset_eco_value_classes.tif")


######## Italian regions and borders #########
# https://www.paulamoraga.com/tutorial-open-spatial-data/#:~:text=For%20example%2C%20the%20worldclim_country(),codes%20of%20the%20world%20countries
#borders <- ne_countries(type = "countries", country = "Italy", scale = "medium", returnclass = "sf")
regions <- rnaturalearth::ne_states("Italy", returnclass = "sf")

ggplot() +
  geom_sf(data = regions, aes(color=region)) +
  theme_minimal()+
  geom_sf_text(data=regions, aes(label=postal), color="black", size=1, check_overlap = TRUE) +
  coord_sf()
ggsave("province_ita.png", width=7, height=14)

######## WORLDCLIM VARIABLES #########
# Tramite Geodata package, res 1 km (0.08333333 minutes of a degree)
worldclim_bio_ita <- worldclim_country("Italy", res="0.5", var="bio", path=tempdir(), version="2.1")

######## Info worldclim and stack #######
#res(worldclim_prec_ita) # 0.08 degrees = 1 km 
#ext(worldclim_prec_ita) # SpatExtent : 6.5, 19, 35, 47.5 (xmin, xmax, ymin, ymax)
#ext(corine_cop_wcmc_stack) #SpatExtent : 6.62687185094971, 18.5223628432604, 36.6485686462241, 47.0941787699659 (xmin, xmax, ymin, ymax)
#res(corine_cop_wcmc_stack) # 0.0008
if (crs(worldclim_bio_ita) == crs(corine_cop_wcmc_stack)) {
  print("same CRS")} else {print("different CRS")} #same crs

# rast dimensions
nrows <- nrow(worldclim_bio_ita)
ncols <- ncol(worldclim_bio_ita)

######## RICARICA OGGETTI SALVATI ########
corine_cop_wcmc_stack <- rast("C:/Users/Sofia/PhD/Analisi/c_stock/corine_cop_wcmc_stack.tif")
corine_cop_wcmc_stack <- rast("/Volumes/USB/Sofia/PhD/Analisi/c_stock/corine_cop_wcmc_stack.tif")

worldclim_bio_ita <- worldclim_country("Italy", res="0.5", var="bio", path=tempdir(), version="2.1")
regions <- rnaturalearth::ne_states("Italy", returnclass = "sf")
resample_worldclim_bio_ita <- rast("C:/Users/Sofia/PhD/Analisi/c_stock/resample_worldclim_bio_ita.tif")

cropped_raster_list <- readRDS("C:/Users/Sofia/PhD/Analisi/c_stock/provinces_level/cropped_raster_list.tif")

######## Test extents provinces level ########
# Creating extents for every province
# Create a list to store SpatExtent objects for each province
province_extent_list <- list()
# Loop through each province
for (i in 1:nrow(regions)) {
  province <- regions[i, ] # selecting rows of regions df
  province_geometry <- st_geometry(province) #extract geometry
  province_bbox <- st_bbox(province_geometry) #return bounding of sf
  province_extent <- ext(province_bbox) #create SpatExtent object for the province
  province_extent_list[[i]] <- province_extent
}

# Assign region names to every SpatRaster object
for (i in 1:nrow(regions)) {
  names(province_extent_list)[i] <- regions$name[i]
}

######## Stack and crop raster stack and Worldclim ######

# Resample raster worldclim on raster stack with cubic spline interpolation method
# resample_worldclim_bio_ita <- resample(worldclim_bio_ita, corine_cop_wcmc_stack, method="cubicspline")
# writeRaster(resample_worldclim_bio_ita, "C:/Users/Sofia/PhD/Analisi/c_stock/resample_worldclim_bio_ita.tif")
# resample_worldclim_bio_ita <- rast("C:/Users/Sofia/PhD/Analisi/c_stock/resample_worldclim_bio_ita.tif")

# Stack of the 2 rasters
rasters_stack <- c(resample_worldclim_bio_ita, corine_cop_wcmc_stack)
plot(rasters_stack$`urban-coverfraction`)

######## Create a list of cropped rasters_stack #######
cropped_raster_list <- list()
# Loop through each SpatExtent object
for (i in 1:length(province_extent_list)) {
  cropped_raster <- crop(corine_cop_wcmc_stack, province_extent_list[[i]])
  cropped_raster_list[[i]] <- cropped_raster
}

# assign name regions to the list
for (i in 1:nrow(regions)) {
  names(cropped_raster_list)[i] <- regions$name[i]
}

######## Save SpatRaster objects cropped_raster_list ########
output_dir <- "/media/r_projects/sofia.prandelli2/wcmc_rail/provinces_raster_list/"
#output_dir <- "/Volumes/USB/Sofia/PhD/Analisi/c_stock/provinces_level/rasters_prov_cor_cop_wcmc/"

for (i in 1:length(cropped_raster_list)) {
  file_name <- paste0(output_dir, "", names(cropped_raster_list)[i], ".tif")
  writeRaster(cropped_raster_list[[i]], filename=file_name, overwrite=T)
}

######## Import SpatRaster cropped_raster_list as a List ####
output_dir <- "/media/r_projects/sofia.prandelli2/wcmc_rail/provinces_raster_list/"
#output_dir <- "/Volumes/USB/Sofia/PhD/Analisi/c_stock/provinces_level/rasters_prov_cor_cop_wcmc/"
raster_files <- list.files(output_dir, pattern = "\\.tif$", full.names=T)
prov_raster_list <- list()
for (file in raster_files) {
  spat_raster <- rast(file)
  base_name <- tools::file_path_sans_ext(basename(file))
  prov_raster_list[[base_name]] <- spat_raster
}


######## Create a list of df with Mean C values for Corine LC for provinces ########

# Remove NA values from the list cropped_raster_list
# cleaned_cropped_raster_list <- list()
# 
# for (i in 1:length(cropped_raster_list)) {
#   current_raster <- cropped_raster_list[[i]]
#   carbon_layer <- current_raster[["carbon_tonnes_per_ha"]]
#   carbon_layer_cleaned <- carbon_layer[!is.na(carbon_layer)]
#   cleaned_raster <- current_raster
#   cleaned_raster[["carbon_tonnes_per_ha"]] <- carbon_layer_cleaned
#   cleaned_cropped_raster_list[[i]] <- cleaned_raster
# }
# 
# for (i in 1:nrow(regions)) {
#   names(cleaned_cropped_raster_list)[i] <- regions$name[i]
# }

# Calculate random values of C (of values btw Q1 and Q3 quartiles) per Corine land cover
mean_c_values_corine <- list()
for (i in 1:length(cropped_raster_list)) {
  carbon_raster <- cropped_raster_list[[i]][["carbon_tonnes_per_ha"]]
  landcover_raster <- cropped_raster_list[[i]][["landcover"]]
  q_values <- quantile(values(carbon_raster), probs = c(0.25, 0.75), na.rm=T) #calculate first and third quartiles of carbon_tonnes_per_ha
  values <- zonal(carbon_raster, landcover_raster, fun = function(x) { 
    mean(x[x >= q_values[1] & x <= q_values[2]], na.rm=T) #calculate mean of C values between second and third quantile
    #sample(x[x >= q_values[1] & x <= q_values[2]], size=1)
  })
  mean_c_values_corine[[i]] <- values
}
for (i in 1:nrow(regions)) {
  names(mean_c_values_corine)[i] <- regions$name[i]
}

# Aggiunta di errore standard (meglio rispetto alla deviazione standard)
mean_c_values_corine <- list()

sem <- function(x) {
  sd(x, na.rm=T)/sqrt(length(x[!is.na(x)]))
}

for (i in seq_along(cropped_raster_list)) {
  carbon_raster <- cropped_raster_list[[i]][["carbon_tonnes_per_ha"]]
  landcover_raster <- cropped_raster_list[[i]][["landcover"]]
  q_values <- quantile(values(carbon_raster), probs = c(0.25, 0.75), na.rm = TRUE)
  results <- zonal(carbon_raster, landcover_raster, fun = function(x) {
    x_subset <- x[x >= q_values[1] & x <= q_values[2]]
    mean_value <- mean(x_subset, na.rm=T)
    sem_value <- sem(x_subset)
    return(c(mean_value=mean_value, sem_value=sem_value))
  })
  mean_c_values_corine[[i]] <- results
}
# assign provinces names
for (i in 1:nrow(regions)) {
  names(mean_c_values_corine)[i] <- regions$name[i]
}

###### Add description of Corine 44 classes
lc_descriptions <- data.frame(
  landcover = c('111','112',"121","122","123","124","131","132","133","141","142","211","212","213","221","222","223","231","241","242","243","244","311","312","313","321","322","323","324","331","332","333","334","335","411","412","421","422","423","511","512", "521", "522", "523"
  ),
  description = c('Continuous urban fabric','Discontinuous urban fabric','Industrial or commercial units','Road and rail networks and associated land','Port areas','Airports','Mineral extraction sites','Dump sites','Mine, dump, and construction sites','Green urban areas','Sport and leisure facilities','Non-irrigated arable land',
                  'Permanently irrigated land','Rice fields','Vineyards','Fruit trees and berry plantations','Olive groves','Pastures','Heterogeneous agricultural areas: Annual crops associated with permanent crops',
                  'Complex cultivation patterns','Land principally occupied by agriculture, with significant areas of natural vegetation','Agro-forestry areas','Forests, Broad-leaved forest','Forests, Coniferous forest','Forests, Mixed forest','Scrub and/or herbaceous vegetation: Natural grasslands',
                  'Scrub and/or herbaceous vegetation: Moors and heathland','Scrub and/or herbaceous vegetation: Sclerophyllous vegetation','Scrub and/or herbaceous vegetation: Transitional woodland-shrub','Open spaces with little or no vegetation > Beaches, dunes, sands','Open spaces with little or no vegetation > Bare rocks',
                  'Open spaces with little or no vegetation > Sparsely vegetated areas','Open spaces with little or no vegetation > Burnt areas','Open spaces with little or no vegetation > Glaciers and perpetual snow','Inland wetlands: Inland marshes','Inland wetlands: Peat bogs','Maritime wetlands: Salt marshes',
                  'Maritime wetlands: Salines','Maritime wetlands: Intertidal flats','Water courses','Inland waters','Coastal lagoons','Estuaries','Sea and ocean'
  ),
  category = c('Urban fabric','Urban fabric',
               'Industrial, commercial and transport units', 'Industrial, commercial and transport units', 'Industrial, commercial and transport units', 'Industrial, commercial and transport units',
               'Mine, dump and construction sites','Mine, dump and construction sites','Mine, dump and construction sites',
               'Artificial, non-agricultural vegetated areas','Artificial, non-agricultural vegetated areas',
               'Arable land','Arable land','Arable land',
               'Permanent crops','Permanent crops','Permanent crops',
               'Pastures',
               'Heterogeneous agricultural areas','Heterogeneous agricultural areas','Heterogeneous agricultural areas','Heterogeneous agricultural areas',
               'Forests','Forests','Forests',
               'Scrub and/or herbaceous vegetation associations','Scrub and/or herbaceous vegetation associations','Scrub and/or herbaceous vegetation associations','Scrub and/or herbaceous vegetation associations',
               'Open spaces with little or no vegetation','Open spaces with little or no vegetation','Open spaces with little or no vegetation','Open spaces with little or no vegetation','Open spaces with little or no vegetation',
               'Inland wetlands','Inland wetlands',
               'Maritime wetlands','Maritime wetlands','Maritime wetlands',
               'Inland waters','Inland waters',
               'Marine waters','Marine waters','Marine waters'
  )
)
export(lc_descriptions, "/media/r_projects/sofia.prandelli2/wcmc_rail/lc_corine_descriptions.xlsx")

# fn to add descriptions to df
add_descriptions <- function(df) {
  df$description <- lc_descriptions$description[match(df$landcover, lc_descriptions$landcover)]
  df$category <- lc_descriptions$category[match(df$landcover, lc_descriptions$landcover)]
  return(df)
}

# loop through each df in c_values_corine list and add descriptions
for (i in 1:length(mean_c_values_corine)) {
  mean_c_values_corine[[i]] <- add_descriptions(mean_c_values_corine[[i]])
}

view(mean_c_values_corine$Bologna)

#### sistemazione dataframe per esportazione
extract_val <- function(df) {
  mean_value <- df$carbon_tonnes_per_ha[, "mean_value"]
  sem_value <- df$carbon_tonnes_per_ha[, "sem_value"]
  df_extracted <- data.frame(
    landcover=df$landcover,
    description=df$description,
    category=df$category,
    mean_value=mean_value,
    sem_value=sem_value
  )
  return(df_extracted)
}

mean_c_val_corine <- lapply(mean_c_values_corine, extract_val)

##### REMOVE ALL THE CLASSES 1.Artifical Surfaces (1.1, 1.2, 1.3, 1.4) EXCEPT FOR 141 (Green urban areas) #####
exclude_values <- c(111, 112, 121, 122, 123, 124, 131, 132, 133, 142)

mean_carbon_soil_clc <- lapply(mean_c_val_corine, function(df) {
  df <- df[!df$landcover %in% exclude_values, ]
  return(df)
})

export(mean_carbon_soil_clc, "/media/r_projects/sofia.prandelli2/wcmc_rail/mean_carbon_soil_clc.xlsx") # library(rio)
export(mean_c_val_corine, "/media/r_projects/sofia.prandelli2/wcmc_rail/mean_carbon_soil_clc_withurban.xlsx") # library(rio)

####### Create list of df with only the mean of principal classes (associated to the ones created in google earth engine supervised classification) from the object mean_c_val_corine ########
mean_c_val_principal_classes <- list()

for (df_name in names(mean_c_val_corine)) {
  df <- mean_c_val_corine[[df_name]]
  trees_mean <- mean(df$mean_value[df$landcover %in% c(311, 312, 313)], na.rm = TRUE)
  shrubs_mean <- mean(df$mean_value[df$landcover == 324], na.rm = TRUE)
  grass_mean <- mean(df$mean_value[df$landcover == 321], na.rm = TRUE)
  bare_mean <- mean(df$mean_value[df$landcover %in% c(332, 333)], na.rm = TRUE)

  new_df <- data.frame(
    class = c("trees", "shrubs", "grass", "bare"),
    mean_value = c(trees_mean, shrubs_mean, grass_mean, bare_mean)
  )

  mean_c_val_principal_classes[[df_name]] <- new_df
}

export(mean_c_val_principal_classes, "/media/r_projects/sofia.prandelli2/wcmc_rail/mean_c_val_principal_classes.xlsx")
# sul pc è salvato in C:\Users\Sofia\PhD\Analisi\c_stock\provinces_level\dataset excel come mean_c_val_principal_classes_macrofunction
# aggiungendo la funzione macro per il calcolo immediato del carbonio stoccato in tot m2

##### PLOT mean_carbon_soil_clc landcover distribution across provinces ######
# library(ggplot2)
# library(dplyr)
setwd("/media/r_projects/sofia.prandelli2/wcmc_rail/plots_c_stock_per_clc")

aggregated_data <- list()

# Loop through each df in the list
for (province_name in names(mean_carbon_soil_clc)) {
  df <- mean_carbon_soil_clc[[province_name]]
  df$province <- province_name
  aggregated_data[[province_name]] <- df
}

# Combine all df of the list into one
combined_data <- bind_rows(aggregated_data)
# get unique descriptions
unique_descriptions <- unique(combined_data$description)

for (desc in unique_descriptions) {
  # filter the combined data for the current description
  plot_data <- combined_data %>% filter(description == desc)
  if (nrow(plot_data) > 0) {
    p <- ggplot(plot_data, aes(x = province, y = mean_value)) +
      geom_point(size =0.5, color = "black") +         
      geom_line(aes(group=1), size=0.5, color = "black") +                    
      labs(title = desc, x = "Italian province", y = "Carbon stock mean value(T/ha)") +
      theme_minimal() +                              
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=4)) +
      scale_y_continuous(seq(0, max(plot_data$mean_value, na.rm = TRUE), by = 10))
     # print(p)
      ggsave(filename = paste0("c_stock_", gsub(" ", "_", desc), ".png"), plot=p, width=10, height=6)
  }
}

########## single plot for clc 313 ####
plot_data <- combined_data %>% filter(landcover=="313")

ggplot(plot_data, aes(x = province, y = mean_value)) +
  geom_point(size = 0.5, color = "black") +         
  geom_line(aes(group=2), size=0.5, color = "darkgrey") +                    
  labs(title = "Agro-forestry areas", x = "Italian province", y = "Carbon stock mean value(T/ha)") +
  #theme_minimal() + 
  theme_bw() +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
  theme(axis.text.x = element_text(angle=45, hjust=1, size=4), plot.margin = margin(t=10, r=10, b=30, l=10)) +
  scale_y_continuous(seq(0, max(plot_data$mean_value, na.rm = TRUE), by=10))


######## PLOT c_stock in function of lc classes for provinces ##########
setwd("C:/Users/Sofia/PhD/Analisi/c_stock/provinces_level/plots/")

category_colors <- c('Urban fabric' = 'grey', 'Industrial, commercial and transport units' ='lightgray',
                     'Mine, dump and construction sites'='darkkhaki', 'Artificial, non-agricultural vegetated areas'='lightgreen',
                     'Arable land'='gold2' , 'Permanent crops'='darkolivegreen3', 'Pastures'='palegreen3', 'Heterogeneous agricultural areas'='greenyellow','Forests'='darkgreen', 
                     'Scrub and/or herbaceous vegetation associations'='darkolivegreen2', 'Open spaces with little or no vegetation'='yellow',
                     'Inland wetlands'='darkslategray1', 'Maritime wetlands'='darkslategray2', 'Inland waters'='darkslategray3', 'Marine waters'='royalblue')

# Convert category_colors to a data frame
legend_data <- data.frame(Category = names(category_colors), Color = category_colors, stringsAsFactors = FALSE)

######## Bologna
bologna <- mean_c_val_corine[["Bologna"]]
length(unique(bologna$category))
view(bologna)

ggplot(bologna, aes(x = as.factor(landcover), y = mean_value, fill=category)) +
  geom_col() +
  labs(y ="C t/ha", x="Corine LC", title = paste("Distribution of carbon t/ha by Corine Land Cover classes - Bologna province")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #    scale_fill_viridis_d()+ # non ha abbastanza colori per le mie classi 
  scale_fill_manual(values=category_colors, name="Macro-classes of Land Cover")+
  #  scale_fill_gradient(low="yellow", high="darkgreen") +
  #  scale_fill_gradientn(colors=terrain.colors(num_categ), name = "Values of C t/ha")
  scale_y_continuous(limits=c(0,60))
ggsave("c_stock_corine_bologna.png", width=12, height=7)

sapply(lapply(bologna, unique), length)


############# LOOP FOR ALL PROVINCES
# print max values
library(purrr)
highest_values <- map_dbl(mean_c_val_corine, ~ max(.$carbon_tonnes_per_ha))
highest_values <- as.data.frame(highest_values)
max(highest_values, na.rm=T) # 84.13283 C t/ha

plots_list <- list()

for (i in seq_along(mean_C_values_corine)) {
  df <- mean_C_values_corine[[i]]
  
  plot <- ggplot(df, aes(x=as.factor(landcover), y=carbon_tonnes_per_ha, fill=category)) +
    geom_col() +
    labs(x ="Corine LC", y="C t/ha", title = paste("Distribution of carbon_tonnes_per_ha by landcover -", names(mean_C_values_corine)[i])) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values=category_colors, name="Macro-classes of Land Cover")+
    scale_y_continuous(limits=c(0,90))
  plots_list[[i]] <- plot
  print(plot)
  #  ggsave(paste0("c_stock_classes_", names(mean_C_values_corine)[i], ".png"), width=12, height=7)
}


######## Plot df without Urban land cover categories ########
# Removing Urban fabrics and industrial, commercial and transport units
cat_to_filter <- c('Continuous urban fabric', 'Discontinuous urban fabric', 'Industrial or commercial units','Port areas','Airports') # Road and rail networks and associated land, 'Mineral extraction sites','Dump sites','Mine, dump and construction sites','Sport and leisure facilities', 'Green urban areas'

filt_mean_C_values_corine <- list()

for (i in seq_along(mean_C_values_corine)) {
  filt_mean_C_values_corine[[i]] <- mean_C_values_corine[[i]][!mean_C_values_corine[[i]]$description %in% cat_to_filter, ]
}

for (i in 1:nrow(regions)) {
  names(filt_mean_C_values_corine)[i] <- regions$name[i]
}

# Plot
setwd("C:/Users/Sofia/PhD/Analisi/c_stock/provinces_level/plots/without_urban")
plots_list <- list()

for (i in seq_along(filt_mean_C_values_corine)) {
  df <- filt_mean_C_values_corine[[i]]
  plot <- ggplot(df, aes(x=as.factor(landcover), y=carbon_tonnes_per_ha, fill=category)) +
    geom_col() +
    labs(x ="C t/ha", y="Corine LC", title = paste("Distribution of carbon_tonnes_per_ha by landcover -", names(mean_C_values_corine)[i])) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values=category_colors, name="Macro-classes of Land Cover")+
    scale_y_continuous(limits=c(0,85))
  ggsave(paste0("c_stock_classes_", names(filt_mean_C_values_corine)[i], ".png"), width=12, height=7)
}

export(filt_mean_C_values_corine, "c_stock_corine_ita_without_urban.xlsx")



######## COVER FRACTION COPERNICUS ANALYSIS #############

####### Convert SpatRaster List to DF List
cropped_rast_list_df <- lapply(cropped_raster_list, function(raster) {
  df <- as.data.frame(raster, xy = TRUE)   
  return(df)
})

# Clean list of df
func_clean_df <- function(df){
  colnames(df) <- gsub("-", "_", colnames(df))
  colnames(df) <- gsub("wc2.1_30s_", "", colnames(df)) # questa riga se ci sono anche le Worldclim bioclimatic variables
  df$landcover <- as.factor(df$landcover)
  df$forest_type <- round(as.numeric(df$forest_type), 0)
  return(df)
}
cropped_rast_list_df <- lapply(cropped_rast_list_df, func_clean_df)
read(cropped_rast_list_df, "cropped_rast_list_df.xlsx")

######## Create a list of nested dataframes #########
# Define the coverfraction types and the corresponding column names
coverfraction_types <- c("bare", "urban", "crops", "grass", "shrub", "tree")
coverfraction_columns <- list(
  c("carbon_tonnes_per_ha", "bare_coverfraction", "forest_type"),
  c("carbon_tonnes_per_ha", "urban_coverfraction", "forest_type"),
  c("carbon_tonnes_per_ha", "crops_coverfraction", "forest_type"),
  c("carbon_tonnes_per_ha", "grass_coverfraction", "forest_type"),
  c("carbon_tonnes_per_ha", "shrub_coverfraction", "forest_type"),
  c("carbon_tonnes_per_ha", "tree_coverfraction", "forest_type")
)

# Function to extract variables and create dataframes
extract_variables <- function(df, columns) {
  dfs <- list()
  for (i in seq_along(columns)) {
    # Extract variables for current coverfraction type
    df_subset <- df[, columns[[i]]]
    # Name the dataframe based on coverfraction type
    df_name <- paste("df_", coverfraction_types[i], sep = "")
    # Store the named dataframe in the list
    dfs[[df_name]] <- df_subset
  }
  return(dfs)
}

# Apply the function to each df in cropped_rast_list_df
final_list <- lapply(cropped_rast_list_df, function(df) {
  extract_variables(df, coverfraction_columns)
})

# assign name regions to the list
for (i in 1:nrow(regions)) {
  names(final_list)[i] <- regions$name[i]
}

view(final_list$Bologna$df_tree)


# NON FUNZIONA
# save(final_list, file="C:/Users/Sofia/PhD/Analisi/c_stock/provinces_level/c_stock_coverfractions.xlsx")
# final_list <- load("C:/Users/Sofia/PhD/Analisi/c_stock/provinces_level/c_stock_coverfractions.xlsx")

# Save each list as a separate RDS file
# DA SISTEMARE
setwd("C:/Users/Sofia/PhD/Analisi/c_stock/provinces_level/cover_fractions_df")
for (i in seq_along(final_list)) {
  saveRDS(final_list[[i]], file = paste0("list_", i, ".xlsx"))
}

######## Filter only values above 80 for the _coverfraction variables #######

# SI BLOCCA !!!
# forse perché non per tutti i df ho valori che superano l'80%
# lavoro con i percentili? in questo caso però non avremmo valori di cs per fraction cover massima, ci interessa?

# Function to filter dataframes considering values above 80 for cover fraction columns

# Define a function to filter values above 80 for the coverfraction columns
filter_coverfraction_above_80 <- function(df) {
  df_filtered <- df %>%
    mutate(across(
      starts_with(c("bare", "urban", "crops", "grass", "shrub", "tree")),  # Select columns
      ~ ifelse(. >= 80, ., NA)  # Replace values < 80 with NA
    ))
  return(df_filtered)
}

# Apply the function to each nested dataframe in final_list
final_list_filtered <- lapply(final_list, function(nested_df_list) {
  lapply(nested_df_list, filter_coverfraction_above_80)
})



######## plot boxplots or scatterplot depending on the variable coverfraction









######## Intersection with railway shapefile #################
# Simplify geometry of railway shapefile
# Buffer of 100 m ?
# Check crs 






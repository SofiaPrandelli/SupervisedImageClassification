# Supervised Image Classification
### Summary
Contains Google Earth Engine (GEE) JavaScript code to classify and quantify areas adjacent to rail infrastructures. The Supervised Image Classification was performed through the Machine Learning Random Forest Model, considering an area equal to the intersection of high resolution areial photos in raster format with Italian railways infrastructures in shapefile format. 

### Preliminary operations
Make sure to have activated the account on Google Earth Engine and to have correctly uploaded and renamed all files needed for the analysis: copy the "Java Script Code" file from the folder to the "Script" section in your GEE account, by creating a new file; while in the "Assets" section, upload the raster images of the area of interest and the railroad shapefile.

### Test area
A mosaic of aerial images of a railway section in Guidonia Montecelio (Rome) was classified; the images were uploaded in the GEE platform and their resolution has been decreased to 1 meter. The region of interest (ROI) was defined by considering the geometry of the mosaic created. The images used are availabe upon request due to the large dimensons of the files.
It is possible to change the test area by uploading in the GEE platform the georeferenced tiff aereal images of a new area of interest and by changing the parameters in line 17 of the provided script.

### Railways infrastructures
Data used for this analysis were downloaded from Geofabrik (http://download.geofabrik.de/europe/italy.html) (Geofabrik, D.S., 2017), a consulting and software development firm based in Karlsruhe, Germany, specialized in OpenStreetMap services: data are accurated and up-to-date (2024). Data has been pre-processed in the R software, by aggregating and filtering them by considering only suitable categories, corresponding to "rail", "light rail" and "narrow gauge" infrastructures. Sections on tunnel and bridges were also removed. Thereafter, the processed shapefile was uploaded in the GEE platform and filtered in the test area created, equivalent to the extent of the images used. 
The R script containing the changes made to the original shapefile and of the descriptive statistics regarding the types of rail infrastructure at the Italian level is available in the document “railway_filter_and_description.R”.

### Buffer zones
The Supervised Classification was performed more than once, by filtering the mosaic with buffers around the railway infrastructure of different width. It is possible to change this parameter at line 50 of the script, that is currently setted at 10 meters.

### Groun Control Points
Once the image is ready, we need to create the points that represent each landcover class in the classification. In this case, the classes created were: rails, trees, grass, bare and shrubs vegetation.
Using the imagery as guidance (or the Google Earth Engine Satellite image), hover over the ‘Geometry Imports’ box next to the geometry drawing tools and click ‘+ new layer.’ Each new layer represents one class within the training data. Let the first new layer represent "rails". Locate points in the new layer in railway infrastructure area. When finished collecting points, click ‘Exit’ and configure the import (top of the script) as follows. Name the layer "rails" and click the icon to configure it. ‘Import as’ FeatureCollection. ‘Add property’ landcover and set its value to 0. (Subsequent classes will be 1 for trees, 2 for grass, etc.). 
To simplify the procedure, it is possible to uncomment from line 84 to line 88 of the script to automatically create the classes. Then, start adding the Ground Control Points per each layer. Note that the number of points should be homogeneous for all land cover categories and the pixel values ​​of a given land cover class should be approximately all represented, encompassing all spectral variances that are in a class.

<!-- Keep in mind that you want to encompass all spectral variances that are in a class. We will look into assessing the quality of polygons in the next section. In general, each class should be normally distributed and not overlapping substantially with any other class. da spiegare meglio -->

### Supervised Classification
Using ee.Classifier.randomForest() function, it is possible to instantiate a classifier and train it on the training data specifying the features to use (training), the landcover categories as the classProperty we want to categorize the imagery into and the reflectance in Bands of the raster image as the inputProperties.
After classifiyng the image (Script Lines from 122 to 133) and display the results (using the Inspector tool, it is possible to visualize the caracteristics of each classified pixel), we assess the Accuracy of the model and calculate the Confusion Matrix. 
In this particular example, we are just looking at the trainAccuracy, which basically describes how well the classifier was able to correctly label resubstituted training data (i.e. data that the classifier had already seen).

### Calculation of area by class
In the section of the script from line 197 to line 273, the sum of the areas in terms of square meters and percentage for each class of our classification is calculated; results can be displayed through bar chart and pie chart.

### Copernicus Corine Land Cover overview
In the last section of the script, an additional overview of the test area is given by considering the Copernicus Corine Land Cover Level III datasets, that contains 44 categorical values of Land Use (Copernicus, E.U., 2018). Total areas for each of the 44 classes are calculated and displayed by charts. 

### Data Paper
A brief theoretical and pratical report to explain the specific study case used as line guides in GEE for the Supervised Image Classification is attached under the name "Report_Guidonia".

# Quantification of Biomass Carbon Density  
### Provincial-level dataset
Following the recognition of the different soil types found in a given area of interest for railroad planning, we hypothesized to quantify in terms of carbon potentially stored in soil the differences that can be found in different land use classes, considering both the different environmental and geographic conditions across the Italian territory and the political administrative divisions. The data made available in this section are contained in the excel files "mean_carbon_soil_clc_macrofunction.xlsm" and "mean_c_val_principalclasses_macrofunction.xlsm". After a literature review of existing datasets on biomass carbon in terrestrial ecosystems, analysis was conducted through R Data manipulation of the open access Database of Above and Below Ground Biomass Carbon Density from UNEP-WCMC (UN Environment Programme World Conservation Monitoring Centre), representing the ground terrestrial carbon storage (C t/ha) for circa 2010. The whole dataset was filtered through Google Earth Engine, downloaded and uploaded in R. The initial raster resolution was refined by Cubic spline interpolation method and interpolated with Copernicus Corine Land Cover Level III Database, obtaining a Dataset with mean values  and Standard Error of the mean of Biomass Carbon Density for each class of Land Use, divided per italian provinces. Results show substantial differences in soil stored Carbon between land use types for some provinces, but differences can be observed also between the average values of the same land use classes across provinces due to differences in altitudinal, latitudinal, and climatic conditions. Note that approximation errors may be present due to the operations that increased the resolution of the raster from 300 meters to 100 meters (by interpolating values between adjacent raster cells to create new values), and associating the new values with the corresponding land cover values, so it is possible that they were associated with carbon quantities of a different class of land cover in the reality. The "mean_carbon_soil_clc_macrofunction.xlsm" file consists in a list of dataset, each one containing the average carbon stock values of each of the land use types found in a given province.
The file "mean_c_val_principalclasses_macrofunction.xlsm" contains the average values ​​of the four land use classes (tree, grass, bare and shrub vegetation) identified by the Supervised Image Classification, so as to be able to have a larger scale overview and quantify the soil carbon stock in a specific area of ​​interest. These values were obtained by calculating the mean value of Tonnes of soil carbon stock per hectares of different classes of Corine Land Cover, in particular:
- Values of the class "tree" were obtained by calculating the mean value between values of the classes "Forests, Broad-leaved forest" and "Forests, Coniferous forest";
- Values of the class "grass" were obtained by considering the mean value of the class "Scrub and/or herbaceous vegetation: Natural grasslands";
- Values of the class "shrub" were obtained by considering the mean value of the class "Scrub and/or herbaceous vegetation: Transitional woodland-shrub";
- Values of the class "bare" were obtained by calculating the mean value between values of the classes"Open spaces with little or no vegetation > Bare rocks" and "Open spaces with little or no vegetation > Sparsely vegetated areas".
The "rail" class was excluded by this calculation, since a carbon stock soil for this kind of land cover was considered equal to 0 C T/Ha.

The R script containing all analyses for the creation of carbon stock soil datasets is available among the attachments as "c_stock_prov_level.R".

These findings are in line with broader patterns observed in Italian landscapes. According to Chiti et al. (2012), Soil Organic Carbon (SOC) stocks vary significantly across different climatic regions. The highest SOC stocks were found in the temperate sub-oceanic regions' lowland depressions and in the alluvial soils of the warm temperate continental climate, while the lowest stocks were observed in the Mediterranean sub-continental region. Agroforestry systems show SOC stocks ranging from 40.1 ± 12.3 Mg C ha−1 in the Mediterranean sub-tropical region to 70.1 ± 23.3 Mg C ha−1 in the temperate mountain climate, reflecting significant differences based on climate type. Generally, mature forests typically have substantially larger soil carbon stocks than agricultural soils, that generally have lower carbon stocks due to cultivation practices, with values often between 30 and 80 Mg C ha^-1 in the top 30 cm (Lal, 2004).

<!-- DA SISTEMARE ULTIMI PARAGRAFI? -->
In temperate forests, SOC stocks can range from 60 to 130 Mg C ha^-1 in the top 30 cm of soil (Lal, 2005). However, these values can be much higher in certain forest types, with some temperate deciduous forests storing up to 200-300 Mg C ha^-1 in the soil (Jobbagy and Jackson, 2000). Grasslands in temperate regions also have substantial soil carbon stocks, typically ranging from 60 to 150 Mg C ha^-1 in the top 30 cm (Conant et al., 2017). Land use change, particularly the conversion of natural ecosystems to agriculture, can lead to significant losses of soil carbon. For instance, the conversion of temperate grasslands to croplands has been estimated to result in a 25-30% reduction in soil carbon stocks in the top 30 cm over several decades (Guo and Gifford, 2002).

Above-ground biomass carbon stocks also differ significantly across land uses. In temperate forests, above-ground biomass carbon can range from 60 to 250 Mg C ha^-1, depending on forest type, age, and management history (Keith et al., 2009). Old-growth temperate forests can accumulate even higher amounts, with some studies reporting above-ground carbon stocks exceeding 500 Mg C ha^-1 in mature Douglas-fir forests of the Pacific Northwest (Harmon et al., 2004). In contrast, temperate grasslands typically have much lower above-ground biomass carbon, usually less than 5 Mg C ha^-1 (Gibson, 2009). Agricultural systems show wide variation, with annual croplands having minimal persistent above-ground carbon, while perennial crops and agroforestry systems can accumulate substantial biomass carbon over time.
Below-ground biomass carbon is often overlooked but can represent a significant portion of total ecosystem carbon, particularly in grasslands and shrublands. In temperate grasslands, below-ground biomass carbon can range from 3 to 20 Mg C ha^-1, often exceeding above-ground biomass (Jackson et al., 1996). In temperate forests, below-ground biomass carbon typically represents about 20-30% of the above-ground biomass carbon (Cairns et al., 1997), translating to roughly 15-75 Mg C ha^-1 for most temperate forest types.

### Contact
For any queries, contact sofia.prandelli2@unibo.it

### References
- Chiti, T., Gardin, L., Perugini, L., Quaratino, R., Vaccari, F. P., Miglietta, F., & Valentini, R. (2012). Soil organic carbon stock assessment for the different cropland land uses in Italy. Biology and Fertility of Soils, 48, 9-17.
- Cairns, M.A., Brown, S., Helmer, E.H., Baumgardner, G.A., (1997). Root biomass allocation in the world's upland forests. Oecologia 111, 1-11.
- Copernicus, E. U. (2018). Copernicus land monitoring service. Corine Land Cover data.
- Cunningham, S. C., Cavagnaro, T. R., Mac Nally, R., Paul, K. I., Baker, P. J., Beringer, J., ... & Thompson, R. M. (2015). Reforestation with native mixed‐species plantings in a temperate continental climate effectively sequesters and stabilizes carbon within decades. Global Change Biology, 21(4), 1552-1566.
- Geofabrik, D. S. (2017). Open Street Map.
- Gibson, D.J., (2009). Grasses and Grassland Ecology. Oxford University Press, Oxford.
- Guo, L.B., Gifford, R.M., (2002). Soil carbon stocks and land use change: a meta analysis. Global Change Biology 8, 345-360.
- Harmon, M.E., Moreno, A., Domingo, J.B., (2009). Effects of partial harvest on the carbon stores in Douglas-fir/western hemlock forests: a simulation study. Ecosystems 12, 777-791.
- Jackson, R.B., Canadell, J., Ehleringer, J.R., Mooney, H.A., Sala, O.E., Schulze, E.D., (1996). A global analysis of root distributions for terrestrial biomes. Oecologia 108, 389-411.
- Jobbagy, E.G., Jackson, R.B., (2000). The vertical distribution of soil organic carbon and its relation to climate and vegetation. Ecological Applications 10, 423-436.
- Keith, H., Mackey, B.G., Lindenmayer, D.B., (2009). Re-evaluation of forest biomass carbon stocks and lessons from the world's most carbon-dense forests. Proceedings of the National Academy of Sciences 106, 11635-11640.
- Lal, R., (2004). Soil carbon sequestration to mitigate climate change. Geoderma, 123(1-2), 1-22.
- Lal, R., (2005). Forest soils and carbon sequestration. Forest Ecology and Management 220, 242-258.







# SupervisedImageClassification
### Summary
Contains Google Earth Engine (GEE) JavaScript code to classify and quantify areas adjacent to rail infrastructures. The Supervised Image Classification was performed through the Machine Learning Random Forest Model, considering an area equal to the intersection of high resolution areial photos in raster format with Italian railways infrastructures in shapefile format. 

### Test area
A mosaic of aerial images of a railway section in Guidonia Montecelio (Rome) was classified; the images were uploaded in the GEE platform and their resolution has been decreased to 1 meter. The region of interest (ROI) was defined by considering the geometry of the mosaic created. The images used are availabe upon request due to the large dimensons of the files. It is possible to change the test area by uploading in the GEE platform the georeferenced tiff aereal images of a new area of interest and by changing the parameters in line 17 of the provided script.

### Railways infrastructures
Data used for this analysis were downloaded from Geofabrik (http://download.geofabrik.de/europe/italy.html), a consulting and software development firm based in Karlsruhe, Germany, specialized in OpenStreetMap services: data are accurated and up-to-date (2024). Data has been pre-processed in the R software, by aggregating and filtering them by considering only suitable categories, corresponding to "rail", "light rail" and "narrow gauge" infrastructures. Sections on tunnel and bridges were also removed. Thereafter, the processed shapefile was uploaded in the GEE platform and filtered in the test area created, equivalent to the extent of the images used.

### Buffer zones
The Supervised Classification was performed more than once, by filtering the mosaic with buffers around the railway infrastructure of different width. It is possible to change this parameter at line 50 of the script, that is actually setted at 10 meters.

### Groun Control Points
Once the image is ready, we need to create the points that represent each landcover class in the classification. In this case, the classes created were: rails, trees, grass, bare and shrubs vegetation.
Using the imagery as guidance (or the Google Earth Engine Satellite image), hover over the ‘Geometry Imports’ box next to the geometry drawing tools and click ‘+ new layer.’ Each new layer represents one class within the training data. Let the first new layer represent "rails". Locate points in the new layer in railway infrastructure area. When finished collecting points, click ‘Exit’ and configure the import (top of the script) as follows. Name the layer "rails" and click the icon to configure it. ‘Import as’ FeatureCollection. ‘Add property’ landcover and set its value to 0. (Subsequent classes will be 1 for trees, 2 for grass, etc.). 
To simplify the procedure, it is possible to uncomment from line 84 to line 88 of the script to automatically create the classes. Then, start adding the Ground Control Points per each layer. Note that the number of points should be homogeneous for all land cover categories and the pixel values ​​of a given land cover class should be approximately all represented, encompassing all spectral variances that are in a class.

<!-- Keep in mind that you want to encompass all spectral variances that are in a class. We will look into assessing the quality of polygons in the next section. In general, each class should be normally distributed and not overlapping substantially with any other class. da spiegare meglio -->

### Supervised classification
Using ee.Classifier.randomForest() function, it is possible to instantiate a classifier and train it on the training data specifying the features to use (training), the landcover categories as the classProperty we want to categorize the imagery into and the reflectance in Bands of the raster image as the inputProperties.
After classifiyng the image (Script Lines from 122 to 133) and display the results (using the Inspector tool, it is possible to visualize the caracteristics of each classified pixel), we assess the Accuracy of the model and calculate the Confusion Matrix. 
In this particular example, we are just looking at the trainAccuracy, which basically describes how well the classifier was able to correctly label resubstituted training data, i.e. data the classifier had already seen.

### Calculation of area by class
In the last section of the script, the sum of the areas for each class of Land Cover is calculated in square meters, and results are displayed through bar chart or pie chart.

### Contact
For any queries, contact sofia.prandelli2@unibo.it

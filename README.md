# SupervisedImageClassification
### Summary
Contains Google Earth Engine (GEE) JavaScript code to classify and quantify areas adjacent to rail infrastructures. The Supervised Image Classification was performed through the Machine Learning Random Forest Model, considering an area equal to the intersection of high resolution areial photos in raster format with Italian railways infrastructures in shapefile format. 
### Test area
A mosaic of aerial images of a railway section in Guidonia Montecelio (Rome) was classified; the images were uploaded in the GEE platform and their resolution has been decreased to 1 meter. The region of interest (ROI) was defined by considering the geometry of the mosaic created. The images used are availabe upon request due to the large dimensons of the files. It is possible to change the test area by uploading in the GEE platform the georeferenced tiff aereal images of a new area of interest and by changing the parameters in line 17 of the provided script.
### Railways infrastructures
Data used for this analysis were downloaded from Geofabrik (http://download.geofabrik.de/europe/italy.html), a consulting and software development firm based in Karlsruhe, Germany, specialized in OpenStreetMap services: data are accurated and up-to-date (2024). Data has been pre-processed in the R software, by aggregating and filtering them by considering only suitable categories, corresponding to "rail", "light rail" and "narrow gauge" infrastructures. Sections on tunnel and bridges were also removed. Thereafter, the processed shapefile was uploaded in the GEE platform and filtered in the test area created equivalent to the extent of the images used.
### Buffer zones
The Supervised Classification was performed more than once, by filtering the mosaic with buffers around the railway infrastructure of different width. It is possible to change this parameter at line 50 of the script, that is actually setted at 10 meters. 
### Groun Control Points
The second step is to collect training data. Using the imagery as guidance, hover over the ‘Geometry Imports’ box next to the geometry drawing tools and click ‘+ new layer.’ Each new layer represents one class within the training data. Let the first new layer represent ‘urban.’ Locate points in the new layer in urban or built up areas (buildings, roads, parking lots, etc.). When finished collecting points, click ‘Exit’ and configure the import (top of the script) as follows. Name the layer ‘urban’ and click the icon to configure it. ‘Import as’ FeatureCollection. ‘Add property’ landcover and set its value to 0. (Subsequent classes will be 1 for water, 2 for forest, etc.) when finished, click ‘OK’ as shown:

### Supervised classification
The Classifier package handles supervised classification by traditional ML algorithms

In this particular example, we are just looking at the trainAccuracy, which basically describes how well the classifier was able to correctly label resubstituted training data, i.e. data the classifier had already seen. 

### Contact
For any queries, contact sofia.prandelli2@unibo.it

# SupervisedImageClassification
### Summary
Contains Google Earth Engine (GEE) JavaScript code to classify and quantify areas adjacent to rail infrastructures. The Supervised Image Classification was performed through the Machine Learning Random Forest Model, using high resolution areial photos in raster format. 
### Test area
As test area, a mosaic of aerial images of a railway section in Guidonia Montecelio (Rome) was classified; the images were uploaded in the GEE platform and their resolution has been decreased to 1 meter. The image used are availabe upon request due to the large dimensons of the files.
It is possible to change the test area by uploading in the GEE platform the georeferenced tiff aereal images of a new area of interest and by changing the parameters in line 17 of the script provided.
### Railways shapefile
pre-processing:
### Buffer zones
The classification was performed more than once, by filtering the mosaic with buffers around the railway infrastructure of different width. It is possible to change this parameter at line 50 of the script, that is actually setted as 10 meters. 
### Groun Control Points

### Supervised classification
The Classifier package handles supervised classification by traditional ML algorithms

### Contact
For any queries, contact sofia.prandelli2@unibo.it

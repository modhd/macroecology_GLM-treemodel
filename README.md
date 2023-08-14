# Species Distribution Model for five European tree species
A simple Species Distribution Model (SDM) applying a generalized linear model (GLM) in the R environment. Point occurences of five European tree species are linked to spatial data of three environmental variables (annuam mean temperature plus precipitation, elevation above sealevel). \
Some conceptual problems are present due to computational reasons: (i) the area limiting shapefile (retreived form the coastline shapefile) is larger than the area covered by the actual point occurences, leading to false peseudo-absences; (ii) the pseudo-absences were generated bbased on grid centroids, but centroid resolution was 100x100 km insteadt of 1x1 km given in the point occurence dataset.

# Datasets
**Point occurences** \
Mauri, A., Strona, G. & San-Miguel-Ayanz, J. (2017). EU-Forest, a high-resolution tree occurrence dataset for Europe. Scientific data, 4, 1–8.\
**Precipitation and temperature data** Must be added to the "Results" folder\
[Chelsa - Climate data http://chelsa-climate.org/]
Temperature data: [https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/climatologies/bio/CHELSA_bio10_01.tif]
Precipitation data: [https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/climatologies/bio/CHELSA_bio10_12.tif]


**Elevation data**\
GMTED - Digital elevation data [https://topotools.cr.usgs.gov/gmted_viewer/viewer.html] 

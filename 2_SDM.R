library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)
library(raster)
library(performance)
library(patchwork)


# 1. Read data----
## Read occurence data and filter for top 5 most abundant species----
spec_occ <- read.csv("data/EUForestspecies.csv")
head(spec_occ)
top5 <- c("Pinus sylvestris", "Picea abies", "Fagus sylvatica",
          "Quercus robur", "Betula pubescens")
species.5 <- spec_occ[spec_occ$SPECIES.NAME %in% top5, ] |>
  dplyr::select(X, Y, COUNTRY, SPECIES.NAME)
(species.5.sf <- st_as_sf(species.5, coords = c("X", "Y")))
st_crs(species.5.sf) <- 3035
st_bbox(species.5.sf, crs = 3035)

?st_bbox()
#species.5.sf
#head(species.5)
#table(species.5$SPECIES.NAME)


## Read boundary data----
coastline <- st_read("./data/Europe_coastline/Europe_coastline_poly.shp")
st_transform(coastline, crs = 3035)
#plot(coastline, col = "red")
#st_geometry(coastline)
#st_area(coastline)
coastline$area <- st_area(coastline)

# simplifcation: drop african polygon and islands
coastline.rough <- coastline[-1471,] # drop 2nd largest Polygon (African continent)
coastline.rough <- coastline.rough[order(coastline.rough$area, decreasing = T),]
coastline.rough <- coastline.rough[1:10,] # only 10 largest Polygons
coastline.rough <- st_intersection(coastline.rough, species.5.sf)
coastline.rough$area <- units::set_units(coastline.rough$area, "km^2")
head(coastline.rough)
plot(coastline.rough)  

# ggplot(data = coastline.rough) + geom_sf() + 
#   geom_sf(data = species.5.sf) + facet_wrap(~ SPECIES.NAME) +
#   labs(title = "Dataset point occurences") +
#   scale_fill_continuous(name = "Tree species") -> plot.species

## Grid polygon creation----

# Define required functions
gpat_create_grid = function(x, brick = FALSE){
  header = gpat_header_parser(x)

  x1 = header$start_x
  y1 = header$start_y
  x2 = header$start_x + header$res_x * header$n_cols
  y2 = header$start_y + header$res_y * header$n_rows

  single_cell_creator = function(x1, y1, x2, y2){
    list(rbind(c(x1, y1), c(x1, y2), c(x2, y2), c(x2, y1), c(x1, y1)))
  }

  my_bb = single_cell_creator(x1, y1, x2, y2) %>%
    st_polygon() %>%
    st_sfc()

  if (header$proj_4 != ""){
    my_bb = my_bb %>%
      st_set_crs(value = header$proj_4)
  }

  my_grid = gpat_st_make_grid(my_bb,
                              n = c(header$n_cols, header$n_rows),
                              brick = brick)
  my_grid
}

gpat_st_make_grid = function(x,
                             n = c(10, 10),
                             brick = FALSE){

  offset = st_bbox(x)[c(1, 4)]
  bb = st_bbox(x)
  n = rep(n, length.out = 2)
  nx = n[1]
  ny = n[2]
  xc = seq(offset[1], bb[3], length.out = nx + 1)
  yc = seq(offset[2], bb[2], length.out = ny + 1)

  ret = vector("list", nx * ny)
  square = function(x1, y1, x2, y2){
    st_polygon(list(matrix(c(x1, x2, x2, x1, x1, y1, y1, y2, y2, y1), 5)))
  }
  for (i in 1:nx) {
    for (j in 1:ny) {
      ret[[(j - 1) * nx + i]] = square(xc[i], yc[j], xc[i + 1], yc[j + 1])
    }
  }

  my_grid = st_sf(geometry = st_sfc(ret, crs = st_crs(x)))

  if (brick){
    ids_y = rep(c(1, 1, 2, 2), length.out = ny) # rows groups
    ids = numeric(length = nrow(my_grid)) # local ids

    spatial_ids = seq_len(nx %/% 2 + 1) # starting cell numbers in columns
    for (i in seq_len(ny)){
      ids_y_i = ids_y[i] # which row group
      if (ids_y_i == 1){
        ids[seq_len(nx) + (i - 1) * nx] = rep(spatial_ids, each = 2, length.out = nx) # add cell numbers
      } else if (ids_y_i == 2){
        if (ids_y[i-1] != ids_y[i]){ # if row group changed (yes)
          spatial_ids = spatial_ids + (nx %/% 2 + 1) # new cell numbers in columns
          ids[seq_len(nx) + (i - 1) * nx] = c(spatial_ids[1], rep(spatial_ids[-1], each = 2, length.out = nx-1)) # add cell numbers
        } else {
          ids[seq_len(nx) + (i - 1) * nx] = c(spatial_ids[1], rep(spatial_ids[-1], each = 2, length.out = nx-1)) # add cell numbers
          spatial_ids = spatial_ids + (nx %/% 2 + 1) # new cell numbers in columns
        }
      }
    }
    n = c(length(spatial_ids), max(ids) / length(spatial_ids))

    my_grid = aggregate(my_grid, by = list(ids), mean) %>%
      st_cast(to = "POLYGON", warn = FALSE) # aggregate by ids and convert to POLYGON
    my_grid$Group.1 = NULL
  }

  df_ids = create_ids(n[1], n[2])

  my_grid = st_sf(data.frame(my_grid, df_ids))

  return(my_grid)
}

create_ids = function(num_c, num_r){
  m = 1
  m_row = 1
  m_col = 1

  result = data.frame(col = integer(length = num_r * num_c),
                      row = integer(length = num_r * num_c))

  for (i in seq_len(num_r)){
    for (j in seq_len(num_c)){
      result[m, "col"] = m_col
      result[m, "row"] = m_row
      m = m + 1
      m_col = m_col + 1
    }
    m_col = 1
    m_row = m_row + 1
  }
  return(result)
}

# Create grid
#grid <- gpat_st_make_grid(coastline.rough, n = c(100, 100)) |> st_intersection(coastline.rough) 
#st_write(grid, "100x100grid.shp")
grid <- st_read("100x100grid.shp")

## Pseudo-absences----
grid_centroids <- st_centroid(grid)
grid_centroids$present <- 0
grid_centroids <- st_transform(grid_centroids, crs = 3035)

ggplot() + geom_sf(data = coastline.rough) +
  geom_sf(shape = 3, data = st_geometry(grid_centroids)) +
  theme(text = element_text(size = 18)) +
  labs(title = "Generated grid centroids") -> plot.centroids
#plot.centroids


## Pinus sylvestris
pin_syl <- st_as_sf(
  species.5[which(species.5.sf$SPECIES.NAME == "Pinus sylvestris"),
            c("X", "Y")],
  coords = c("X", "Y"),
  crs = 3035)
pin_syl$present <- 1

# add pseudo-absences
set.seed(100)
pin_syl_abs <- rbind(grid_centroids[sample(nrow(grid_centroids),
                                           nrow(pin_syl), replace = TRUE),
                                    "present"], pin_syl) |> st_transform(4326)

#table(pin_syl_abs$present)

## Fagus sylvatica
fag_syl <- st_as_sf(
  species.5[which(species.5$SPECIES.NAME == "Fagus sylvatica"),
            c("X", "Y")],
  coords = c("X", "Y"),
  crs = 3035)
fag_syl$present <- 1

# add pseudo-absences
set.seed(100)
fag_syl_abs <- rbind(grid_centroids[sample(nrow(grid_centroids),
                                           nrow(fag_syl), replace = TRUE),
                                  "present"], fag_syl) |> st_transform(4326)


## Picea abies
pic_abi <- st_as_sf(
  species.5[which(species.5$SPECIES.NAME == "Picea abies"),
            c("X", "Y")],
  coords = c("X", "Y"),
  crs = 3035)
pic_abi$present <- 1

# add pseudo-absences
set.seed(100)
pic_abi_abs <- rbind(grid_centroids[sample(nrow(grid_centroids),
                                           nrow(pic_abi), replace = TRUE),
                                    "present"], pic_abi)  |> st_transform(4326)

## Quercus robur
que_rob <- st_as_sf(
  species.5[which(species.5$SPECIES.NAME == "Quercus robur"),
            c("X", "Y")],
  coords = c("X", "Y"),
  crs = 3035)
que_rob$present <- 1
# add pseudo-absences
set.seed(100)
que_rob_abs <- rbind(grid_centroids[sample(nrow(grid_centroids),
                                           nrow(que_rob), replace = TRUE),
                                    "present"], que_rob)  |> st_transform(4326)
## Betula pubescens
bet_pub <- st_as_sf(
  species.5[which(species.5$SPECIES.NAME == "Betula pubescens"),
            c("X", "Y")],
  coords = c("X", "Y"),
  crs = 3035)
bet_pub$present <- 1
# add pseudo-absences
set.seed(100)
bet_pub_abs <- rbind(grid_centroids[sample(nrow(grid_centroids),
                                           nrow(bet_pub), replace = TRUE),
                                    "present"], bet_pub)  |> st_transform(4326)

pseudo_abs <- c("pin_syl", "pic_abi", "fag_syl", "que_rob", "bet_pub")

# 2. Environmental data----
# align the crs
species.5.sf <- st_transform(species.5.sf, 4326)
coastline.rough <- st_transform(coastline.rough, 4326)
grid84 <- st_transform(grid, 4326)

## Precipitation data----
prec <- rast("data/CHELSA_bio10_12.tif")
prec <- crop(prec, coastline.rough) 
prec <- mask(prec, coastline.rough)
prec <- aggregate(prec, fact = 4, fun = mean)

ggplot() + geom_spatraster(data = prec) + geom_sf(data = coastline.rough, fill = NA) + 
  labs(fill = "[mm]", title = "Annual mean precipitation data") +
  theme(text = element_text(size = 18)) +
  scale_fill_viridis_c(direction = -1, na.value = NA) -> plot.prec
#plot.prec

## Elevation data ----
elev <- rast("data/mn30_grd/") + 0 # treat as numeric
elev <- crop(elev, coastline.rough)
elev <- mask(elev, coastline.rough)
elev <- aggregate(elev, fact = 4, fun = mean)

ggplot() + geom_spatraster(data = elev) + geom_sf(data = coastline.rough, fill = NA) + 
  labs(fill = "[m]", title = "Elevation above sealevel data") +
  theme(text = element_text(size = 18)) +
  scale_fill_viridis_c(option = "inferno", direction = -1, na.value = NA) -> plot.elev
# plot.elev

## Temperature data----
temp <- rast("data/CHELSA_bio10_01.tif")/10 
temp <- crop(temp, coastline.rough)
temp <- mask(temp, coastline.rough)
temp <- aggregate(temp, fact = 4, fun = mean)

ggplot() + geom_spatraster(data = temp) + geom_sf(data = coastline.rough, fill = NA) + 
  labs(fill = "[°C]", title = "Annual mean temperature data") +
  theme(text = element_text(size = 18)) +
  scale_fill_viridis_c(option = "plasma", na.value = NA) -> plot.temp
#plot.temp

#plot.prec + plot.elev + plot.temp


##Extract env. data----
  pin_syl_abs$temp <- extract(temp,pin_syl_abs)[, "CHELSA_bio10_01"]
  pin_syl_abs$elev <- extract(elev, pin_syl_abs)[, "mn30_grd"]
  pin_syl_abs$prec <- extract(prec,pin_syl_abs)[, "CHELSA_bio10_12"]
  pin_syl_abs <- na.omit(pin_syl_abs)
  
  pic_abi_abs$elev <- extract(elev, pic_abi_abs)[, "mn30_grd"]
  pic_abi_abs$temp <- extract(temp,pic_abi_abs)[, "CHELSA_bio10_01"]
  pic_abi_abs$prec <- extract(prec,pic_abi_abs)[, "CHELSA_bio10_12"]
  pic_abi_abs <- na.omit(pic_abi_abs)
  
  fag_syl_abs$elev <- extract(elev, fag_syl_abs)[, "mn30_grd"]
  fag_syl_abs$temp <- extract(temp,fag_syl_abs)[, "CHELSA_bio10_01"]
  fag_syl_abs$prec <- extract(prec,fag_syl_abs)[, "CHELSA_bio10_12"]
  fag_syl_abs <- na.omit(fag_syl_abs)
  
  que_rob_abs$elev <- extract(elev, que_rob_abs)[, "mn30_grd"]
  que_rob_abs$temp <- extract(temp,que_rob_abs)[, "CHELSA_bio10_01"]
  que_rob_abs$prec <- extract(prec,que_rob_abs)[, "CHELSA_bio10_12"]
  que_rob_abs <- na.omit(que_rob_abs)
  
  bet_pub_abs$elev <- extract(elev, bet_pub_abs)[, "mn30_grd"]
  bet_pub_abs$temp <- extract(temp,bet_pub_abs)[, "CHELSA_bio10_01"]
  bet_pub_abs$prec <- extract(prec,bet_pub_abs)[, "CHELSA_bio10_12"]
  bet_pub_abs <- na.omit(bet_pub_abs)


# 3. Logical regression models----
pseudo_abs <- c("pin_syl", "pic_abi", "fag_syl", "que_rob", "bet_pub")

performance <- data.frame(model = pseudo_abs,
                          RMSE = 1:5,
                          R2_adj = 1:5)

for (i in 1:5) {
  temp.name <- paste("sdm", pseudo_abs[i], sep = "_")
  
  assign(temp.name, 
         glm(present ~ elev + prec + temp, data = get(paste(pseudo_abs[i], "abs", sep = "_")),
                 family = "binomial"))
  print(paste("Model from", temp.name, ":", sep = " "))
  print(summary(get(temp.name)))
  print(performance::model_performance(get(temp.name), metrics = c("RMSE", "R2_adj")))
  
  performance$RMSE[i] <- performance::model_performance(get(temp.name),
                                                       metrics = c("RMSE", "R2_adj"))["RMSE"]
  performance$R2_adj[i] <- performance::model_performance(get(temp.name), 
                                                          metrics = c("RMSE", "R2_adj"))["R2_Tjur"]

}

performance


#4. Predictions----

## New data----
new_data <- data.frame(elev = as.numeric(values(elev)),
                       temp = as.numeric(values(temp)),
                       prec = as.numeric(values(prec)))

## Predictions loop----
for (i in 1:5) {

  temp.mod <- paste("sdm", pseudo_abs[i], sep = "_")
  temp.pred <- paste(pseudo_abs[i], "predict", sep = "_")
  
  # predictions
  assign(temp.pred,
         predict(get(temp.mod), newdata = new_data, type = "response"))

  # create rasters
  assign(paste(pseudo_abs[i], "raster", sep = "_"),
         rast(ext(coastline.rough), nrows = nrow(elev), ncols = ncol(elev),
              vals = get(temp.pred)))

}

## AUCs----

sdm_auc <- data.frame(model = pseudo_abs,
                      auc = 1:5)

for (i in 1:5) {
  temp.mod <- paste("sdm", pseudo_abs[i], sep = "_")
  temp.pred <- paste(pseudo_abs[i], "predict", sep = "_")
  
  # calc AUC
  p <- predict(get(temp.mod), get(paste(pseudo_abs[i], "abs", sep = "_")), 
               type = "response")
  pr <- ROCR::prediction(p, get(paste(pseudo_abs[i], "abs", sep = "_"))$present)
  auc <- ROCR::performance(pr, measure = "auc")
  sdm_auc$auc[i] <- auc@y.values[[1]]
}
sdm_auc

# 5. Comparative plots----

## Betula pubescens
ggplot() + geom_spatraster(data = bet_pub_raster) +
  geom_sf(data = coastline.rough, fill = NA) +
  labs(title = "Betula pubescens predictions", fill = "Occurence probability") +
  theme(text = element_text(size = 8)) +
    scale_fill_gradientn(colours = terrain.colors(10, rev = T), na.value = NA) -> bet_pub.plot2

ggplot() + geom_sf(data = coastline.rough) +
  geom_sf(color = "darkolivegreen3", 
          data = species.5.sf[which(species.5.sf$SPECIES.NAME == "Betula pubescens"),], 
          alpha = 0.2) +
  theme(text = element_text(size = 8)) +
  labs(title = "Betula pubescens point occurences") -> bet_pub.plot1

# bet_pub.plot1 + bet_pub.plot2

## Fagus sylvatica
ggplot() + geom_spatraster(data = fag_syl_raster) +
  geom_sf(data = coastline.rough, fill = NA) +
  theme(text = element_text(size = 8)) +
  labs(title = "Fagus sylvatica predictions", fill = "Occurence probability") +
  scale_fill_gradientn(colours = terrain.colors(10, rev = T), na.value = NA) -> fag_syl.plot2
#fag_syl.plot2

ggplot() + geom_sf(data = coastline.rough) +
  geom_sf(color = "darkolivegreen3", 
          data = species.5.sf[which(species.5.sf$SPECIES.NAME == "Fagus sylvatica"),], 
          alpha = 0.2) +
  theme(text = element_text(size = 8)) +
  labs(title = "Fagus sylvatica point occurences") -> fag_syl.plot1

#fag_syl.plot1 + fag_syl.plot2

## Picea abies
ggplot() + geom_spatraster(data = pic_abi_raster) +
  geom_sf(data = coastline.rough, fill = NA) +
  labs(title = "Picea abies predictions", fill = "Occurence probability") +
  theme(text = element_text(size = 8)) +
  scale_fill_gradientn(colours = terrain.colors(10, rev = T), na.value = NA) -> pic_abi.plot2

ggplot() + geom_sf(data = coastline.rough) +
  geom_sf(color = "darkolivegreen3", 
          data = species.5.sf[which(species.5.sf$SPECIES.NAME == "Picea abies"),], 
          alpha = 0.2) +
  theme(text = element_text(size = 8)) +
  labs(title = "Picea abies point occurences") -> pic_abi.plot1

# pic_abi.plot1 + pic_abi.plot2

## Pinus sylvestris
ggplot() + geom_spatraster(data = pin_syl_raster) +
  geom_sf(data = coastline.rough, fill = NA) +
  labs(title = "Pinus sylvestris predictions", fill = "Occurence probability") +
  scale_fill_gradientn(colours = terrain.colors(10, rev = T), na.value = NA) -> pin_syl.plot2

ggplot() + geom_sf(data = coastline.rough) +
  geom_sf(color = "darkolivegreen3", 
          data = species.5.sf[which(species.5.sf$SPECIES.NAME == "Pinus sylvestris"),], 
          alpha = 0.2) +
  theme(text = element_text(size = 8)) +
  labs(title = "Pinus sylvestris point occurences") -> pin_syl.plot1

# pin_syl.plot1 + pin_syl.plot2

## Quercus robur
ggplot() + geom_spatraster(data = que_rob_raster) +
  geom_sf(data = coastline.rough, fill = NA) +
  labs(title = "Quercus robur predictions", fill = "Occurence probability") +
  theme(text = element_text(size = 8)) +
  scale_fill_gradientn(colours = terrain.colors(10, rev = T), na.value = NA) -> que_rob.plot2

ggplot() + geom_sf(data = coastline.rough) +
  geom_sf(color = "darkolivegreen3", 
          data = species.5.sf[which(species.5.sf$SPECIES.NAME == "Quercus robur"),], 
          alpha = 0.2) +
  theme(text = element_text(size = 8)) +
  labs(title = "Quercus robur point occurences") -> que_rob.plot1

#que_rob.plot1 + que_rob.plot2

bet_pub.plot1 + bet_pub.plot2 +
  fag_syl.plot1 + fag_syl.plot2 +
  pic_abi.plot1 + pic_abi.plot2 +
  pin_syl.plot1 + pin_syl.plot2 +
  que_rob.plot1 + que_rob.plot2 +
  plot_layout(nrow = 5) +
  plot_annotation("Occurences vs. SDM predictions") -> plot.summary

#plot.summary




# Dump----
# Coastline data from Eurostat
# https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units/nuts#nuts21
grid <- st_read("data/NUTS_RG_10M_2021_3035.shp")

# Grid creation)
##(A)
# grid of 100 x 100 km
"Problem here: point geometry, hardöy convertable into Polygons. Thus, another approach"
# grid <- coastline.rough |> # caution comutional power!
#   st_make_grid(cellsize = c(100000, 100000), what = "centers")
#   st_intersection(coastline.rough)

# grid10 <- coastline.rough |> 
#   st_make_grid(cellsize = c(10000, 10000), what = "centers") |> 
#   st_intersection(coastline.rough)

# grid100 <- st_read("./data/100x100_grid/grid100x100.shp")
# grid100 <- st_transform(grid100, crs = st_crs(coastline.rough)) # not working somehow, too largely scaled
# st_geometry(grid100)

##(B)
"Try 3"
# grid.temp <- grid
# st_crs(grid.temp)
#  grid.temp$grid_id <- 1:length(grid)
# grid
# grid.temp
# sf::st_write(grid, paste0(tempdir(), "/", "grid100x100.shp"))

#plot(grid, pch = 3, main = "100x100 km^2 grid") 
#m plot(coastline.rough, add = T, col = NA)

# (grid.poly <- grid)
# grid.poly <- as.data.frame(od::sfc_point_to_matrix(grid.poly)) |> rename(X = V1, Y = V2)
# #grid.poly <- 
# grid::grid.polygon(x = grid.poly$X, y = grid.poly$Y) |> plot()
# #sfheaders::sfc_polygon(grid.poly, x = "X", y = "Y")
# plot(grid.poly)
# grid.poly <- sfheaders::st_sfc(grid, "MULTYPOLYGON") 



# Make prediction grid overlapping strata and NAFO divisions
library(sf)
library(stars)

strata_shp <- read_sf("shapefiles/2HJ3KLNOP_Strata_Polygons_WGS84.shp") %>%
  st_transform(crs = 32622)

grid_all <- sf::st_make_grid(strata_shp, cellsize = c(5000, 5000), what = "polygons")

ind_all <- st_within(grid_all, strata_shp) %>%
  sapply(function(i) if (!length(i)) NA else i[1])

grid <- grid_all[!is.na(ind_all)] %>%
  st_as_sf()
grid_coord <- st_centroid(grid) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(ind = ind_all[!is.na(ind_all)],
         strat = strata_shp$STRAT[ind],
         div = strata_shp$DIV[ind],
         area_km2 = 25) %>%
  select(!ind)
#readr::write_csv(grid_coord, file = "shapefiles/prediction_grid.csv")

# Add depth to grid from bathymetry data sets
####
# Got map (ETOPO1 ice) (exportImage.tif) from: https://www.ncei.noaa.gov/maps/bathymetry/
# Coordinates from st_bbox(strata_shp), 60 arc-second resolution
#xmin      ymin      xmax      ymax
#-60.95551  42.65603 -46.46244  57.71356
bathy <- stars::read_stars("shapefiles/bathymetry_60.tif") %>%
  st_as_sf()

bind_all <- st_centroid(grid) %>%
  st_transform(crs = st_crs(bathy)) %>%
  st_intersects(bathy) %>%
  sapply(function(i) if (!length(i)) NA else i[1])

depth <- bathy$`bathymetry_60.tif`[bind_all]

bathy_grid_coord <- cbind(
  grid_coord,
  depth
)
readr::write_csv(bathy_grid_coord, file = "shapefiles/prediction_grid.csv")

# Figures and summary

g <- ggplot(grid_coord, aes(X, Y, fill = factor(strat))) +
  geom_raster()

g <- ggplot(grid_coord, aes(X, Y, fill = div)) +
  geom_raster()

g <- ggplot(bathy_grid_coord, aes(X, Y, fill = -1 * depth)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log")

summarise(grid_coord, n = n(), .by = strat)

summarise(grid_coord, area = sum(area_km2), .by = strat)
summarise(grid_coord, n = n(), .by = div)



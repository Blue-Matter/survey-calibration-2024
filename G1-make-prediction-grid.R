

library(raster)
library(sf)

# Identify division (for figures, see code below once spatial grid for distribution model is set up)
poly_1CD <- sf::st_read("calibration_ghl/NAFO SA0 Survey Strata/Divisions/Divisions.shp") %>%
  filter(ZONE %in% c("1C", "1D")) %>%
  select(AREA, ZONE, geometry) %>%
  mutate(Code = "All")
st_crs(poly_1CD) <- "WGS84"

poly_0A_off <- sf::st_read(
  "calibration_ghl/NAFO SA0 Survey Strata/NewNAFOstrata2021/NAFO_0A_Offshore_Strata_2021.shp"
) %>%
  filter(Code != "Shallow") %>%
  select(AREA, ZONE, geometry, Code)

poly_0B_off <- sf::st_read(
  "calibration_ghl/NAFO SA0 Survey Strata/NewNAFOstrata2021/NAFO_0B_Offshore_Strata_2021.shp"
) %>%
  filter(Code != "Shallow") %>%
  select(AREA, ZONE, geometry, Code)

poly <- rbind(
  poly_1CD,
  poly_0A_off,
  poly_0B_off
)

poly_utm <- st_transform(poly, 32621)

plot(poly_utm)

grid <- expand.grid(
  X = seq(-80000, 850000, by = 10000),
  Y = seq(6700000, 8700000, by = 10000)
)

grid_sf <- sf::st_as_sf(grid, coords = c("X", "Y"), crs = 32621)
get_div <- sf::st_intersects(grid_sf, poly_utm, sparse = TRUE) %>%
  lapply(function(i) if (length(i) > 1) i[1] else i) %>%
  structure(names = 1:nrow(grid))
div_key <- do.call(c, get_div)

grid$div <- rep(NA_character_, nrow(grid))
grid$div[as.numeric(names(div_key))] <- poly_utm$ZONE[as.numeric(div_key)]

grid_out <- filter(grid, !is.na(div)) %>%
  mutate(area = 10000^2)
readr::write_csv(grid_out, "calibration_ghl/grid_0AB1CD.csv")

g <- ggplot(grid_out, aes(X, Y, fill = div)) +
  geom_sf(data = poly_utm, inherit.aes = FALSE) +
  coord_sf() +
  geom_tile(width = 10000, height = 10000, colour = "grey40")

poly %>%
  as.data.frame() %>%
  summarise(Area = sum(AREA), .by = ZONE)

grid_out %>%
  summarise(area = sum(area)/1e3/1e3, .by = div) %>%
  arrange(div)



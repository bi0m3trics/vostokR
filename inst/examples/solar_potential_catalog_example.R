# Example script for using vostokR package

library(lidR)
library(vostokR)
library(future)

# Read example LiDAR data
ctg <- readLAScatalog("D:/lidar_temp/NAUCampus_2019/")
ctg <- catalog_select(ctg)

opt_chunk_size(ctg) <- 100
opt_chunk_buffer(ctg) <- 150
opt_stop_early(ctg) <- TRUE
opt_wall_to_wall(ctg) <- TRUE
opt_chunk_alignment(ctg) <- c(0,0)

catalog_solarpot <- function(chunk, start_date = NULL, end_date = NULL)
{
  # read the chunk as a LAS
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)

  # Add normal vectors if not present
  las <- add_normals(las)

  # Extract the location and timezone information
  las_info <- extract_crs_info(las)

  # Converts start and end dates to day of year numbers
  time <- vostokR:::date_to_day_numbers("2019/06/01", "2019/06/30", 2019)

  # Calculate solar potential
  # Using Seattle, WA coordinates as an example
  las_solar <- calculate_solar_potential(las,
                                         year = 2019,
                                         day_start = time$day_start,
                                         day_end = time$day_end,
                                         day_step = 30,  # Calculate every 30 days
                                         minute_step = 30,  # 30-minute intervals
                                         min_sun_angle = 5,
                                         voxel_size = 1,
                                         lat = las_info$lat,  # Seattle latitude
                                         lon = las_info$lon,  # Seattle longitude
                                         timezone = las_info$timezone)  # Pacific Time

  # Extracts ground points from solar potential results and converts to terra SpatRaster
  solar_raster <- solar_ground_raster(las_solar, res = 1)

  bbox  <- terra::ext(chunk)
  solar_raster <- terra::crop(solar_raster, bbox)
  return(solar_raster)
}

plan(multisession, workers = parallel::detectCores(logical=FALSE)-1)
set_lidr_threads(1L)

output <- catalog_apply(ctg, catalog_solarpot)
output = do.call(terra::merge, output)

# After calculating solar potential
plot(output, main = "Ground Solar Potential", col = heat.colors(100))

# =============================================================================
library(leaflet)
library(terra)

# Get values
vals <- terra::values(output, na.rm = TRUE)

# Define your palette from full sun -> shadow
sun_palette <- c(
  "#FDB813",  # full sun
  "#FFD966",  # moderate
  "#A3C4BC",  # low
  "#2E3440"   # no sun / shadow
)

# Make a numeric color palette
pal <- colorNumeric(
  palette = sun_palette,
  domain = vals
)

# Example use in leaflet
leaflet() %>%
  addTiles() %>%
  addRasterImage(output, colors = pal) %>%
  addLegend(pal = pal, values = vals, title = "Solar Potential")


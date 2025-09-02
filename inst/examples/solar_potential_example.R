# Example script for using vostokR package

library(lidR)
library(vostokR)

# Read example LiDAR data
LASfile <- system.file("extdata", "test.laz", package="vostokR")
las <- readLAS(LASfile)

# Add normal vectors if not present
las <- add_normals(las)

# Calculate solar potential
# Using Seattle, WA coordinates as an example
las_solar <- calculate_solar_potential(las,
                                    year = 2025,
                                    day_start = 1,
                                    day_end = 365,
                                    day_step = 30,  # Calculate every 30 days
                                    minute_step = 30,  # 30-minute intervals
                                    min_sun_angle = 5,
                                    voxel_size = 1,
                                    lat = 47.6062,  # Seattle latitude
                                    lon = -122.3321,  # Seattle longitude
                                    timezone = -8)  # Pacific Time

# Plot the results as a raster
plot(las_solar, color = "solar_potential", pal = heat.colors(100))

# Print summary statistics
summary(las_solar$solar_potential)

# After calculating solar potential
solar_raster <- solar_ground_raster(las_solar, res = 1)# Get values
vals <- terra::values(solar_raster, na.rm = TRUE)

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

plot(solar_raster, main = "Ground Solar Potential", col = pal)

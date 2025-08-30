# Example script for using vostokR package

library(lidR)
library(vostokR)

# Read example LiDAR data
LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
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
solar_raster <- plot_solar_raster(las_solar)

# Print summary statistics
summary(las_solar$solar_potential)

# ================================================================
# Solar Potential Workflow with vostokR
# Combines: single LAS example + catalog example + custom sun palette
# ================================================================

# ---- Load required libraries ----
library(lidR)      # For LiDAR data reading and catalog processing
library(vostokR)   # Solar potential calculations (VOSTOK interface)
library(terra)     # Raster handling (SpatRaster)
library(leaflet)   # Interactive map visualization
library(future)    # Parallel processing (for catalog operations)

# ---- Define custom sun palette ----
sun_palette <- c("#414C65", # No sun
                   "#4A587A",
                   "#6F6256", # low sun (muted cool gray-green)
                   "#A39568",
                   "#C29E25", # moderate sun (pale yellow)
                   "#C2A64B",
                   "#CEB255", # full sun (bright golden yellow)
                   "#E6C65C",
                   "#FFF4CE"  # Full sun
)

# Create a color ramp (100 steps) for terra::plot()
sun_ramp <- colorRampPalette(sun_palette)

# ================================================================
# A) SINGLE LAS EXAMPLE
# ================================================================

# ---- Load example LAS ----
# Replace with your own LAS file path as needed.
LASfile <- system.file("extdata", "test.laz", package = "vostokR")
las <- readLAS(LASfile)
if (is.empty(las)) stop("LAS is empty.")  # Safety check

# ---- Add normal vectors ----
# Required for solar potential calculation (determines surface orientation)
las <- add_normals(las)

# ---- Calculate solar potential ----
# Full-year, every 30 days, 30-minute steps
# Using Seattle, WA as the location (adjust lat/lon/timezone for your site)
las_solar <- calculate_solar_potential(
  las,
  year          = 2019,
  day_start     = 121,     # Start on day 1 (Jan 1)
  day_end       = 212,   # End on day 365 (Dec 31)
  day_step      = 1,    # Simulate every 30 days
  minute_step   = 60,    # Solar position every 30 minutes
  min_sun_angle = 15,     # Minimum sun angle considered (degrees)
  voxel_size    = 1,     # Voxel size (meters)
  lat           = 35.27041,   # Latitude
  lon           = -111.6905, # Longitude
  timezone      = -7         # UTC offset (Pacific Time)
)

plot(las_solar, color = "solar_potential", pal = sun_ramp(100))

# ---- Convert to ground raster ----
# Aggregates solar potential to ground surface
solar_raster_single <- solar_ground_raster(las_solar, res = 1)

# ---- Fill holes (NA cells) with neighborhood mean ----
# Define a 3x3 moving window (adjust size if you want smoother fill)
window <- matrix(1, nrow = 3, ncol = 3)

# Compute local mean, ignoring NA
solar_raster_single <- focal(solar_raster_single, w = window, fun = mean, na.policy = "only", na.rm = TRUE)

# ---- Visualize as static raster ----
plot(
  solar_raster_single,
  main = "Ground Solar Potential (Single LAS)",
  col  = sun_ramp(100)  # Use custom palette ramp
)

# ---- Print summary statistics ----
# Helpful for quick check of value ranges
print(summary(values(solar_raster_single, na.rm = TRUE)))

# ================================================================
# B) CATALOG (MULTI-TILE) EXAMPLE
# ================================================================

# ---- Read catalog ----
# Replace with your directory containing multiple LAS/LAZ tiles
# ctg <- readLAScatalog("D:/lidar_temp/NAUCampus_2019/")
ctg <- "D:/lidar_temp/NAUCampus_2019/USGS_LPC_AZ_Coconino_2019_B19_w1408n1465.laz"
ctg <- readLAScatalog(ctg) # Make sure files are selected

# ---- Configure catalog options ----
opt_chunk_size(ctg)      <- 100    # Size of processing chunks (meters)
opt_chunk_buffer(ctg)    <- 150    # Buffer size to avoid edge artifacts
opt_stop_early(ctg)      <- TRUE   # Stops processing if no data in chunk
opt_wall_to_wall(ctg)    <- TRUE   # Forces full coverage processing
opt_chunk_alignment(ctg) <- c(0, 0) # Aligns chunks to coordinate grid

# ---- Define per-chunk processing function ----
catalog_solarpot <- function(chunk,
                             start_date = "2019/01/01",
                             end_date   = "2019/12/31") {

  # Read chunk
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)  # Skip if no points

  # Add normals (required for solar calculation)
  las <- add_normals(las)

  # Extract CRS-based location info (lat, lon, timezone)
  info <- extract_crs_info(las)

  # Convert date range to day-of-year numbers
  dnums <- vostokR:::date_to_day_numbers(start_date, end_date, 2019)

  # Calculate solar potential for this chunk
  las_solar <- calculate_solar_potential(
    las,
    year          = 2019,
    day_start     = dnums$day_start,
    day_end       = dnums$day_end,
    day_step      = 1,
    minute_step   = 60,
    min_sun_angle = 15,
    voxel_size    = 1,
    lat           = info$lat,
    lon           = info$lon,
    timezone      = info$timezone
  )

  # Convert solar potential to ground raster
  r <- solar_ground_raster(las_solar, res = 1)

  # Crop raster to chunk extent (avoids overhangs from buffer)
  r <- terra::crop(r, terra::ext(chunk))
  return(r)
}

# ---- Parallelization settings ----
plan(multisession, workers = max(1, parallel::detectCores(logical = FALSE) - 1))
set_lidr_threads(1L)  # Use single thread per chunk (parallel handled by future)

# ---- Apply to catalog ----
# Produces a list of SpatRasters
rasters <- catalog_apply(ctg, catalog_solarpot)

# ---- Merge chunk rasters into one seamless raster ----
output <- do.call(terra::merge, rasters)

# ---- Fill holes (NA cells) with neighborhood mean ----
# Compute local mean, ignoring NA
output <- focal(output, w = window, fun = mean, na.policy = "only", na.rm = TRUE)

# ---- Plot merged raster ----
plot(
  output,
  main = "Ground Solar Potential (Catalog)",
  col  = sun_ramp(100)
)

# ================================================================
# C) INTERACTIVE LEAFLET MAP
# ================================================================

# ---- Prepare color palette for leaflet ----
library(leaflet)
library(terra)

# Values + palette
vals <- terra::values(output, na.rm = TRUE)
pal <- colorNumeric(
  palette  = sun_palette,
  domain   = vals,
  na.color = "transparent"
)

# Optional: keep raster above basemaps
leaflet(options = leafletOptions(preferCanvas = TRUE)) %>%
  # --- Base layers ---
  addProviderTiles(providers$OpenTopoMap,         group = "OpenTopoMap") %>%
  addProviderTiles(providers$Esri.WorldImagery,   group = "Esri World Imagery") %>%
  addProviderTiles(providers$Esri.WorldShadedRelief, group = "Shaded Relief") %>%  # hillshade
  # --- Overlay: solar raster ---
  addRasterImage(output, colors = pal, project = TRUE, opacity = 0.5,
                 group = "Solar Potential") %>%
  addLegend(pal = pal, values = vals, title = "Solar Potential") %>%
  # --- Layer switcher ---
  addLayersControl(
    baseGroups = c("OpenTopoMap", "Esri World Imagery", "Shaded Relief"),
    overlayGroups = c("Solar Potential"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  # Start with imagery as the active base (order matters: last added base is not auto-selected),
  # so we explicitly show/hide groups:
  hideGroup("OpenTopoMap") %>%
  hideGroup("Shaded Relief")

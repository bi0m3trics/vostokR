#' Calculate Solar Potential for LiDAR Point Cloud
#'
#' This function takes a LiDAR point cloud from the lidR package and calculates
#' the solar potential for each point, taking into account shadowing effects from
#' surrounding points. Location information (lat/lon/timezone) is auto-detected 
#' from the CRS when possible.
#'
#' @param las LAS or LAScatalog object from lidR package containing point cloud data
#' @param year Numeric. Year for solar calculation (default: 2025)
#' @param day_start Numeric. Start day of year (1-365). Ignored if start_date is provided.
#' @param day_end Numeric. End day of year (1-365). Ignored if end_date is provided.
#' @param start_date Date or character. Start date (e.g., "2025-06-01"). Overrides day_start.
#' @param end_date Date or character. End date (e.g., "2025-08-31"). Overrides day_end.
#' @param day_step Numeric. Step size in days for calculation (default: 30)
#' @param minute_step Numeric. Time step in minutes for calculation within each day (default: 30)
#' @param min_sun_angle Numeric. Minimum sun angle in degrees for calculation (default: 5)
#' @param voxel_size Numeric. Size of voxels for octree shadow calculation (default: 1)
#' @param lat Numeric. Latitude of the study area (auto-detected from CRS if NULL)
#' @param lon Numeric. Longitude of the study area (auto-detected from CRS if NULL)
#' @param timezone Numeric. Time zone offset from UTC (auto-detected from CRS if NULL)
#' @param ... Additional arguments passed to other methods
#' @return LAS object with added column for solar potential in Wh/m²/day
#' @export
#' @examples
#' \dontrun{
#' library(lidR)
#' library(vostokR)
#' library(terra)
#' 
#' # Single LAS file example using included test data
#' LASfile <- system.file("extdata", "test.laz", package="vostokR")
#' las <- readLAS(LASfile)
#' 
#' # Calculate normals if not present
#' las <- add_normals(las, k = 10)
#' 
#' # Example 1: Auto-detect location from CRS
#' las_solar <- calculate_solar_potential(las, 
#'                                      year = 2025,
#'                                      start_date = "2025-06-01",
#'                                      end_date = "2025-08-31")
#' 
#' # Example 2: Manual location specification  
#' las_solar <- calculate_solar_potential(las, 
#'                                      year = 2025,
#'                                      day_start = 150,
#'                                      day_end = 250,
#'                                      lat = 33.4484,
#'                                      lon = -112.0740,
#'                                      timezone = -7)
#' 
#' # Plot solar potential on point cloud
#' plot(las_solar, color = "solar_potential", pal = heat.colors(100))
#' 
#' # Create ground raster and plot
#' ground_raster <- solar_ground_raster(las_solar, res = 1)
#' plot(ground_raster, main = "Ground Solar Potential")
#' }
#' Extract Geographic Information from LAS CRS
#'
#' Attempts to extract latitude, longitude, and timezone from LAS coordinate reference system
#'
#' @param las LAS object
#' @return List with lat, lon, timezone (or NA if cannot extract)
#' @keywords internal
extract_crs_info <- function(las) {
    crs_info <- list(lat = NA, lon = NA, timezone = NA)
    
    tryCatch({
        # Get CRS from LAS object
        crs_obj <- sf::st_crs(las)
        
        if (!is.na(crs_obj)) {
            # Get center of bounding box
            bbox <- sf::st_bbox(las)
            center_x <- mean(c(bbox["xmin"], bbox["xmax"]))
            center_y <- mean(c(bbox["ymin"], bbox["ymax"]))
            
            # Convert to WGS84 if not already
            if (sf::st_is_longlat(crs_obj)) {
                # Already in geographic coordinates
                crs_info$lon <- center_x
                crs_info$lat <- center_y
            } else {
                # Transform to WGS84
                center_point <- sf::st_sfc(sf::st_point(c(center_x, center_y)), crs = crs_obj)
                center_wgs84 <- sf::st_transform(center_point, crs = 4326)
                coords <- sf::st_coordinates(center_wgs84)
                crs_info$lon <- coords[1]
                crs_info$lat <- coords[2]
            }
            
            # Estimate timezone from longitude (rough approximation)
            crs_info$timezone <- round(crs_info$lon / 15)
        }
    }, error = function(e) {
        # Return NA values if extraction fails
    })
    
    return(crs_info)
}

#' Convert Date Range to Day Numbers
#'
#' Converts start and end dates to day of year numbers
#'
#' @param start_date Date or character. Start date (e.g., "2025-06-01" or as.Date("2025-06-01"))
#' @param end_date Date or character. End date (e.g., "2025-08-31" or as.Date("2025-08-31"))
#' @param year Numeric. Year for the calculation
#' @return List with day_start and day_end
#' @keywords internal
date_to_day_numbers <- function(start_date, end_date, year) {
    # Convert to Date objects if not already
    if (is.character(start_date)) start_date <- as.Date(start_date)
    if (is.character(end_date)) end_date <- as.Date(end_date)
    
    # Ensure dates are in the correct year
    start_date <- as.Date(paste(year, format(start_date, "%m-%d"), sep = "-"))
    end_date <- as.Date(paste(year, format(end_date, "%m-%d"), sep = "-"))
    
    # Convert to day of year
    day_start <- as.numeric(format(start_date, "%j"))
    day_end <- as.numeric(format(end_date, "%j"))
    
    return(list(day_start = day_start, day_end = day_end))
}

#' @export
calculate_solar_potential <- function(las, ...) {
  UseMethod("calculate_solar_potential", las)
}

#' @export
calculate_solar_potential.LAS <- function(las,
                                       year = 2025,
                                       day_start = NULL,
                                       day_end = NULL,
                                       start_date = NULL,
                                       end_date = NULL,
                                       day_step = 30,
                                       minute_step = 30,
                                       min_sun_angle = 5,
                                       voxel_size = 1,
                                       lat = NULL,
                                       lon = NULL,
                                       timezone = NULL,
                                       ...) {
    
    # Auto-detect CRS information if not provided
    if (is.null(lat) || is.null(lon) || is.null(timezone)) {
        crs_info <- extract_crs_info(las)
        
        if (is.null(lat)) lat <- crs_info$lat
        if (is.null(lon)) lon <- crs_info$lon  
        if (is.null(timezone)) timezone <- crs_info$timezone
        
        # Check if we still have missing values
        if (is.na(lat) || is.na(lon) || is.na(timezone)) {
            stop("Could not extract lat/lon/timezone from CRS. Please provide these parameters explicitly.")
        }
        
        message("Auto-detected from CRS: lat=", round(lat, 4), ", lon=", round(lon, 4), ", timezone=", timezone)
    }
    
    # Handle date range conversion
    if (!is.null(start_date) && !is.null(end_date)) {
        if (!is.null(day_start) || !is.null(day_end)) {
            warning("Both date range and day numbers provided. Using date range.")
        }
        day_nums <- date_to_day_numbers(start_date, end_date, year)
        day_start <- day_nums$day_start
        day_end <- day_nums$day_end
        message("Date range ", start_date, " to ", end_date, " converted to days ", day_start, "-", day_end)
    } else {
        # Use default values if not provided
        if (is.null(day_start)) day_start <- 1
        if (is.null(day_end)) day_end <- 365
    }
    
    # Check if normals exist in the point cloud
    if (!all(c("nx", "ny", "nz") %in% names(las@data))) {
        stop("Point cloud must contain normal vectors (nx, ny, nz). Use add_normals() first.")
    }

    # Extract coordinates and normals
    coords <- as.matrix(data.table::data.table(
        X = las@data$X,
        Y = las@data$Y,
        Z = las@data$Z
    ))

    normals <- as.matrix(data.table::data.table(
        nx = las@data$nx,
        ny = las@data$ny,
        nz = las@data$nz
    ))

    # Defensive checks
    if (nrow(coords) == 0) stop("No points in LAS object.")
    if (nrow(coords) != nrow(normals)) stop("Mismatch between coordinates and normals.")

    # Call Rcpp function for solar potential calculation
    tryCatch({
        solar_potential <- calculate_solar_potential_cpp(
            coords = coords,
            normals = normals,
            year = as.integer(year),
            day_start = as.integer(day_start),
            day_end = as.integer(day_end),
            day_step = as.integer(day_step),
            minute_step = as.integer(minute_step),
            min_sun_angle = as.numeric(min_sun_angle),
            voxel_size = as.numeric(voxel_size),
            lat = as.numeric(lat),
            lon = as.numeric(lon),
            timezone = as.numeric(timezone)
        )
    }, error = function(e) {
        stop("Error in native solar potential calculation: ", e$message)
    })

    # Add solar potential to LAS object
    las@data$solar_potential <- solar_potential

    return(las)
}

#' @export
calculate_solar_potential.LAScatalog <- function(las,
                                              year = 2025,
                                              day_start = 1,
                                              day_end = 365,
                                              day_step = 30,
                                              minute_step = 30,
                                              min_sun_angle = 5,
                                              voxel_size = 1,
                                              lat,
                                              lon,
                                              timezone,
                                              ...) {
    # Set processing options for the LAScatalog
    opt_chunk_buffer(las) <- voxel_size * 3  # Buffer size based on voxel size
    opt_chunk_size(las) <- 250  # Process in 250x250 chunks
    opt_output_files(las) <- paste0(tempfile(), "_solar_{ID}")
    
    # Create processing engine function
    solar_engine <- function(chunk) {
        if (is.empty(chunk)) return(NULL)
        
        # Calculate lat/lon for this chunk if not provided
        if (missing(lat) || missing(lon)) {
            chunk_center <- st_coordinates(st_centroid(chunk@bbox))
            chunk_lat <- chunk_center[2]
            chunk_lon <- chunk_center[1]
        } else {
            chunk_lat <- lat
            chunk_lon <- lon
        }
        
        # Add normals if missing
        if (!all(c("nx", "ny", "nz") %in% names(chunk@data))) {
            chunk <- add_normals(chunk)
        }
        
        # Calculate solar potential for this chunk
        chunk_processed <- calculate_solar_potential.LAS(chunk,
                                                      year = year,
                                                      day_start = day_start,
                                                      day_end = day_end,
                                                      day_step = day_step,
                                                      minute_step = minute_step,
                                                      min_sun_angle = min_sun_angle,
                                                      voxel_size = voxel_size,
                                                      lat = chunk_lat,
                                                      lon = chunk_lon,
                                                      timezone = timezone)
        return(chunk_processed)
    }
    
    # Process the catalog
    output <- catalog_apply(las, solar_engine)
    return(output)
}

#' Plot Solar Potential Point Cloud
#'
#' Plots solar potential values directly on the point cloud using lidR's native plotting
#'
#' @param las LAS object with solar potential values
#' @param ... Additional arguments passed to lidR::plot()
#' @export
#' @examples
#' \dontrun{
#' # Plot solar potential on point cloud
#' plot(las_solar, color = "solar_potential", pal = heat.colors(100))
#' }
plot_solar_potential <- function(las, ...) {
    if (!"solar_potential" %in% names(las@data)) {
        stop("LAS object must contain solar_potential values")
    }
    
    # Use lidR's native plot function
    lidR::plot(las, color = "solar_potential", pal = grDevices::heat.colors(100), ...)
}

#' Convert Solar Potential Ground Points to Raster
#'
#' Extracts ground points from solar potential results and converts to terra SpatRaster
#'
#' @param las LAS object with solar potential values
#' @param res Numeric. Resolution of output raster (default: 1)
#' @param ground_class Numeric. Classification code for ground points (default: 2)
#' @param use_all_points Logical. If TRUE and no ground points found, use all points (default: FALSE)
#' @return SpatRaster object from terra
#' @export
#' @examples
#' \dontrun{
#' library(lidR)
#' library(vostokR)
#' library(terra)
#' 
#' # After calculating solar potential
#' solar_raster <- solar_ground_raster(las_solar, res = 0.5)
#' plot(solar_raster, main = "Ground Solar Potential")
#' }
solar_ground_raster <- function(las, res = 1, ground_class = 2, use_all_points = FALSE) {
    if (!"solar_potential" %in% names(las@data)) {
        stop("LAS object must contain solar_potential values")
    }
    
    # Check what classification values are actually present
    class_values <- unique(las@data$Classification)
    cat("Available classification values:", paste(class_values, collapse = ", "), "\n")
    
    # Filter ground points
    ground_points <- lidR::filter_poi(las, Classification == ground_class)
    
    if (nrow(ground_points@data) == 0) {
        if (use_all_points) {
            warning("No ground points found with classification ", ground_class, 
                   ". Using all points instead.")
            ground_points <- las
        } else {
            warning("No ground points found with classification ", ground_class, 
                   ". Try use_all_points=TRUE or check classification values.")
            return(NULL)
        }
    }
    
    cat("Using", nrow(ground_points@data), "points for raster creation\n")
    cat("Solar potential range:", 
        round(range(ground_points@data$solar_potential), 2), "kWh/m²/year\n")
    
    # Create raster using grid_metrics which is better for point attributes
    solar_raster <- lidR::grid_metrics(ground_points, 
                                      ~mean(solar_potential, na.rm = TRUE), 
                                      res = res)
    
    # Convert to terra SpatRaster
    solar_raster <- terra::rast(solar_raster)
    
    # Check if raster has values
    if (all(is.na(terra::values(solar_raster)))) {
        warning("Raster contains no values. Check resolution or point distribution.")
    } else {
        cat("Raster created successfully with", terra::ncell(solar_raster), "cells\n")
        cat("Raster value range:", round(range(terra::values(solar_raster), na.rm = TRUE), 2), "\n")
    }
    
    return(solar_raster)
}

#' @importFrom lidR opt_chunk_buffer opt_chunk_size opt_output_files is.empty st_centroid catalog_apply readLAS rasterize_canopy p2r filter_poi plot grid_metrics
#' @importFrom sf st_coordinates st_crs st_bbox st_is_longlat st_sfc st_point st_transform
#' @importFrom data.table data.table
#' @importFrom methods is
#' @importFrom stats na.omit
#' @importFrom terra rast plot
#' @importFrom grDevices hcl.colors heat.colors



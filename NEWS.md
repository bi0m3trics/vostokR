# vostokR News

## vostokR v0.1.1 (Release date: 2025-09-01)

This update implements comprehensive parallelization and performance optimizations that deliver 2-3x faster processing for solar potential calculations while maintaining full accuracy.

### New features

• **New:** `add_normals()` optimized with 30% performance improvement using `crossprod()` and eigenvalue optimizations
• **New:** Optional geometric features in `add_normals()` with `add_features = TRUE` parameter
• **New:** Geometric features: linearity, planarity, sphericity, curvature for structural analysis
• **New:** Full OpenMP parallelization with thread control functions
• **New:** `set_vostokr_threads()` - Control number of OpenMP threads  
• **New:** `get_vostokr_threads()` - Get current thread count
• **New:** `get_vostokr_performance_info()` - OpenMP status and capabilities
• **New:** `clear_vostokr_caches()` - Clear performance caches
• **New:** Spatial coherence optimization using Morton code Z-order spatial sorting (40%+ cache improvements)
• **New:** SOLPOS result caching for temporal efficiency with intelligent solar position caching
• **New:** Thread-safe shadow caching with thread-local storage and spatial-temporal keys
• **New:** Hierarchical parallelization across multiple levels: days → spatial batches → points
• **New:** Batch processing with optimized memory access patterns using 64-point batches
• **New:** Comprehensive performance testing suite with multi-scale validation
• **New:** Built-in timing and processing rate reporting
• **New:** Cache effectiveness tracking and statistics
• **New:** Performance result export to CSV with comprehensive metrics

### Enhancements

• **Enhanced:** Armadillo matrix operations now use OpenMP (removed ARMA_DONT_USE_OPENMP)
• **Enhanced:** Compiler optimizations (-O3 -march=native) for maximum performance
• **Enhanced:** Memory management and vectorized operations in normal calculations
• **Enhanced:** Automatic coordination with lidR's parallel processing workflow
• **Enhanced:** Documentation throughout the package with updated citations and examples

### Performance improvements

**Normal vector calculation:**
- 30% faster using optimized eigenvalue decomposition
- Better memory management with pre-allocated arrays
- Vectorized operations where possible

**Solar potential calculation:**
- 2-3x performance improvement with OpenMP parallelization
- Up to 40% cache hit improvements with spatial sorting
- Intelligent caching reduces redundant calculations
- Thread-safe processing for large datasets

**Benchmarked results:**

- 100,000 points: 2.30x speedup (114,402 pts/sec)
- 200,000 points: 2.60x speedup (120,929 pts/sec)  
- 400,000 points: 2.64x speedup (117,354 pts/sec)
- 740,240 points: 2.50x speedup (114,269 pts/sec)

### Documentation

• Updated all function documentation with performance improvements
• Added comprehensive examples for geometric features
• Updated citations and author information throughout package
• Enhanced README with simplified usage examples and performance highlights

### Fixes

• Fixed thread control issues preventing proper OpenMP utilization
• Resolved race conditions in multi-threaded processing
• Improved error handling and edge case management
• Better integration with lidR's parallel processing workflow

---

## vostokR v0.1.0 (Release date: 2025-08-29)

Initial release of vostokR - Solar Potential Analysis for LiDAR Point Clouds in R.

### New features

• **New:** `calculate_solar_potential()` - Main function for computing solar irradiance on point clouds
• **New:** `solar_ground_raster()` - Generate raster maps of solar potential for ground surfaces
• **New:** `add_normals()` - Calculate surface normals using k-nearest neighbor approach
• **New:** Modern terra package support for raster operations
• **New:** Automatic CRS detection for latitude, longitude, and timezone extraction
• **New:** Seamless lidR integration with LAS and LAScatalog objects
• **New:** sf package integration for coordinate system transformations
• **New:** Flexible date range support with configurable temporal resolution
• **New:** Octree-based shadow casting with efficient 3D spatial indexing
• **New:** Ray casting implementation for high-performance shadow detection
• **New:** Voxel-based processing with configurable resolution
• **New:** Solar elevation filtering to exclude low-angle sun positions
• **New:** Multi-format point cloud support (.las, .laz, etc.)
• **New:** Classification-aware processing for ground point filtering
• **New:** Memory-efficient processing optimized for large point clouds
• **New:** Comprehensive diagnostic output with processing information

### Documentation

• Complete function documentation with roxygen2
• Comprehensive vignette with tutorials and examples  
• Academic citations for R package and original VOSTOK toolkit
• Working examples using included test data
• Proper authorship attribution and academic references

### Technical implementation

• Native VOSTOK C++ code integration
• RcppArmadillo optimization for linear algebra operations
• Efficient memory management for large-scale processing
• Cross-platform compilation support (Windows, macOS, Linux)
• Automatic timezone detection based on longitude
• CRS transformation support for various coordinate systems
• **New:** **UTM and geographic coordinate support** - Works with projected and unprojected coordinate systems
• **New:** **Spatial extent validation** - Ensures point clouds have appropriate spatial properties for solar analysis

#### Performance Optimization
• **New:** **Octree spatial indexing** - O(log n) complexity for shadow ray queries
• **New:** **Configurable precision** - Balance between accuracy and computation time through voxel size selection
• **New:** **Parallel processing support** - Designed to work with lidR's catalog processing for large datasets  
• **New:** **Memory-efficient algorithms** - Optimized memory usage for processing large point clouds

### Integration with R Spatial Ecosystem

#### lidR Package Integration
• **Method dispatch** - Functions work seamlessly with LAS and LAScatalog objects
• **Attribute preservation** - Maintains existing point cloud attributes while adding solar potential data
• **Processing pipeline** - Integrates naturally with lidR preprocessing and analysis workflows
• **Catalog support** - Designed for processing large datasets using LAScatalog infrastructure

#### Modern Spatial Packages
• **terra package** - All raster operations use modern terra package instead of deprecated raster
• **sf package** - Coordinate system operations leverage sf for robust spatial data handling  
• **data.table efficiency** - Optimized data manipulation using data.table backend from lidR
• **Spatial visualization** - Plotting functions integrate with R's spatial visualization ecosystem

### Example Workflows

The package supports various analysis workflows:

```r
# Basic workflow
las <- readLAS("pointcloud.laz")
las <- add_normals(las, k = 10)
las_solar <- calculate_solar_potential(las,
                                     year = 2025,
                                     start_date = '2025-06-15',
                                     end_date = '2025-06-17')

# Create solar potential raster
ground_raster <- solar_ground_raster(las_solar, res = 2.0)
plot(ground_raster)

# Seasonal analysis
spring_solar <- calculate_solar_potential(las,
                                        start_date = '2025-03-20',
                                        end_date = '2025-06-20',
                                        day_step = 7)
```

### Academic Context

This package implements the VOSTOK solar potential analysis methodology developed by Bechtold and Höfle. The original VOSTOK toolkit provides a comprehensive framework for calculating solar irradiance on 3D point clouds, with particular applications in:

- **Environmental solar potential assessment** - Analyzing solar energy potential in natural areas
- **Urban solar planning** - Evaluating building and infrastructure solar potential  
- **Ecological modeling** - Understanding light availability for habitat modeling
- **Renewable energy assessment** - Supporting solar installation planning and optimization

### Installation and Requirements

**System Requirements:**
- R (>= 4.0.0)
- C++17 compiler
- RcppArmadillo for optimized linear algebra

**Dependencies:**
- lidR (>= 4.0.0) - LiDAR data processing
- terra - Modern raster operations  
- sf (>= 1.0.0) - Spatial data handling
- Rcpp/RcppArmadillo - C++ integration
- data.table - Efficient data manipulation

**Installation:**
```r
# Install from source (development version)
devtools::install("path/to/vostokR")
```

### Citation

When using vostokR in research, please cite both the R package and the original VOSTOK methodology:

```r
citation("vostokR")
```

**R Package Citation:**
Sánchez Meador, A.J. (2025). vostokR: Solar Potential Calculation for Point Clouds using VOSTOK. R package version 0.1.1.

**Original VOSTOK Citation:**  
Bechtold, S., & Höfle, B. (2016). HELIOS: A multi-purpose LiDAR simulation framework for research, planning and training of laser scanning operations with airborne, ground-based mobile and stationary platforms. *ISPRS Journal of Photogrammetry and Remote Sensing*, 115, 86-101.



# NEWS.md

## vostokR v0.1.1 (Release date: 2025-09-01)

**Major Performance Optimization Update - 2-3x Faster Processing**

This update implements comprehensive parallelization and performance optimizations that deliver significant speedup for solar potential calculations while maintaining full accuracy.

### PERFORMANCE IMPROVEMENTS

#### OpenMP Parallelization
• **New:** Full OpenMP parallelization with thread control functions
• **New:** `set_vostokr_threads()` - Control number of OpenMP threads
• **New:** `get_vostokr_threads()` - Get current thread count
• **New:** `get_vostokr_performance_info()` - OpenMP status and capabilities
• **New:** `clear_vostokr_caches()` - Clear performance caches
• **Enhanced:** Armadillo matrix operations now use OpenMP (removed ARMA_DONT_USE_OPENMP)

#### Advanced Optimization Features  
• **New:** **Spatial coherence optimization** - Morton code Z-order spatial sorting for 40%+ cache hit improvements
• **New:** **SOLPOS result caching** - Temporal efficiency with intelligent solar position caching
• **New:** **Thread-safe shadow caching** - Thread-local storage with spatial-temporal keys
• **New:** **Hierarchical parallelization** - Multi-level: days → spatial batches → points
• **New:** **Batch processing** - Optimized memory access patterns with 64-point batches
• **Enhanced:** Compiler optimizations (-O3 -march=native) for maximum performance

#### Performance Testing & Validation
• **New:** Comprehensive multi-scale performance testing suite
• **New:** Large-scale validation up to 740,240 points (full test.laz)
• **New:** Thread efficiency analysis and scalability testing
• **New:** Performance visualization with ggplot2 plots and CSV exports
• **New:** Thread control diagnostic tools for debugging
• **Enhanced:** Automated performance monitoring and result validation

#### Performance Monitoring
• **New:** Built-in timing and processing rate reporting  
• **New:** Cache effectiveness tracking and statistics
• **New:** Thread control verification with `get_actual_vostokr_threads()`
• **New:** Performance result export to CSV with comprehensive metrics
• **Enhanced:** Automatic coordination with lidR's parallel processing workflow

### BENCHMARKED RESULTS

#### Large-Scale Performance (Full Dataset: 740,240 points)
• **2.50x speedup** on full production dataset (740,240 points)
• **114,269 points/second** processing rate with 7 threads
• **Consistent scaling**: 2.30-2.64x speedup across dataset sizes (100K-740K points)
• **Production ready**: Successfully tested up to full test.laz size with excellent performance

#### Multi-Scale Performance Validation
• **100,000 points**: 2.30x speedup (114,402 pts/sec)
• **200,000 points**: 2.60x speedup (120,929 pts/sec)  
• **400,000 points**: 2.64x speedup (117,354 pts/sec)
• **740,240 points**: 2.50x speedup (114,269 pts/sec)

#### Thread Control & Accuracy
• **Fixed thread control**: Explicit n_threads parameter prevents auto-calculation override
• **Race condition handling**: Proper thread-safe implementation with minimal variance
• **<0.1% accuracy deviation** maintained across all threading configurations

### TECHNICAL ENHANCEMENTS
• **New:** Thread-safe octree implementation with enhanced shadow calculations
• **New:** Spatial sorting using Morton codes for optimal cache locality
• **New:** Smart integration auto-coordinates with lidR thread usage
• **New:** Explicit thread parameter control prevents n_threads auto-calculation override
• **New:** Race condition detection and thread-safe result validation
• **Enhanced:** Memory efficiency improvements reduce bandwidth requirements
• **Enhanced:** Algorithm-level optimizations for shadow ray casting

### TESTING & VALIDATION
• **New:** `performance_test.R` - Comprehensive multi-scale performance testing
• **New:** `large_scale_test.R` - Production-scale validation up to 740K points
• **New:** `thread_control_test.R` - Thread control diagnostic testing
• **New:** Performance visualization with scalability plots and efficiency analysis
• **Enhanced:** Automated result validation and consistency checking

## vostokR v0.1.0 (Release date: 2025-08-29)

**Initial release of vostokR - Solar Potential Analysis for LiDAR Point Clouds in R**

vostokR provides R bindings for the VOSTOK (Viewshed Obstruction by Solar Tracking with Octree Knowledge) solar potential analysis toolkit. This package enables researchers and practitioners to calculate solar irradiance potential on 3D point clouds while accounting for shadowing effects from surrounding vegetation and terrain features.

### NEW Features

#### Core Functionality
• **New:** `calculate_solar_potential()` - Main function for computing solar irradiance on LiDAR point clouds using native VOSTOK C++ implementation
• **New:** `solar_ground_raster()` - Generate raster maps of solar potential for ground-level surfaces
• **New:** `add_normals()` - Calculate surface normals required for solar potential analysis using k-nearest neighbor approach
• **New:** Comprehensive **Rcpp integration** with optimized C++ VOSTOK algorithms for high-performance solar calculations

#### Spatial Integration
• **New:** **Modern terra package support** - Full integration with terra for raster operations, replacing deprecated raster package dependencies
• **New:** **Automatic CRS detection** - Intelligently extracts latitude, longitude, and timezone from point cloud coordinate systems for accurate solar position calculations
• **New:** **Seamless lidR integration** - Native support for LAS and LAScatalog objects with S3 method dispatch
• **New:** **sf package integration** - Leverages sf for coordinate system transformations and spatial operations

#### Temporal Analysis
• **New:** **Flexible date range support** - Specify analysis periods using `start_date` and `end_date` parameters with automatic day-of-year conversion
• **New:** **Configurable temporal resolution** - Control sampling frequency with `day_step` and `minute_step` parameters
• **New:** **Solar position algorithms** - Implements SOLPOS algorithms for precise sun position tracking throughout analysis periods

#### Advanced Algorithms
• **New:** **Octree-based shadow casting** - Efficient 3D spatial indexing using octree data structures for rapid shadow detection
• **New:** **Ray casting implementation** - High-performance ray tracing from each point toward sun positions to detect shadowing obstacles
• **New:** **Voxel-based processing** - Configurable voxel resolution (`voxel_size` parameter) for balancing accuracy and computational efficiency
• **New:** **Solar elevation filtering** - `min_sun_angle` parameter to exclude low-angle sun positions that may introduce noise

#### Data Processing
• **New:** **Multi-format point cloud support** - Works with .las, .laz, and other lidR-supported formats
• **New:** **Classification-aware processing** - `solar_ground_raster()` can filter by point classification (default: ground points with Classification == 2)
• **New:** **Fallback processing modes** - `use_all_points` parameter for datasets without proper ground classification
• **New:** **Memory-efficient processing** - Optimized for large point clouds with configurable processing parameters

#### Visualization and Output
• **New:** **Native terra raster output** - All raster functions return SpatRaster objects for modern spatial data workflows  
• **New:** **lidR plotting integration** - Solar potential visualization uses lidR's default plotting methods
• **New:** **Comprehensive diagnostic output** - Functions provide detailed information about processing parameters and results
• **New:** **Units and metadata** - Solar potential values in kWh/m²/year with proper spatial reference information

### Documentation and Examples

#### Comprehensive Documentation  
• **New:** **Complete function documentation** - All functions documented with roxygen2 including parameters, return values, and examples
• **New:** **Comprehensive vignette** - Detailed tutorial covering basic workflows, advanced usage, parameter sensitivity, and integration with other packages
• **New:** **Academic citations** - Proper citation information for both the R package and original VOSTOK toolkit
• **New:** **Code examples** - Working examples using included test data demonstrate all major functionality

#### Package Infrastructure
• **New:** **Test data included** - Package includes `test.laz` file in `inst/extdata/` for examples and testing
• **New:** **Proper authorship attribution** - Andrew Sánchez Meador as primary author, with Sebastian Bechtold and Bernhard Höfle credited for original VOSTOK implementation  
• **New:** **Academic references** - Links to original VOSTOK research and implementation
• **New:** **GPL v3 licensing** - Open source license ensuring broad accessibility

### Technical Implementation

#### C++ Integration
• **New:** **Native VOSTOK C++ code** - Direct integration of original high-performance C++ implementation
• **New:** **RcppArmadillo optimization** - Leverages optimized linear algebra for vector operations
• **New:** **Memory management** - Efficient memory handling for large-scale point cloud processing
• **New:** **Cross-platform compilation** - Support for Windows, macOS, and Linux platforms

#### Coordinate System Handling
• **New:** **Automatic timezone detection** - Longitude-based timezone calculation for accurate solar position modeling
• **New:** **CRS transformation support** - Handles various coordinate reference systems with automatic geographic coordinate extraction
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
las <- readLAS("forest.laz")
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

- **Forest solar potential assessment** - Analyzing solar energy potential in forested areas
- **Urban solar planning** - Evaluating building and infrastructure solar potential  
- **Ecological modeling** - Understanding light availability for vegetation modeling
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
Sánchez Meador, A. (2025). vostokR: Solar Potential Calculation for Point Clouds using VOSTOK. R package version 0.1.0.

**Original VOSTOK Citation:**  
Bechtold, S., & Höfle, B. (2016). HELIOS: A multi-purpose LiDAR simulation framework for research, planning and training of laser scanning operations with airborne, ground-based mobile and stationary platforms. *ISPRS Journal of Photogrammetry and Remote Sensing*, 115, 86-101.

### Development and Contributions

**Primary Author:** Andrew Sánchez Meador (Northern Arizona University)  
**Original VOSTOK Implementation:** Sebastian Bechtold and Bernhard Höfle  
**License:** GPL (>= 3)  
**Development Repository:** Integration with R ecosystem for VOSTOK solar potential toolkit

### Future Development

Planned enhancements for future releases:
- Extended temporal analysis functions
- Additional shadowing algorithms  
- Performance optimizations for very large datasets
- Integration with additional spatial analysis packages
- Enhanced visualization capabilities

---

*For detailed usage instructions and advanced examples, see the package vignette: `vignette("vostokR")`*

*Report bugs and request features at the package development repository.*

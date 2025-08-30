#' Add Normal Vectors to Point Cloud
#'
#' This function adds normal vectors to a LiDAR point cloud using k-nearest neighbors.
#' The normals are computed using the covariance matrix of neighboring points.
#'
#' @param las LAS object from lidR package
#' @param k Number of neighbors to use for normal estimation (default: 10)
#' @return LAS object with added normal vectors (nx, ny, nz)
#' @export
#' @importFrom lidR add_attribute knn readLAS
#' @importFrom stats na.omit
#' @examples
#' \dontrun{
#' library(lidR)
#' library(vostokR)
#' 
#' # Load test data
#' LASfile <- system.file("extdata", "test.laz", package="vostokR")
#' las <- readLAS(LASfile)
#' 
#' # Add normals using 10 nearest neighbors
#' las <- add_normals(las, k = 10)
#' 
#' # Check that normals were added
#' names(las@data)
#' }
#' @seealso \code{\link{calculate_solar_potential}} for solar potential calculation
add_normals <- function(las, k = 10) {
    if (!inherits(las, "LAS")) {
        stop("Input must be a LAS object from lidR package")
    }
    
    # Initialize matrices for normals
    n <- nrow(las@data)
    normals <- matrix(0, n, 3)
    
    # Function to compute normal from covariance matrix
    compute_normal <- function(points) {
        if (nrow(points) < 3) return(c(0, 0, 1))  # Default if not enough points
        
        # Center the points
        centered <- scale(points, center = TRUE, scale = FALSE)
        
        # Compute covariance matrix
        cov_mat <- t(centered) %*% centered / (nrow(points) - 1)
        
        # Compute eigenvectors
        eig <- eigen(cov_mat)
        
        # The normal is the eigenvector corresponding to the smallest eigenvalue
        normal <- eig$vectors[, which.min(eig$values)]
        
        # Ensure normal points "up" (positive Z)
        if (normal[3] < 0) normal <- -normal
        
        return(normal)
    }
    
    # Find k nearest neighbors using lidR's knn function
    nn <- lidR::knn(las, k = k)

    # Get coordinates for all points using lidR convention
    coords <- as.matrix(las@data[, c("X", "Y", "Z")])

    # Compute normals for each point
    for (i in 1:n) {
        # Get neighbor indices from knn result, robust to vector/matrix
        if (is.matrix(nn)) {
            idx <- nn[i, !is.na(nn[i,])]
        } else {
            idx <- nn[!is.na(nn)]
        }
        if (length(idx) >= 3) {
            # Get neighbor coordinates and compute normal
            neighbor_coords <- coords[idx, , drop = FALSE]
            normals[i,] <- compute_normal(neighbor_coords)
        } else {
            normals[i,] <- c(0, 0, 1)  # Default normal if not enough neighbors
        }
    }
    
    # Add normals to LAS object
    las <- add_attribute(las, normals[,1], "nx")
    las <- add_attribute(las, normals[,2], "ny")
    las <- add_attribute(las, normals[,3], "nz")
    
    return(las)
}

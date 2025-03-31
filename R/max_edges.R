#' Maximum number of edges in your grid
#' 
#' The maximum number of edges in your grid, determining what type of transition function you can use in \link{shortest_paths}.
#' 
#' @param rst SpatRaster (terra), RasterLayer (raster), matrix, or list of matrices object. RasterLayers are converted to SpatRasters.
#' @param contiguity \code{"queen"} (default) for queen's case contiguity or \code{"rook"} for rook's case contiguity. In the latter case, the algorithm 
#' only moves between horizontally or vertically adjacent cells. In the former case, it is also travels between diagonally adjacent cells.
#' @param spherical Logical specifying whether coordinates are unprojected, i.e. lonlat, and refer to a sphere, e.g., a planet, if \code{rst} is a matrix 
#' or a list of matrices. It defaults to \code{TRUE}. If \code{FALSE}, the function assumes the coordinates to originate from a planar projection. This 
#' argument has no effect when \code{rst} is a SpatRaster or RasterLayer.
#' @param extent Vector of length 4, specifying the extent of \code{rst}, if \code{rst} is a matrix or a list of matrices. It must contain xmin, xmax, 
#' ymin, and ymax, in that order. The argument has no effect when \code{rst} is a SpatRaster or RasterLayer.
#' 
#' @details An edge is a connection between adjacent grid cells. With queen's case contiguity, each cell has up to 8 edges. With rook's case contiguity, 
#' they have up 4 edges. It is up to 8 and 4, and not exactly 8 or 4, because \code{rst}'s outer pixels do not have neighbors in all directions.
#' 
#' If the data is unprojected and spans from 180 degrees West to 180 degrees East, the easternmost cells are connected to the westernmost cells. Like in 
#' any other geo-spatial software, this is the only scenario in which the algorithm connects cells across the grid boundary, i.e. in which, e.g., a 
#' shortest path leaves the grid on one side and enters it on the opposite side. In all other cases, cells are only connected to their direct neighbors 
#' within the grid.
#' 
#' Another source of edge removal are NA cells. There are no edges to or from NA cells in \code{rst}.
#' 
#' @md
#' @returns Returns a numeric value denoting the maximum number of edges in \code{rst} that \code{shortest_paths} may store in an adjacency list in R at 
#' some point. If there are any NA cells, the returned number is greater than the final graph's number of edges for two reasons. First, 
#' \link{shortest_paths} assembles the adjacency list in multiple steps. Second, for efficiency reasons, \code{max_edges} does not evaluate where in the 
#' grid the NA cells are and assumes the most conservative impact.
#' 
#' The returned value is the same number upon which \link{shortest_paths} decides to either construct the adjacency list in R or in C++. If the result is 
#' larger than `r format(.Machine$integer.max, big.mark = ",")` (the maximum number of elements native R objects can store), it chooses C++. Otherwise, it 
#' selects R. In the C++ case, any \code{tr_fun} transition function must be an Rcpp C++ function with various restrictions (see the 
#' \href{../doc/transition_functions.html}{transition functions vignette}).
#' 
#' @seealso \link{shortest_paths}.
#' 
#' @examples
#' # Generate example data
#' set.seed(2L)
#' input_grid <- terra::rast(crs = "epsg:4326", resolution = 2, vals = sample(c(1L, NA_integer_),
#'   16200L, TRUE, c(0.8, 0.2)))
#' 
#' # Obtain maximum number of edges
#' max_edges(input_grid)
#' 
#' @export
max_edges <- function(rst, contiguity = c("queen", "rook"), spherical = TRUE, extent = NULL) {
  
  # Convert RasterLayer to SpatRaster
  class_rst <- class(rst)
  if("RasterLayer" %chin% class_rst) {
    rst <- terra::rast(rst)
    rst_terra <- TRUE
  } else if("SpatRaster" %chin% class_rst) {
    rst_terra <- TRUE
  } else if(is.matrix(rst)) {
    rst_terra <- FALSE
    rst_list <- FALSE
    n_grids <- 1L
  } else if(is.list(rst)) {
    rst_terra <- FALSE
    n_grids <- length(rst)
    if(n_grids == 0L) stop("rst must not be empty")
    if(!all(vapply(rst, is.matrix, logical(1L), USE.NAMES = FALSE))) stop("Not all rst list elements are matrices")
    rst_list <- TRUE
  } else {
    stop("rst must be a SpatRaster, RasterLayer, matrix, or list of matrices")
  }
  
  # Check CRS
  if(rst_terra) {
    r_crs <- terra::crs(rst)
    if(r_crs == "") stop("rst without specified CRS")
    lonlat <- terra::is.lonlat(rst)
    n_cells <- terra::ncell(rst)
    if(n_cells > .Machine$integer.max) stop("rst exceeds with ", n_cells, " cells this function's current limit of ", .Machine$integer.max, " cells")
    n_cells <- as.integer(n_cells)
    rst_ncol <- as.integer(terra::ncol(rst))
    rst_xmin <- terra::xmin(rst)
    rst_xmax <- terra::xmax(rst)
  } else {
    if(length(spherical) != 1L || !is.logical(spherical) || is.na(spherical)) stop("spherical must be logical and of length one")
    lonlat <- spherical
    if(!is.vector(extent) || length(extent) != 4L || !is.numeric(extent) || anyNA(extent)) stop("extent must be a numeric vector of length 4")
    rst_xmin <- extent[1L]
    rst_xmax <- extent[2L]
    rst_ymin <- extent[3L]
    rst_ymax <- extent[4L]
    if(rst_xmin >= rst_xmax) stop("xmax (extent[2]) must be larger than xmin (extent[1])")
    if(rst_ymin >= rst_ymax) stop("ymax (extent[4]) must be larger than ymin (extent[3])")
    if(lonlat) {
      if(rst_xmin < -180) stop("xmin (extent[1]) must not be smaller than -180, if spherical = TRUE")
      if(rst_xmax > 180) stop("xmax (extent[2]) must not be larger than 180, if spherical = TRUE")
      if(rst_ymin < -90) stop("ymin (extent[3]) must not be smaller than -90, if spherical = TRUE")
      if(rst_ymax > 90) stop("ymax (extent[2]) must not be larger than 90, if spherical = TRUE")
    }
    if(rst_list) {
      n_cells <- length(rst[[1L]])
      rst_nrow <- nrow(rst[[1L]])
      rst_ncol <- ncol(rst[[1L]])
      if(!all(vapply(rst[2:n_grids], function(x) nrow(x) == rst_nrow && ncol(x) == rst_ncol, logical(1L), USE.NAMES = FALSE))) {
        stop("rst matrices differ in dimensions")
      }
    } else {
      n_cells <- length(rst)
      rst_ncol <- ncol(rst)
    }
  }
  global <- lonlat && abs(rst_xmax - rst_xmin - 360) < 0.0001
  
  contiguity <- match.arg(contiguity)
  queen <- contiguity == "queen"
  
  if(rst_terra) {
    n_cells_na <- n_cells - sum(stats::complete.cases(terra::values(rst)))
  } else if(rst_list) {
    n_cells_na <- n_cells - sum(stats::complete.cases(lapply(rst, function(m) as.vector(m))))
  } else {
    n_cells_na <- sum(as.vector(is.na(rst)))
  }
  
  max_neighbors <- get_max_neighbors(n_cells, n_cells_na, queen, rst_ncol, n_cells / rst_ncol, global)
  
  return(max_neighbors)
}



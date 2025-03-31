#' Shortest paths and/ or distances between locations
#' 
#' The shortest paths and/ or distance between locations in a grid according to Dijkstra's (1959) algorithm.
#' 
#' @param rst SpatRaster (terra), RasterLayer (raster), matrix, or list of matrices object. RasterLayers are converted to SpatRasters. Pixels with non-NA 
#' values in all layers mark the cells through which the algorithm may pass. The values of \code{rst} can be accessed in \code{tr_fun}.
#' @param origins Origin points at which the shortest paths start. If \code{rst} is a SpatRaster or RasterLayer object, these points can be passed as a 
#' single SpatVector (terra), sf (sf), or Spatial* (sp) object. sf and sp objects are converted to SpatVectors. Polygons and lines are converted to points 
#' using their centroid. If \code{rst} is a matrix or list of matrices, \code{origins} must be a single matrix, data.frame, or data.table of coordinates 
#' with columns named \code{"x"} and \code{"y"}. The coordinates must refer to points in the reference system that \code{rst} utilizes. Lines and 
#' polygons are thus not accepted in this case. Details on which points the function connects are outlined below.
#' @param destinations Destination points to which the shortest paths are derived. It defaults to \code{NULL}, resulting in the function to compute 
#' shortest paths between the \code{origins} points. Otherwise, the same input rules as for \code{origins} apply. Details on which points the function 
#' connects are outlined below.
#' @param output \code{"distances"} (default), \code{"lines"}, or \code{"both"}. \code{"distances"} lists the total transition costs along the shortest 
#' paths. By default, it is the distance between origin and destination in meters, if \code{rst} is an unprojected SpatRaster or RasterLayer or if 
#' \code{dist_comp = "terra"}. Otherwise, it is denoted in the projection's units. If you pass another function to \code{tr_fun}, the total transition cost 
#' is measured in the units of \code{tr_fun}'s results. \code{"lines"} returns the shortest paths as spatial lines. \code{"both"} returns both distances 
#' and lines. \code{"distances"} is faster and requires less RAM than \code{"lines"} or \code{"both"}.
#' @param output_class Class of the returned object. With \code{output = "distances"}, the options are \code{"data.table"} (default) and 
#' \code{"data.frame"}. With \code{output = "lines"} or \code{output = "both"}, the options are \code{"SpatVector"} (default when \code{rst} is a 
#' SpatRaster or RasterLayer) and \code{"list"} (default when \code{rst} is a matrix or a list of matrices). In the case of \code{"list"}, the attributes, 
#' the line coordinates, and the CRS are returned as individual list elements. The first element in the list of line coordinates refers to the first row in 
#' the attributes table etc. \code{"SpatVector"} is only available with a SpatRaster or RasterLayer \code{rst}. \code{output_class} changes the format of 
#' the returned object, not the information it contains.
#' @param origin_names Character specifying the name of the column in the \code{origins} object used to label the origins in the output object. It 
#' defaults to row numbers.
#' @param destination_names Character specifying the name of the column in the \code{destinations} object used to label the destinations in the output 
#' object. It defaults to row numbers.
#' @param pairwise Logical specifying whether to compute pairwise paths, if \code{origins} and \code{destinations} have equally many rows. If \code{TRUE}, 
#' the function computes the shortest path between the first origin and the first destination, the second origin and the second destination, etc. 
#' Otherwise, it derives the shortest paths from all origins to all destinations. \code{pairwise = TRUE} can alter the order in which results are returned. 
#' Check the output's origins and destinations columns for the respective order.
#' @param contiguity \code{"queen"} (default) for queen's case contiguity or \code{"rook"} for rook's case contiguity. In the latter case, the algorithm 
#' only moves between horizontally or vertically adjacent cells. In the former case, it is also travels between diagonally adjacent cells. \code{"rook"} is 
#' more efficient than \code{"queen"} as it implies fewer edges.
#' @param spherical Logical specifying whether coordinates are unprojected, i.e. lonlat, and refer to a sphere, e.g., a planet, if \code{rst} is a matrix 
#' or a list of matrices. It defaults to \code{TRUE}. If \code{FALSE}, the function assumes the coordinates to originate from a planar projection. This 
#' argument has no effect when \code{rst} is a SpatRaster or RasterLayer.
#' @param radius Radius of the object, e.g. planet, if \code{spherical = TRUE}. This argument has no effect when \code{rst} is a SpatRaster or RasterLayer.
#' @param extent Vector of length 4, specifying the extent of \code{rst}, if \code{rst} is a matrix or a list of matrices. It must contain xmin, xmax, 
#' ymin, and ymax, in that order. The argument has no effect when \code{rst} is a SpatRaster or RasterLayer.
#' @param dist_comp Method to compute distances between adjacent cells. \code{"spaths"} (default) or \code{"terra"}. The default \code{"spaths"} uses 
#' spherical (Haversine) distances in case of lonlat data and planar (Euclidean) distances in case of projected (non-lonlat) data. The functions are 
#' optimized based on the fact that many inter-pixel distances are identical. Modelling the planet as a perfect sphere is in line with e.g. the s2 package, 
#' but is of course an oversimplification. With \code{"terra"}, the function derives distances via \code{terra::distance}. Because this computes all 
#' inter-pixel distances separately, it is slower than the \code{"spaths"} approach. It does take the non-spherical nature of the planet into account 
#' though. With \code{tr_fun}, you can specify a custom distance function that uses neither \code{"spaths"} nor \code{"terra"} distances.
#' @md
#' @param tr_fun The transition function based on which to compute edge weights, i.e. the travel cost between adjacent grid cells. Defaults to the 
#' geographic distance between pixel centroids. Permitted function parameter names are \code{d} (distance between the pixel centroids), \code{x1} (x 
#' coordinate or longitude of the first cell), \code{x2} (x coordinate or longitude of the second cell), \code{y1} (y coordinate or latitude of the first 
#' cell), \code{y2} (y coordinate or latitude of the second cell), \code{v1} (\code{rst} layers' values from the first cell), \code{v2} (\code{rst} layers' 
#' values from the second cell), and \code{nc} (number of CPU cores according to the \code{ncores} argument). If the data is unprojected, i.e. lonlat, or 
#' if \code{dist_comp = "terra"}, \code{d} is measured in meters. Otherwise, it uses the units of the CRS. If \code{rst} has one layer, the values are 
#' passed to \code{v1} and \code{v2} as vectors, otherwise they are passed as a data table where the first column refers to the first layer, the second 
#' column to the second layer etc. Note that data tables are data frames. If \code{rst} produces a graph with more than 
#' `r format(.Machine$integer.max, big.mark = ",")` edges, the adjacency list cannot be stored in R and any \code{tr_fun} needs to be written in C++ with 
#' a few additional requirements. See the details below.
#' @param v_matrix Logical specifying whether to pass values to \code{v1} and \code{v2} in \code{tr_fun} as matrices (\code{TRUE}) instead of data tables 
#' in the multi-layer case and vectors in the single-layer case (\code{FALSE}). It defaults to \code{FALSE}. Setting it to \code{TRUE} might, e.g., be 
#' useful when defining \code{tr_fun} as a C++ Armadillo function.
#' @param tr_directed Logical specifying whether \code{tr_fun} creates a directed graph. In a directed graph, transition costs can be asymmetric. 
#' Traveling from cells A to B may imply a different cost than traveling from B to A. It defaults to \code{TRUE} and only has an effect when 
#' \code{tr_fun} is not \code{NULL}. The default without \code{tr_fun} constructs an undirected graph.
#' @param pre Logical specifying whether to compute the distances between neighboring cells before executing the shortest paths algorithm in C++. 
#' \code{pre} only has an effect, when no \code{tr_fun} is specified and \code{dist_comp = "spaths"}, as the distances are otherwise imported from R. 
#' \code{TRUE} (default) is in the vast majority of cases faster than \code{FALSE}. \code{FALSE} computes distances between neighboring cells while the 
#' shortest paths algorithm traverses the graph. This requires less RAM, but is slower than \code{TRUE}, unless \code{early_stopping = TRUE} and all points 
#' are close to each other. \code{TRUE}'s speed advantage is even larger, when \code{update_rst} is not \code{NULL}.
#' @param early_stopping Logical specifying whether to stop the shortest path algorithm once the target cells are reached. It defaults to \code{TRUE}, 
#' which can be faster than \code{FALSE}, if the points are close to each other compared to the full set of \code{rst} cells. If at least one points 
#' pair is far from each other, \code{FALSE} is the faster setting. \code{FALSE} computes the distance to all cells and then extracts the distance to the 
#' target cells. It, therefore, does not check for each visited cell, whether it is in the set of targets. \code{TRUE} and \code{FALSE} produce the same 
#' result and only differ in terms of computational performance.
#' @param bidirectional Logical specifying whether to produce paths or distances in both directions, if destinations are not specified and no directed 
#' transition function is given. In that case, the distance and the path from point A to point B is the same as the distance and path from point B to point 
#' A. \code{FALSE} (default) only returns distances or paths in one direction. Declaring \code{TRUE} returns distances or paths in both directions. This 
#' parameter's objective is to control the return object's RAM requirement. It only has an effect, if destinations are not specified and no directed 
#' transition function is given.
#' @param update_rst Object updating \code{rst} with moving barriers. It defaults to \code{NULL}, corresponding to \code{rst} not being updated. If 
#' \code{rst} is a SpatRaster or RasterLayer, \code{update_rst} can be a SpatVector (terra), sf (sf), or Spatial* (sp) object, or a list of them. sf and sp 
#' objects are converted to SpatVectors. The function updates \code{rst} by setting any cell intersecting with \code{update_rst} to NA, thereby not 
#' allowing the shortest paths algorithm to pass through that cell. \code{update_rst} only sets non-NA cells to NA, not vice versa. The elements of 
#' \code{update_rst} always update the unmodified \code{rst}. I.e. if \code{update_rst} is a list of two polygons, the shortest paths are derived three 
#' times: once based on the not updated \code{rst}, denoted layer 0 in the output, once based on \code{rst} updated with the first polygon, referred to as 
#' layer 1, and once based on \code{rst} updated with the second polygon, termed layer 2. The second polygon updates the unmodified \code{rst}, not the 
#' \code{rst} updated by the first polygon. If \code{rst} is a matrix or a list of matrices, \code{update_rst} can be a vector of cell numbers, a matrix, 
#' or a list of either. Analogously to the SpatRaster case, these objects mark which cells to set to NA. As in terra, cell numbers start with 1 in the top 
#' left corner and then increase first from left to right and then from top to bottom. The cell numbers in the vector and the NA cells in the matrix 
#' identify the pixels to set to NA. Accordingly, the matrix is of equal dimensions as \code{rst}.
#' @param touches Logical specifying the \code{touches} argument of \code{terra::extract} used when \code{update_rst} is a SpatVector, sf, or Spatial* 
#' object. It defaults to \code{TRUE}. If \code{FALSE}, the function only removes cells on the line render path or with the center point inside a polygon.
#' @param ncores An integer specifying the number of CPU cores to use. It defaults to the number of cores installed on the machine. A value of 1 induces a 
#' single-threaded process.
#' @param par_lvl \code{"points"} or \code{"update_rst"}, indicating the level at which to parallelize when using multiple cores and \code{update_rst} is a 
#' list. \code{"points"} parallelizes over the origin (and destination) point combinations in both the base grid not updated by \code{update_rst} and the 
#' grids updated with \code{update_rst}. The default \code{"update_rst"} is equivalent to \code{"points"} in the base grid, but parallelizes at the 
#' \code{update_rst} list level in the updated grid stage.
#' @param show_progress Logical specifying whether the function prints messages on progress. It defaults to \code{FALSE}.
#' @param bar_limit Integer specifying until up to how many paths or list elements of \code{update_rst} to display a progress bar, if 
#' \code{show_progress = TRUE}. It defaults to 150, in which case the function prints one \code{=} per computed path, if there are no more than 150 paths 
#' requested. In the grids updated with \code{update_rst}, the function displays one \code{=} per processed \code{update_rst} list element, not per path. 
#' In parallel applications, the progress bar can notably slow the execution as the functions only permit one thread to write to output at a time. Do not 
#' set the argument too high to avoid R crashes from text buffer overflows.
#' @param path_type Data type with which C++ stores cell numbers. \code{"int"} (default) is the 4 byte signed integers that R also uses and is the fastest 
#' option. \code{"unsigned short int"} is a 2 byte unsigned integer which requires less RAM than \code{"int"}, but only works if there are less than 65,535 
#' non-NA cells and is comparatively slower because it requires type conversions.
#' @param distance_type Data type with which C++ stores distances. \code{"double"} (default) is a double precision 8 byte floating point number. It is the 
#' fastest and most precise option and also used by R as its \code{numeric} data type. \code{"float"} is a single precision 4 byte floating point number, 
#' which stores decimal values less precisely than \code{"double"} does and is comparatively slower because it requires type conversions. \code{"int"} and 
#' \code{"unsigned short int"} are the integer types described in the \code{path_type} documentation above. With \code{"int"} and 
#' \code{"unsigned short int"}, distances are rounded to integers. When employing these integers types, the distance between any cells in \code{rst}, not 
#' just the cells of interest, must not exceed 2,147,483,647 and 65,535 respectively. The distance difference caused by rounding double values to another 
#' type can accumulate along the shortest paths and can result in notable distance deviation in the output. The recommendations is to stick with the 
#' default \code{"double"} unless the machine does not have enough RAM to run the function otherwise.
#' 
#' @details This function computes shortest paths and/ or distance between locations in a grid using Dijkstra's algorithm. Examples are a ship navigating 
#' around land masses or a hiker traversing mountains.
#' 
#' Let us explore the ship example to illustrate how \code{shortest_paths} works. To compute shortest paths between ports around the world, you start with 
#' a global SpatRaster, in which all land pixels are set to NA and all ocean pixels are set a non-NA value, e.g. 1. A SpatVector marks port locations as 
#' points on water pixels. Passing these two objects to the parameters \code{rst} and \code{origins} respectively derives the shortest paths from each port 
#' to all other ports conditional on ships solely traversing water pixels and returns the distances, i.e. lengths of these paths. If you are not interested 
#' in the distances, but in the spatial lines themselves, set \code{output} to \code{"lines"}. If you want to obtain both, set it to \code{"both"}.
#' 
#' In a different application, you do not want to compute the paths between all ports, but only the paths between the ports on the northern hemisphere and 
#' the ports on the southern hemisphere, but not the paths between ports within the same hemisphere. To assess this, you split the ports into two 
#' SpatVectors and pass them to \code{origins} and \code{destinations} respectively. What if you do not want to connect all orgins to all destinations? Set 
#' \code{pairwise} to \code{TRUE} to connect the first origin just to the first destination, the second origin to the second destination, etc.
#' 
#' By default, the distance or transition cost between adjacent cells of the input grid is the geographic distance between the cells' centroids. What if 
#' the boat is a sailing vessel that minimizes travel time conditional on wind speed, wind direction, and ocean currents? Construct a SpatRaster with 
#' three layers containing information on the three variables respectively. Define a transition function that combines the three layers into a travel time 
#' measure and pass the SpatRaster to \code{rst} and this function to \code{tr_fun}. \code{tr_fun} makes this package very versatile. With custom 
#' transition functions, you can take this software out of the geo-spatial context and, e.g., apply it to biomedical research.
#' 
#' If \code{rst} produces a graph with no more than `r format(.Machine$integer.max, big.mark = ",")` edges, \code{tr_fun} can either be an R or an Rcpp C++ 
#' function returning a numeric R vector. Beyond that limit, the number of edges exceed the maximum of elements that R's native data structures can store. 
#' The data then has to be kept in C++, and most function inputs and output are \code{Rcpp::XPtr} types. \link{max_edges} tells you how many edges your 
#' \code{rst} could produce and, hence, what type of \code{tr_fun} you may use. Check the 
#' \href{../doc/transition_functions.html}{transition functions vignette} for details.
#' 
#' Further boosting efficiency, \code{shortest_paths} allows you to handle multiple tasks in one function call. In the ship routing example, consider that 
#' there are hurricanes in the Caribbean. Ships traveling from India to Australia do not care, but ships traveling from Mexico to the Netherlands have to 
#' go around the storm and must not take the shortest path through the hurricane. You have ten SpatVector polygons delineating the extent of the hurricane 
#' on ten different days. You want to know what the shortest paths are given that ships must go around the polygon on that day. Calling 
#' \code{shortest_paths} ten times with ten different SpatRasters would be very inefficient. This would assemble the graph ten times and recompute also 
#' paths unaffected by the hurricanes, such as the path between India and Australia, in each iteration. Instead, pass the SpatVector polygons to 
#' \code{update_rst}. \code{shortest_paths} then produces the shortest paths for a hurricane-free route and all ten hurricane days, only reestimating the 
#' paths that are affected by the hurricane polygon on a specific day.
#' 
#' Applications to Earth should always pass a SpatRaster or RasterLayer to \code{rst}. The option to use a matrix or a list of matrices is meant for 
#' applications to other planets, non-geo-spatial settings, and users who cannot install the terra package on their system.
#' 
#' The largest source of runtime inefficiency is the quantity of non-NA pixels in the \code{rst} grid. Limit the \code{rst} argument to the relevant area. 
#' E.g., crop the grid to the North Atlantic when computing shipping routes between Canada and France. And set regions through which the shortest path does 
#' certainly not pass to NA.
#' 
#' \code{shortest_paths} is optimized for computational performance. Most of its functions are written in C++ and it does not use a general purpose graph 
#' library, but comes with its custom graph-theoritical implementation tailored to gridded inputs.
#' 
#' The \href{../doc/spaths_introduction.html}{introduction vignette} provides further details on the package.
#' 
#' @returns If `output = "distances"`, the output is by default returned as a data table. If you want the result to be a data frame only, not a data table, 
#' set `output_class` to `"data.frame"`. If `output` is `"lines"` or `"both"`, the the function returns a SpatVector, if `rst` is a SpatRaster or a 
#' RasterLayer, and a list, if `rst` is a matrix or a list of matrices. Explicitly setting `output_class` to `"list"` returns a list in any case. 
#' `output_class = "SpatVector"`, however, returns a SpatVector only if `rst` is a SpatRaster or a RasterLayer.
#' 
#' If \code{output = "distances"} or \code{output = "both"}, the \code{distances} variable marks which points are connected. Unconnected point pairs have 
#' an \code{Inf} distance, if \code{distance_type = "double"} or \code{distance_type = "float"}, and an NA distance, if \code{distance_type = "int"} or 
#' \code{distance_type = "unsigned short int"}. If \code{output = "lines"}, the \code{connected} variable marks which points are connected. Points are 
#' connected, when it is possible to travel between them via non-NA cells in \code{rst}.
#' 
#' @seealso \link{max_edges}, \link{rnd_locations}.
#' 
#' @examples
#' \donttest{
#' # Generate example data
#' set.seed(2L)
#' input_grid <- terra::rast(crs = "epsg:4326", resolution = 2, vals = sample(c(1L, NA_integer_),
#'   16200L, TRUE, c(0.8, 0.2)))
#' origin_pts <- rnd_locations(5L, output_type = "SpatVector")
#' origin_pts$name <- sample(letters, 5)
#' destination_pts <- rnd_locations(5L, output_type = "SpatVector")
#' 
#' # Compute distances
#' shortest_paths(input_grid, origin_pts)
#' shortest_paths(input_grid, origin_pts, bidirectional = TRUE)
#' shortest_paths(input_grid, origin_pts, destination_pts)
#' shortest_paths(input_grid, origin_pts, origin_names = "name")
#' shortest_paths(input_grid, origin_pts, destination_pts, pairwise = TRUE)
#' 
#' # Compute lines
#' shortest_paths(input_grid, origin_pts, output = "lines")
#' 
#' # Compute distances and lines
#' shortest_paths(input_grid, origin_pts, output = "both")
#' 
#' # Use custom transition function
#' input_grid[input_grid == 1L] <- stats::runif(terra::freq(input_grid, value = 1L)$count,
#'   max = 1000)
#' shortest_paths(input_grid, origin_pts, tr_fun = function(d, v1, v2) sqrt(d^2 + abs(v2 - v1)^2),
#'   tr_directed = FALSE)
#' 
#' # Compute distances with grid updating
#' barrier <- terra::vect("POLYGON ((-179 1, 30 1, 30 0, -179 0, -179 1))", crs = "epsg:4326")
#' shortest_paths(input_grid, origin_pts, update_rst = barrier)
#' barriers <- list(barrier, terra::vect("POLYGON ((0 0, 0 89, 1 89, 1 0, 0 0))",
#'   crs = "epsg:4326"))
#' shortest_paths(input_grid, origin_pts, update_rst = barriers)
#' }
#' 
#' @references Dijkstra, E. W. 1959. "A note on two problems in connexion with graphs." \emph{Numerische Mathematik} 1 (1): 269–71.
#' 
#' @importFrom data.table :=
#' @importFrom data.table .N
#' @importFrom data.table .SD
#' @importFrom data.table %chin%
#' 
#' @useDynLib spaths, .registration = TRUE
#' 
#' @export
shortest_paths <- function(rst, origins, destinations = NULL, output = c("distances", "lines", "both"), output_class = NULL, origin_names = NULL,
  destination_names = NULL, pairwise = FALSE, contiguity = c("queen", "rook"), spherical = TRUE, radius = 6371010, extent = NULL, dist_comp = c("spaths",
  "terra"), tr_fun = NULL, v_matrix = FALSE, tr_directed = TRUE, pre = TRUE, early_stopping = TRUE, bidirectional = FALSE, update_rst = NULL,
  touches = TRUE, ncores = NULL, par_lvl = c("update_rst", "points"), show_progress = FALSE, bar_limit = 150L, path_type = c("int", "unsigned short int"),
  distance_type = c("double", "float", "int", "unsigned short int")) {
  
  if(length(show_progress) != 1L || !is.logical(show_progress) || is.na(show_progress)) stop("show_progress must be logical and of length one")
  if(show_progress) message("Checking arguments")
  
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
  rm(class_rst)
  
  # Check CRS
  if(rst_terra) {
    r_crs <- terra::crs(rst)
    if(r_crs == "") stop("rst without specified CRS")
    lonlat <- terra::is.lonlat(rst)
    n_cells <- terra::ncell(rst)
    if(n_cells > .Machine$integer.max) stop("rst exceeds with ", n_cells, " cells this function's current limit of ", .Machine$integer.max, " cells")
    n_cells <- as.integer(n_cells)
    rst_nrow <- as.integer(terra::nrow(rst))
    rst_ncol <- as.integer(terra::ncol(rst))
    rst_xmin <- terra::xmin(rst)
    rst_xmax <- terra::xmax(rst)
    rst_ymin <- NULL
    rst_ymax <- terra::ymax(rst)
    rst_xres <- terra::xres(rst)
    rst_yres <- terra::yres(rst)
    rst_list <- FALSE
    n_grids <- terra::nlyr(rst)
    radius2 <- 12742020
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
      if(length(radius) != 1L || !is.numeric(radius) || is.na(radius) || is.infinite(radius) || radius <= 0) {
        stop("radius must be a positive, finite numeric")
      }
      radius2 <- radius * 2
      if(rst_xmin < -180) stop("xmin (extent[1]) must not be smaller than -180, if spherical = TRUE")
      if(rst_xmax > 180) stop("xmax (extent[2]) must not be larger than 180, if spherical = TRUE")
      if(rst_ymin < -90) stop("ymin (extent[3]) must not be smaller than -90, if spherical = TRUE")
      if(rst_ymax > 90) stop("ymax (extent[2]) must not be larger than 90, if spherical = TRUE")
    } else {
      radius2 <- 12742020
    }
    if(rst_list) {
      n_cells <- length(rst[[1L]])
      rst_nrow <- nrow(rst[[1L]])
      rst_ncol <- ncol(rst[[1L]])
      if(!all(vapply(rst[2:n_grids], function(x) nrow(x) == rst_nrow && ncol(x) == rst_ncol, logical(1L), USE.NAMES = FALSE))) {
        stop("rst matrices differ in dimensions")
      }
      rst_names <- names(rst)
      rst <- lapply(rst, function(m) as.vector(t(m)))
      data.table::setDT(rst)
      if(!is.null(rst_names)) data.table::setnames(rst, rst_names)
      rm(rst_names)
    } else {
      n_cells <- length(rst)
      rst_nrow <- nrow(rst)
      rst_ncol <- ncol(rst)
      rst <- as.vector(t(rst))
    }
    r_crs <- NULL
    rst_xres <- (rst_xmax - rst_xmin) / rst_ncol
    rst_yres <- (rst_ymax - rst_ymin) / rst_nrow
  }
  global <- lonlat && abs(rst_xmax - rst_xmin - 360) < 0.0001
  rm(extent, radius)
  
  # Check origin_names and destination_names
  origin_nms_specified <- !is.null(origin_names)
  if(origin_nms_specified && (!(is.character(origin_names) && length(origin_names) == 1L) || is.na(origin_names))) {
    stop("origin_names must be NULL or a character object of length one")
  }
  destination_nms_specified <- !is.null(destination_names)
  if(destination_nms_specified && (!(is.character(destination_names) && length(destination_names) == 1L) || is.na(destination_names))) {
    stop("destination_names must be NULL or a character object of length one")
  }
  
  if(length(pairwise) != 1L || !is.logical(pairwise) || is.na(pairwise)) stop("pairwise must be logical and of length one")
  
  # Convert origins
  origins <- convert_points(origins, rst, r_crs, origin_names, origin_nms_specified, pairwise, rst_terra, rst_list, n_grids, rst_xmin, rst_xmax, rst_ymin,
    rst_ymax, rst_xres, rst_yres, rst_ncol)
  
  # Convert destinations
  dest_specified <- !is.null(destinations)
  if(dest_specified) {
    destinations <- convert_points(destinations, rst, r_crs, destination_names, destination_nms_specified, pairwise, rst_terra, rst_list, n_grids, rst_xmin,
      rst_xmax, rst_ymin, rst_ymax, rst_xres, rst_yres, rst_ncol, FALSE)
  }
  
  rst_xmin <- rst_xmin + 0.5 * rst_xres
  rst_ymax <- rst_ymax - 0.5 * rst_yres
  
  # Check other arguments
  if(pairwise) {
    if(dest_specified) {
      if(length(origins[["cls"]]) != length(destinations[["cls"]])) {
        stop("origins and destinations must have the same number of rows, if pairwise is TRUE")
      }
    } else {
      stop("destinations must not be NULL, if pairwise is TRUE")
    }
  } else if(!dest_specified && length(origins[["cls"]]) < 2L) {
    stop("origins must have more than one row, if destinations are not specified")
  }
  
  # Output type
  output <- match.arg(output)
  return_dists <- output == "both"
  output_lines <- return_dists || (output == "lines")
  if(is.null(output_class)) {
    if(output_lines) {
      output_spatvector <- rst_terra
    } else {
      output_datatable <- TRUE
    }
  } else {
    if(length(output_class) != 1L) stop('output_class must be NULL, "SpatVector", "data.table", "list", or "data.frame"')
    if(output_lines) {
      if(!(output_class %chin% c("SpatVector", "list"))) stop('When output = "', output, '", output_class must be "SpatVector" or "list"')
      output_spatvector <- rst_terra && output_class == "SpatVector"
    } else {
      if(!(output_class %chin% c("data.table", "data.frame"))) stop('When output = "', output, '", output_class must be "data.table" or "data.frame"')
      output_datatable <- output_class == "data.table"
    }
  }
  
  contiguity <- match.arg(contiguity)
  queen <- contiguity == "queen"
  
  dist_comp_terra <- match.arg(dist_comp) == "terra"
  if(dist_comp_terra && !rst_terra) stop('dist_comp = "terra" requires rst to be a SpatRaster or RasterLayer')
  
  if(length(pre) != 1L || !is.logical(pre) || is.na(pre)) stop("pre must be logical and of length one")
  if(length(early_stopping) != 1L || !is.logical(early_stopping) || is.na(early_stopping)) stop("early_stopping must be logical and of length one")
  if(length(bidirectional) != 1L || !is.logical(bidirectional) || is.na(bidirectional)) stop("bidirectional must be logical and of length one")
  
  # Transition function parameters
  tr_fun_specified <- !is.null(tr_fun)
  if(tr_fun_specified) {
    if(is.function(tr_fun)) {
      tr_fun_v <- names(formals(tr_fun))                                        # Extract function parameter names
      if(!all(tr_fun_v %chin% c("d", "x1", "x2", "y1", "y2", "v1", "v2", "nc"))) {
        stop("The permitted parameter names in tr_fun are d, x1, x2, y1, y2, v1, v2, and nc")
      }
      args_used <- c("d", "x1", "y1", "x2", "y2", "v1", "v2", "nc") %chin% tr_fun_v
    } else {
      stop("tr_fun must be either NULL or a function")
    }
    if(length(v_matrix) != 1L || !is.logical(v_matrix) || is.na(v_matrix)) stop("v_matrix must be logical and of length one")
    if(length(tr_directed) != 1L || !is.logical(tr_directed) || is.na(tr_directed)) stop("tr_directed must be logical and of length one")
  }
  if(!dest_specified && origin_nms_specified && !(tr_fun_specified && tr_directed) && !is.character(origins[["nms"]]) && !is.numeric(origins[["nms"]])) {
    stop("origin_names must refer to a character, integer, or numeric variable in this setting")
  }
  wweights <- dist_comp_terra || tr_fun_specified
  
  # Grid updating
  upd_rst_specified <- !is.null(update_rst)
  if(upd_rst_specified) {
    upd_rst_list <- any(class(update_rst) == "list")
    if(upd_rst_list) {
      n_layers <- length(update_rst) + 1L
      if(n_layers == 1L) stop("The list passed to update_rst is empty")
      if(rst_terra) {
        update_rst <- lapply(update_rst, function(v) {
          if(all(class(v) != "SpatVector")) v <- terra::vect(v)
          if(terra::crs(v) != r_crs) v <- terra::project(v, r_crs)
          return(v)
        })
      } else {
        upd_rst_matrix <- is.matrix(update_rst[[1L]])
        if(upd_rst_matrix) {
          if(!all(vapply(update_rst, function(m) is.matrix(m) && nrow(m) == rst_nrow && ncol(m) == rst_ncol, logical(1L), USE.NAMES = FALSE))) {
            stop("If update_rst is a list of matrices, all of these matrices must be of the same dimensions as rst")
          }
        } else {
          update_rst <- lapply(update_rst, function(v) {
            if(!(is.vector(v) && is.numeric(v))) {
              stop("If rst is a matrix or list of matrices, update_rst must be NULL, a matrix, a numeric vector, a list of matrices, or a list of numeric ",
                "vectors")
            }
            if(!is.integer(v)) v <- as.integer(v)
            if(min(v, na.rm = TRUE) < 1L) stop("The values in update_rst must be >= 1")
            if(max(v, na.rm = TRUE) > n_cells) stop("The values in update_rst refer to cells beyond rst")
            return(v)
          })
        }
      }
    } else {
      if(rst_terra) {
        if(all(class(update_rst) != "SpatVector")) update_rst <- terra::vect(update_rst)
        if(terra::crs(update_rst) != r_crs) update_rst <- terra::project(update_rst, r_crs)
      } else {
        upd_rst_matrix <- is.matrix(update_rst)
        if(upd_rst_matrix) {
          if(!(is.matrix(update_rst) && nrow(update_rst) == rst_nrow && ncol(update_rst) == rst_ncol)) {
            stop("If update_rst is a list of matrices, all of these matrices must be of the same dimensions as rst")
          }
        } else {
          if(!(is.vector(update_rst) && is.numeric(update_rst))) {
            stop("If rst is a matrix or list of matrices, update_rst must be NULL, a matrix, a numeric vector, a list of matrices, or a list of numeric ",
              "vectors")
          }
          if(!is.integer(update_rst)) update_rst <- as.integer(update_rst)
          if(min(update_rst, na.rm = TRUE) < 1L) stop("The values in update_rst must be >= 1")
          if(max(update_rst, na.rm = TRUE) > n_cells) stop("The values in update_rst refer to cells beyond rst")
        }
      }
      n_layers <- 2L
    }
    if(length(touches) != 1L || !is.logical(touches) || is.na(touches)) stop("touches must be logical and of length one")
  } else {
    n_layers <- 1L
  }
  
  # Parallel setting
  if(is.null(ncores)) {
    ncores <- parallel::detectCores()
  } else if(length(ncores) != 1L || ncores < 1) {
    stop("ncores must be either NULL or a positive integer")
  }
  ncoresg1 <- ncores > 1
  
  # Distance type
  distance_type <- match.arg(distance_type)
  numeric_weights <- distance_type %chin% c("double", "float")
  if(numeric_weights) {
    double_weights <- distance_type == "double"
    signed_weights <- FALSE
  } else {
    signed_weights <- distance_type == "int"
    double_weights <- FALSE
  }
  
  # Progress bar
  if(length(bar_limit) != 1L || !is.numeric(bar_limit) || is.na(bar_limit) || bar_limit > .Machine$integer.max || bar_limit < 0) {
    stop("bar_limit must be an integer between 0 and ", .Machine$integer.max)
  }
  if(!is.integer(bar_limit)) bar_limit <- as.integer(bar_limit + 0.5)
  
  if(show_progress) message("Converting spatial inputs")
  
  # List non-NA cells as data.table
  if(tr_fun_specified && any(args_used[6:7])) {
    if((rst_terra || (rst_list && !is.null(names(rst)))) && any(c("x", "y", "cell_numbers") %chin% names(rst))) {
      p_names <- c("x", "y", "cell_numbers") %chin% names(rst)
      if(sum(p_names) > 1L) {
        if(sum(p_names) > 2L) {
          p_names <- "x, y, and cell_numbers"
        } else {
          p_names <- paste0(c("x", "y", "cell_numbers")[p_names], collapse = " and ")
        }
        p_names <- paste0(p_names, " are not permitted rst layer names")
      } else {
        p_names <- paste0(c("x", "y", "cell_numbers")[p_names], " is not a permitted rst layer name")
      }
      stop(paste0(p_names, ", when using a transition function with parameter v"))
    }
    if(rst_terra) {
      crd <- terra::values(rst, dataframe = TRUE)                              # Extract layer values
      data.table::setDT(crd)
    } else if(rst_list) {
      crd <- rst
    } else {
      crd <- data.table::data.table(v = rst)
    }
    crd[, cell_numbers := 1:.N]
    crd <- stats::na.omit(crd, cols = 1:(NCOL(crd) - 1L))                      # Subset to cells with non-NA values in all layers
  } else {
    if(rst_terra) {
      if(n_grids > 1) {
        crd <- data.table::data.table(cell_numbers = (1:n_cells)[stats::complete.cases(terra::values(rst))]) # Obtain ids of non-NA cells
      } else {
        crd <- data.table::data.table(cell_numbers = (1:n_cells)[!is.na(terra::values(rst, FALSE))]) # Obtain ids of non-NA cells
      }
    } else {
      if(rst_list) {
        rst <- stats::complete.cases(rst)
      } else {
        rst <- !is.na(rst)
      }
      crd <- data.table::data.table(cell_numbers = (1:n_cells)[rst])
    }
  }
  rm(n_grids)
  if(!(rst_terra && upd_rst_specified)) rm(rst)
  data.table::setkey(crd, cell_numbers)
  
  # set origin (and destinations) to shifted cell values
  origins[["cls"]] <- crd[.(origins[["cls"]]), nomatch = NULL, which = TRUE, on = "cell_numbers"] - 1L
  if(dest_specified) destinations[["cls"]] <- crd[.(destinations[["cls"]]), nomatch = NULL, which = TRUE, on = "cell_numbers"] - 1L
  
  if(show_progress) message("Preparing algorithm inputs")
  
  int_path <- match.arg(path_type) == "int" || nrow(crd) > 65535L
  
  # Obtain adjacency matrix (faster than terra::adjacent)
  # Construct adjacency list in R, if it fits into a data table; otherwise do it in C++
  n_cells_na <- n_cells - nrow(crd)
  max_neighbors <- get_max_neighbors(n_cells, n_cells_na, queen, rst_ncol, rst_nrow, global)
  rm(n_cells_na)
  from_to_r <- max_neighbors <= .Machine$integer.max
  if(from_to_r) {
    if(queen) {
      from_to <- 1:8
    } else {
      from_to <- c(2L, 4L, 5L, 7L)
    }
    from_to <- data.table::rbindlist(lapply(from_to, function(n) {
      if(n == 1L) {
        if(global) {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] - rst_ncol - 1L)[(to %% rst_ncol) == 0L,
            to := to + rst_ncol][to > 0L,])
        } else {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] - rst_ncol - 1L)[to > 0L & (to %% rst_ncol) != 0L,])
        }
      }
      if(n == 2L) {
        return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] - rst_ncol)[to > 0L,])
      }
      if(n == 3L) {
        if(global) {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] - rst_ncol + 1L)[(to %% rst_ncol) == 1L,
            to := to - rst_ncol][to > 0L,])
        } else {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] - rst_ncol + 1L)[to > 0L & (to %% rst_ncol) != 1L,])
        }
      }
      if(n == 4L) {
        if(global) {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] - 1L)[(to %% rst_ncol) == 0L, to := to + rst_ncol])
        } else {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] - 1L)[(to %% rst_ncol) != 0L,])
        }
      }
      if(n == 5L) {
        if(global) {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] + 1L)[(to %% rst_ncol) == 1L, to := to - rst_ncol])
        } else {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] + 1L)[(to %% rst_ncol) != 1L,])
        }
      }
      if(n == 6L) {
        if(global) {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] + rst_ncol - 1L)[(to %% rst_ncol) == 0L,
            to := to + rst_ncol][to <= n_cells,])
        } else {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] + rst_ncol - 1L)[to <= n_cells & (to %% rst_ncol) != 0L,])
        }
      }
      if(n == 7L) {
        return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] + rst_ncol)[to <= n_cells,])
      }
      if(n == 8L) {
        if(global) {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] + rst_ncol + 1L)[(to %% rst_ncol) == 1L,
            to := to - rst_ncol][to <= n_cells,])
        } else {
          return(data.table::data.table(from = crd[["cell_numbers"]], to = crd[["cell_numbers"]] + rst_ncol + 1L)[to <= n_cells & (to %% rst_ncol) != 1L,])
        }
      }
    }), use.names = FALSE)
  } else {
    if(dist_comp_terra) stop('rst is too large for dist_comp = "terra"')
    # C++ function includes below step of excluding NA cells and shifting cell numbers
    if(int_path) {
      from_to <- list(from_to = from_to_ptr_i(crd[["cell_numbers"]], queen, global, rst_ncol, n_cells, tr_fun_specified, max_neighbors))
    } else {
      from_to <- list(from_to = from_to_ptr_u(crd[["cell_numbers"]], queen, global, rst_ncol, n_cells, tr_fun_specified, max_neighbors))
    }
  }
  
  # Shift cell numbers
  if(from_to_r) {
    if(tr_fun_specified || dist_comp_terra) {
      crd[, cell_numbers_shifted := 1:.N]
    } else {
      crd[, cell_numbers_shifted := 0:(.N - 1L)]
    }
    from_to <- from_to[crd[, c("cell_numbers", "cell_numbers_shifted")], nomatch = NULL, on = "to==cell_numbers"][, to := NULL] # Subset to non-NA dest cls
    data.table::setnames(from_to, "cell_numbers_shifted", "to")
    from_to <- as.list(data.table::setnames(from_to[crd[, c("cell_numbers", "cell_numbers_shifted")], nomatch = NULL, on = "from==cell_numbers"][,
      from := NULL], "cell_numbers_shifted", "from"))
    crd[, cell_numbers_shifted := NULL]
  }
  
  # Create coords
  coords_cnumbers_included <- FALSE
  if(output_lines || !wweights) {
    coords <- list(xres = rst_xres, yres = rst_yres, ncol = rst_ncol)
    coords[["ymax"]] <- rst_ymax
    if(pre) {
      coords[["nrow"]] <- rst_nrow
    }
    if(output_lines) {
      coords[["xmin"]] <- rst_xmin
    }
    if(output_lines) {
      coords[["double_coords"]] <- (rst_xres %% 1 != 0) || (rst_yres %% 1 != 0) || (rst_xmin %% 1 != 0) || (rst_ymax %% 1 != 0)
    }
  }
  
  # Compute edge weights based on custom transition function or terra::distance
  # from_to uses one-based indices in R and zero-based indices in c++
  if(tr_fun_specified) {
    if(!from_to_r) {
      if(!grepl("[.]Call[(]", deparse(body(tr_fun)))) stop("Because of the input grid's size, the transition function must be implemented via Rcpp (C++)")
    }
    # Collect transition function arguments
    tr_fun_args <- list()
    if(args_used[1L]) {
      crd_cell_numbers <- crd[["cell_numbers"]] - 1L
      if(dist_comp_terra) { # only relevant when from_to is not an XPtr
        crd_xy <- matrix(c(rst_xmin + (crd_cell_numbers %% rst_ncol) * rst_xres, rst_ymax - as.integer(crd_cell_numbers / rst_ncol) * rst_yres), ncol = 2L)
        tr_fun_args$d <- terra::distance(crd_xy[from_to[["from"]],], crd_xy[from_to[["to"]],], lonlat = lonlat, pairwise = TRUE)
        if(args_used[2L]) tr_fun_args$x1 <- crd_xy[from_to[["from"]], 1L]
        if(args_used[3L]) tr_fun_args$y1 <- crd_xy[from_to[["from"]], 2L]
        if(args_used[4L]) tr_fun_args$x2 <- crd_xy[from_to[["to"]], 1L]
        if(args_used[5L]) tr_fun_args$y2 <- crd_xy[from_to[["to"]], 2L]
        rm(crd_xy)
      } else {
        row_number <- as.integer(crd_cell_numbers / rst_ncol) # zero-indexed row numbers
        if(queen || args_used[2L] || args_used[4L]) col_number <- crd_cell_numbers - row_number * rst_ncol # zero-indexed column numbers
        # Haversine distance
        if(lonlat) {
          yres2 <- sin(rst_yres * pi / 360)
          xres2 <- sin(rst_xres * pi / 360)
          yres3 <- yres2 ^ 2
          d_vertical <- radius2 * atan2(yres2, sqrt(1 - yres3))
          yres4 <- cos((rst_ymax - rst_yres * 0:(rst_nrow - 1L)) * pi / 180)
          d <- yres4 * xres2;
          d_horizontal <- radius2 * atan2(d, sqrt(1 - d ^ 2))
          if(queen) {
            d <- yres3 + yres4 * cos((rst_ymax - rst_yres * 1:rst_nrow) * pi / 180) * xres2 ^ 2
            d_diagonal <- radius2 * atan2(sqrt(d), sqrt(1 - d))
            if(from_to_r) {
              tr_fun_args$d <- data.table::fcase(row_number[from_to[["from"]]] == row_number[from_to[["to"]]],
                d_horizontal[row_number[from_to[["from"]]] + 1L],
                col_number[from_to[["from"]]] == col_number[from_to[["to"]]], d_vertical,
                rep.int(TRUE, length(from_to[["from"]])), d_diagonal[row_number[from_to[["from"]]] + 1L])
            } else {
              if(int_path) {
                tr_fun_args$d <- tr_fun_args_d_haversine_queen_i(from_to[["from_to"]], row_number, col_number, d_horizontal, d_vertical, d_diagonal, ncores)
              } else {
                tr_fun_args$d <- tr_fun_args_d_haversine_queen_u(from_to[["from_to"]], row_number, col_number, d_horizontal, d_vertical, d_diagonal, ncores)
              }
            }
            rm(d_diagonal)
            if(!(args_used[2L] || args_used[4L])) rm(col_number)
          } else {
            if(from_to_r) {
              tr_fun_args$d <- data.table::fifelse(row_number[from_to[["from"]]] == row_number[from_to[["to"]]],
                d_horizontal[row_number[from_to[["from"]]] + 1L], d_vertical)
            } else {
              if(int_path) {
                tr_fun_args$d <- tr_fun_args_d_haversine_rook_i(from_to[["from_to"]], row_number, d_horizontal, d_vertical, ncores)
              } else {
                tr_fun_args$d <- tr_fun_args_d_haversine_rook_u(from_to[["from_to"]], row_number, d_horizontal, d_vertical, ncores)
              }
            }
          }
          rm(yres2, xres2, yres3, d_vertical, yres4, d, d_horizontal)
        # Euclidean distance
        } else {
          if(queen) {
            d_diagonal <- sqrt(xres ^ 2 + yres ^ 2)
            if(from_to_r) {
              tr_fun_args$d <- data.table::fcase(row_number[from_to[["from"]]] == row_number[from_to[["to"]]], rst_xres,
                col_number[from_to[["from"]]] == col_number[from_to[["to"]]], rst_yres,
                default = d_diagonal)
            } else {
              if(int_path) {
                tr_fun_args$d <- tr_fun_args_d_euclidean_queen_i(from_to[["from_to"]], row_number, col_number, rst_xres, rst_yres, d_vertical, d_diagonal,
                  ncores)
              } else {
                tr_fun_args$d <- tr_fun_args_d_euclidean_queen_u(from_to[["from_to"]], row_number, col_number, rst_xres, rst_yres, d_vertical, d_diagonal,
                  ncores)
              }
            }
            rm(d_diagonal)
            if(!(args_used[2L] || args_used[4L])) rm(col_number)
          } else {
            if(from_to_r) {
              tr_fun_args$d <- data.table::fifelse(row_number[from_to[["from"]]] == row_number[from_to[["to"]]], rst_xres, rst_yres)
            } else {
              if(int_path) {
                tr_fun_args$d <- tr_fun_args_d_euclidean_rook_i(from_to[["from_to"]], row_number, rst_xres, rst_yres, ncores)
              } else {
                tr_fun_args$d <- tr_fun_args_d_euclidean_rook_u(from_to[["from_to"]], row_number, rst_xres, rst_yres, ncores)
              }
            }
          }
        }
        if(args_used[2L] || args_used[4L]) {
          col_number <- rst_xmin + col_number * rst_xres
          if(from_to_r) {
            if(args_used[2L]) tr_fun_args$x1 <- col_number[from_to[["from"]]]
            if(args_used[4L]) tr_fun_args$x2 <- col_number[from_to[["to"]]]
          } else {
            if(int_path) {
              if(args_used[2L]) tr_fun_args$x1 <- tr_fun_args_coords_i(from_to[["from_to"]], col_number, TRUE)
              if(args_used[4L]) tr_fun_args$x2 <- tr_fun_args_coords_i(from_to[["from_to"]], col_number, FALSE)
            } else {
              if(args_used[2L]) tr_fun_args$x1 <- tr_fun_args_coords_u(from_to[["from_to"]], col_number, TRUE)
              if(args_used[4L]) tr_fun_args$x2 <- tr_fun_args_coords_u(from_to[["from_to"]], col_number, FALSE)
            }
          }
          rm(col_number)
        }
        if(args_used[3L] || args_used[5L]) {
          row_number <- rst_ymax - row_number * rst_yres
          if(from_to_r) {
            if(args_used[3L]) tr_fun_args$y1 <- row_number[from_to[["from"]]]
            if(args_used[5L]) tr_fun_args$y2 <- row_number[from_to[["to"]]]
          } else {
            if(int_path) {
              if(args_used[3L]) tr_fun_args$y1 <- tr_fun_args_coords_i(from_to[["from_to"]], row_number, TRUE)
              if(args_used[5L]) tr_fun_args$y2 <- tr_fun_args_coords_i(from_to[["from_to"]], row_number, FALSE)
            } else {
              if(args_used[3L]) tr_fun_args$y1 <- tr_fun_args_coords_u(from_to[["from_to"]], row_number, TRUE)
              if(args_used[5L]) tr_fun_args$y2 <- tr_fun_args_coords_u(from_to[["from_to"]], row_number, FALSE)
            }
          }
        }
        rm(row_number)
      }
      if(output_lines) {
        coords_cnumbers_included <- TRUE
        coords[["cell_numbers"]] <- crd_cell_numbers
      }
      rm(crd_cell_numbers)
    } else {
      if(args_used[2L] || args_used[4L]) {
        x <- rst_xmin + ((crd[["cell_numbers"]] - 1L) %% rst_ncol) * rst_xres
        if(from_to_r) {
          if(args_used[2L]) tr_fun_args$x1 <- x[from_to[["from"]]]
          if(args_used[4L]) tr_fun_args$x2 <- x[from_to[["to"]]]
        } else {
          if(int_path) {
            if(args_used[2L]) tr_fun_args$x1 <- tr_fun_args_coords_i(from_to[["from_to"]], x, TRUE)
            if(args_used[4L]) tr_fun_args$x2 <- tr_fun_args_coords_i(from_to[["from_to"]], x, FALSE)
          } else {
            if(args_used[2L]) tr_fun_args$x1 <- tr_fun_args_coords_u(from_to[["from_to"]], x, TRUE)
            if(args_used[4L]) tr_fun_args$x2 <- tr_fun_args_coords_u(from_to[["from_to"]], x, FALSE)
          }
        }
        rm(x)
      }
      if(args_used[3L] || args_used[5L]) {
        y <- rst_ymax - as.integer((crd[["cell_numbers"]] - 1L) / rst_ncol) * rst_yres
        if(from_to_r) {
          if(args_used[3L]) tr_fun_args$y1 <- y[from_to[["from"]]]
          if(args_used[5L]) tr_fun_args$y2 <- y[from_to[["to"]]]
        } else {
          if(int_path) {
            if(args_used[3L]) tr_fun_args$y1 <- tr_fun_args_coords_i(from_to[["from_to"]], y, TRUE)
            if(args_used[5L]) tr_fun_args$y2 <- tr_fun_args_coords_i(from_to[["from_to"]], y, FALSE)
          } else {
            if(args_used[3L]) tr_fun_args$y1 <- tr_fun_args_coords_u(from_to[["from_to"]], y, TRUE)
            if(args_used[5L]) tr_fun_args$y2 <- tr_fun_args_coords_u(from_to[["from_to"]], y, FALSE)
          }
        }
        rm(y)
      }
    }
    if(any(args_used[6:7])) {
      if(v_matrix && !from_to_r) {
        warning("v_matrix is ignored because the graph has more than ", .Machine$integer.max, " edges")
        v_matrix <- FALSE
      }
      if(!v_matrix) {
        v_vars <- setdiff(names(crd), "cell_numbers")
        if(!from_to_r) {
          v_vars_classes <- vapply(v_vars, function(v) class(crd[[v]]), character(1L), USE.NAMES = FALSE)
          layer_class_check(v_vars_classes)
        }
      }
    }
    if(args_used[6L]) {
      if(v_matrix) {
        tr_fun_args$v1 <- as.matrix(crd[from_to[["from"]], .SD, .SDcols = setdiff(names(crd), "cell_numbers")])
      } else {
        if(from_to_r) {
          if(length(v_vars) > 1L) {
            tr_fun_args$v1 <- crd[from_to[["from"]], .SD, .SDcols = v_vars]
          } else {
            tr_fun_args$v1 <- crd[from_to[["from"]], .SD, .SDcols = v_vars][[1L]]
          }
        } else {
          if(int_path) {
            tr_fun_args$v1 <- tr_fun_args_v_i(from_to[["from_to"]], crd, v_vars, v_vars_classes, TRUE)
          } else {
            tr_fun_args$v1 <- tr_fun_args_v_u(from_to[["from_to"]], crd, v_vars, v_vars_classes, TRUE)
          }
        }
      }
    }
    if(args_used[7L]) {
      if(v_matrix) {
        tr_fun_args$v2 <- as.matrix(crd[from_to[["to"]], .SD, .SDcols = setdiff(names(crd), "cell_numbers")])
      } else {
        if(from_to_r) {
          if(length(v_vars) > 1L) {
            tr_fun_args$v2 <- crd[from_to[["to"]], .SD, .SDcols = v_vars]
          } else {
            tr_fun_args$v2 <- crd[from_to[["to"]], .SD, .SDcols = v_vars][[1L]]
          }
        } else {
          if(int_path) {
            tr_fun_args$v2 <- tr_fun_args_v_i(from_to[["from_to"]], crd, v_vars, v_vars_classes, FALSE)
          } else {
            tr_fun_args$v2 <- tr_fun_args_v_u(from_to[["from_to"]], crd, v_vars, v_vars_classes, FALSE)
          }
        }
      }
    }
    if(any(args_used[6:7])) {
      crd[, setdiff(names(crd), "cell_numbers") := NULL]
      if(!v_matrix) {
        rm(v_vars)
        if(!from_to_r) rm(v_vars_classes)
      }
    }
    if(args_used[8L]) tr_fun_args$nc <- ncores
    from_to$weights <- do.call(tr_fun, tr_fun_args[tr_fun_v]) # Run transition function
    if(from_to_r) {
      if(!is.vector(from_to[["weights"]])) {
        stop("tr_fun must return a vector")
      }
      if(!is.numeric(from_to[["weights"]])) {
        stop("tr_fun must return a numeric or integer vector")
      }
      if(any(!is.finite(from_to[["weights"]]))) {
        stop("tr_fun must exclusively return finite values")
      }
      if(min(from_to[["weights"]]) < 0) {
        stop("tr_fun must not return negative values")
      }
      if(length(from_to[["weights"]]) != length(from_to[["from"]])) {
        stop("The number of values returned by tr_fun must equal the number of edges")
      }
      if(is.integer(from_to[["weights"]])) {
        if(numeric_weights) {
          warning("distance_type is ", distance_type, ", but tr_fun returns integers, to which spaths responds by converting the tr_fun results to numeric")
          from_to[["weights"]] <- as.numeric(from_to[["weights"]])
        }
      } else {
        if(!numeric_weights) {
          warning("distance_type is ", distance_type, ", but tr_fun returns numerics, to which spaths responds by rounding the tr_fun results to integers")
          from_to[["weights"]] <- as.integer(from_to[["weights"]] + 0.5)
        }
      }
    } else {
      if(int_path) {
        if(numeric_weights) {
          if(double_weights) {
            check_weights_i_d(from_to[["from_to"]], from_to[["weights"]])
          } else {
            check_weights_i_f(from_to[["from_to"]], from_to[["weights"]])
          }
        } else {
          if(signed_weights) {
            check_weights_i_i(from_to[["from_to"]], from_to[["weights"]])
          } else {
            check_weights_i_u(from_to[["from_to"]], from_to[["weights"]])
          }
        }
      } else {
        if(numeric_weights) {
          if(double_weights) {
            check_weights_u_d(from_to[["from_to"]], from_to[["weights"]])
          } else {
            check_weights_u_f(from_to[["from_to"]], from_to[["weights"]])
          }
        } else {
          if(signed_weights) {
            check_weights_u_i(from_to[["from_to"]], from_to[["weights"]])
          } else {
            check_weights_u_u(from_to[["from_to"]], from_to[["weights"]])
          }
        }
      }
    }
    rm(tr_fun_args, args_used, tr_fun_v)
  } else {
    if(dist_comp_terra) {
      crd_cell_numbers <- crd[["cell_numbers"]] - 1L
      crd_xy <- matrix(c(rst_xmin + (crd_cell_numbers %% rst_ncol) * rst_xres, rst_ymax - as.integer(crd_cell_numbers / rst_ncol) * rst_yres), ncol = 2L)
      from_to[["weights"]] <- terra::distance(crd_xy[from_to[["from"]],], crd_xy[from_to[["to"]],], lonlat = lonlat, pairwise = TRUE)
      if(output_lines) {
        coords_cnumbers_included <- TRUE
        coords[["cell_numbers"]] <- crd_cell_numbers
      }
      rm(crd_cell_numbers, crd_xy)
      if(!numeric_weights) {
        from_to[["weights"]] <- as.integer(from_to[["weights"]] + 0.5)
      }
    }
    tr_directed <- FALSE
  }
  if(from_to_r && (tr_fun_specified || dist_comp_terra)) {
    from_to[["from"]] <- from_to[["from"]] - 1L
    from_to[["to"]] <- from_to[["to"]] - 1L
  }
  rm(tr_fun, tr_fun_specified, v_matrix, dist_comp_terra, rst_xmin, rst_xres, rst_ymax, rst_yres)
  
  par_lvl_upd <- ncoresg1 && upd_rst_specified && upd_rst_list && match.arg(par_lvl) == "update_rst"
  
  # List cells masked by update_rst (stores non-NA cell ids, not overall cell numbers)
  if(upd_rst_specified) {
    if(show_progress) message("Checking which cells update_rst masks")
    if(rst_terra) {
      if(upd_rst_list) {
        update_rst <- lapply(update_rst, function(V) crd[.(terra::extract(rst, V, cells = TRUE, ID = FALSE, touches = touches)$cell), nomatch = NULL,
          which = TRUE, on = "cell_numbers"] - 1L)
      } else {
        update_rst <- list(crd[.(terra::extract(rst, update_rst, cells = TRUE, ID = FALSE, touches = touches)$cell), nomatch = NULL, which = TRUE,
          on = "cell_numbers"] - 1L)
      }
      rm(rst)
    } else if(upd_rst_matrix) {
      if(upd_rst_list) {
        update_rst <- lapply(update_rst, function(m) crd[.((1:n_cells)[is.na(t(m))]), nomatch = NULL, which = TRUE, on = "cell_numbers"] - 1L)
      } else {
        update_rst <- list(crd[.((1:n_cells)[is.na(t(update_rst))]), nomatch = NULL, which = TRUE, on = "cell_numbers"] - 1L)
      }
    } else {
      if(upd_rst_list) {
        update_rst <- lapply(update_rst, function(v) crd[.(v), nomatch = NULL, which = TRUE, on = "cell_numbers"] - 1L)
      } else {
        update_rst <- list(crd[.(update_rst), nomatch = NULL, which = TRUE, on = "cell_numbers"] - 1L)
      }
    }
  }
  
  n_cells <- nrow(crd)
  
  if((output_lines || !wweights) && !coords_cnumbers_included) {
    coords[["cell_numbers"]] <- crd[["cell_numbers"]] - 1L
  }
  rm(crd)
  
  # Check direction
  if(dest_specified) {
    flip_direction <- !tr_directed
    if(flip_direction) {
      if(origin_nms_specified) {
        flip_direction <- data.table::uniqueN(origins[["cls"]]) > data.table::uniqueN(destinations[["cls"]])
      } else {
        flip_direction <- data.table::uniqueN(origins[["cls"]]) > data.table::uniqueN(destinations[["cls"]])
      }
    }
  } else {
    flip_direction <- FALSE
  }
  
  # Create starts_targets
  if(pairwise) {
    if(origin_nms_specified) {
      origin_names <- origins[["nms"]]
    }
    if(destination_nms_specified) {
      destination_names <- destinations[["nms"]]
    }
    if(flip_direction) {
      starts_targets <- data.table::data.table(starts = destinations[["cls"]], targets = origins[["cls"]])
    } else {
      starts_targets <- data.table::data.table(starts = origins[["cls"]], targets = destinations[["cls"]])
    }
    rm(origins, destinations)
    starts_targets[, path := 1:.N]
    data.table::setkey(starts_targets, starts)
    path_order <- starts_targets[["path"]]
    starts_targets <- c(as.list(starts_targets[, list(n_paths_per_start = .N), by = "starts"]), list(targets = starts_targets[["targets"]]), n_starts = 0L,
      n_targets = length(starts_targets[["targets"]]))
  } else {
    if(origin_nms_specified) {
      origin_names <- origins[["nms"]]
    }
    n_origins <- length(origins[["cls"]])
    if(flip_direction) {
      starts_targets <- list(starts = destinations[["cls"]])
    } else {
      starts_targets <- list(starts = origins[["cls"]])
      starts_targets[["n_starts"]] <- n_origins
    }
    if(dest_specified) {
      n_destinations <- length(destinations[["cls"]])
      if(destination_nms_specified) {
        destination_names <- destinations[["nms"]]
      }
      if(flip_direction) {
        starts_targets[["targets"]] <- origins[["cls"]]
        starts_targets[["n_starts"]] <- n_destinations
        starts_targets[["n_targets"]] <- n_origins
      } else {
        starts_targets[["targets"]] <- destinations[["cls"]]
        starts_targets[["n_targets"]] <- n_destinations
      }
    } else {
      starts_targets[["n_targets"]] <- 0L
    }
    rm(origins, destinations)
  }
  
  # Paths
  if(output_lines) {
    if(upd_rst_specified) {
      if(wweights) {
        paths_lines <- r_upd_paths_wweights(from_to, starts_targets, coords, n_cells, update_rst, early_stopping, ncores, pairwise, tr_directed,
          par_lvl_upd, int_path, numeric_weights, double_weights, signed_weights, return_dists, show_progress, bar_limit, from_to_r)
      } else {
        paths_lines <- r_upd_paths_woweights(from_to, starts_targets, coords, n_cells, update_rst, early_stopping, lonlat, queen, ncores, pairwise,
          par_lvl_upd, pre, int_path, numeric_weights, double_weights, signed_weights, return_dists, show_progress, bar_limit, radius2, from_to_r)
      }
    } else {
      if(wweights) {
        if(from_to_r) {
          paths_lines <- r_paths_wweights(from_to, starts_targets, coords, n_cells, early_stopping, ncores, tr_directed, pairwise, int_path,
            numeric_weights, double_weights, signed_weights, return_dists, show_progress, bar_limit, from_to_r)
        } else {
          
        }
      } else {
        paths_lines <- r_paths_woweights(from_to, starts_targets, coords, n_cells, early_stopping, lonlat, queen, ncores, pairwise, pre, int_path,
          numeric_weights, double_weights, signed_weights, return_dists, show_progress, bar_limit, radius2, from_to_r)
      }
    }
    paths <- paths_lines[setdiff(names(paths_lines), "paths")]
    if(!return_dists) {
      n_paths <- length(paths_lines[["paths"]])
    }
    if(output_spatvector) {
      paths_lines <- terra::vect(paths_lines[["paths"]], "lines", r_crs)
    } else {
      paths_lines <- paths_lines[["paths"]]
    }
    if(return_dists) {
      paths <- data.table::data.table(distances = paths[["distances"]])
    } else {
      paths <- data.table::data.table(connected = rep.int(TRUE, n_paths))[paths[["unconnected_indices"]], connected := FALSE]
    }
  # Distances
  } else {
    # With grid updating
    if(upd_rst_specified) {
      # With precomputed weights
      if(wweights) {
        if(numeric_weights) {
          paths <- data.table::data.table(distances = r_upd_dists_wweights_d(from_to, starts_targets, n_cells, update_rst, early_stopping, ncores, pairwise,
            tr_directed, par_lvl_upd, int_path, double_weights, show_progress, bar_limit, from_to_r))
        } else {
          paths <- data.table::data.table(distances = r_upd_dists_wweights_i(from_to, starts_targets, n_cells, update_rst, early_stopping, ncores, pairwise,
            tr_directed, par_lvl_upd, int_path, signed_weights, show_progress, bar_limit, from_to_r))
        }
      # Without precomputed weights
      } else {
        if(numeric_weights) {
          paths <- data.table::data.table(distances = r_upd_dists_woweights_d(from_to, starts_targets, coords, n_cells, update_rst, lonlat, queen, ncores,
            par_lvl_upd, pairwise, pre, early_stopping, int_path, double_weights, show_progress, bar_limit, radius2, from_to_r))
        } else {
          paths <- data.table::data.table(distances = r_upd_dists_woweights_i(from_to, starts_targets, coords, n_cells, update_rst, lonlat, queen, ncores,
            par_lvl_upd, pairwise, pre, early_stopping, int_path, signed_weights, show_progress, bar_limit, radius2, from_to_r))
        }
      }
    # Without grid updating
    } else {
      if(wweights) {
        if(numeric_weights) {
          paths <- data.table::data.table(distances = r_dists_wweights_d(from_to, starts_targets, n_cells, early_stopping, ncores, tr_directed, pairwise,
            int_path, double_weights, show_progress, bar_limit, from_to_r))
        } else {
          paths <- data.table::data.table(distances = r_dists_wweights_i(from_to, starts_targets, n_cells, early_stopping, ncores, tr_directed, pairwise,
            int_path, signed_weights, show_progress, bar_limit, from_to_r))
        }
      } else {
        if(numeric_weights) {
          paths <- data.table::data.table(distances = r_dists_woweights_d(from_to, starts_targets, coords, n_cells, lonlat, queen, ncores, pairwise, pre,
            early_stopping, int_path, double_weights, show_progress, bar_limit, radius2, from_to_r))
        } else {
          paths <- data.table::data.table(distances = r_dists_woweights_i(from_to, starts_targets, coords, n_cells, lonlat, queen, ncores, pairwise, pre,
            early_stopping, int_path, signed_weights, show_progress, bar_limit, radius2, from_to_r))
        }
      }
    }
  }
  rm(from_to, starts_targets, update_rst)
  if(output_lines || !wweights) rm(coords)
  
  # Add names
  # Pairwise
  if(pairwise) {
    n_paths_layer <- length(path_order)
    path_order <- rep.int(path_order, n_layers) + rep(seq.int(0L, by = n_paths_layer, length.out = n_layers), each = n_paths_layer)
    if(origin_nms_specified) {
      if(destination_nms_specified) {
        paths[path_order, c("origins", "destinations") := list(rep.int(origin_names, n_layers), rep.int(destination_names, n_layers))]
      } else {
        paths[path_order, c("origins", "destinations") := list(rep.int(origin_names, n_layers), rep.int(1:n_paths_layer, n_layers))]
      }
    } else {
      if(destination_nms_specified) {
        paths[path_order, c("origins", "destinations") := list(rep.int(1:n_paths_layer, n_layers), rep.int(destination_names, n_layers))]
      } else {
        paths[path_order, c("origins", "destinations") := list(rep.int(1:n_paths_layer, n_layers), rep.int(1:n_paths_layer, n_layers))]
      }
    }
    rm(path_order)
    if(flip_direction) {
      data.table::setnames(paths, c("origins", "destinations"), c("destinations", "origins"))
      data.table::setcolorder(paths, c("origins", "destinations"))
    }
  # Not pairwise
  } else if(dest_specified) {
    if(flip_direction) {
      if(origin_nms_specified) {
        if(destination_nms_specified) {
          paths[, c("origins", "destinations") := list(rep.int(rep.int(origin_names, n_destinations), n_layers), rep.int(rep(destination_names,
            each = n_origins), n_layers))]
        } else {
          paths[, c("origins", "destinations") := list(rep.int(rep.int(origin_names, n_destinations), n_layers), rep.int(rep(1:n_destinations,
            each = n_origins), n_layers))]
        }
      } else {
        if(destination_nms_specified) {
          paths[, c("origins", "destinations") := list(rep.int(rep.int(1:n_origins, n_destinations), n_layers), rep.int(rep(destination_names,
            each = n_origins), n_layers))]
        } else {
          paths[, c("origins", "destinations") := list(rep.int(rep.int(1:n_origins, n_destinations), n_layers), rep.int(rep(1:n_destinations,
            each = n_origins), n_layers))]
        }
      }
    } else {
      if(origin_nms_specified) {
        if(destination_nms_specified) {
          paths[, c("origins", "destinations") := list(rep.int(rep(origin_names, each = n_destinations), n_layers), rep.int(rep.int(destination_names,
            n_origins), n_layers))]
        } else {
          paths[, c("origins", "destinations") := list(rep.int(rep(origin_names, each = n_destinations), n_layers), rep.int(rep.int(1:n_destinations,
            n_origins), n_layers))]
        }
      } else {
        if(destination_nms_specified) {
          paths[, c("origins", "destinations") := list(rep.int(rep(1:n_origins, each = n_destinations), n_layers), rep.int(rep.int(destination_names,
            n_origins), n_layers))]
        } else {
          paths[, c("origins", "destinations") := list(rep.int(rep(1:n_origins, each = n_destinations), n_layers), rep.int(rep.int(1:n_destinations,
            n_origins), n_layers))]
        }
      }
    }
  # No targets
  } else {
    if(tr_directed) {
      if(origin_nms_specified) {
        paths[, c("origins", "destinations") := list(rep.int(rep(origin_names, each = n_origins - 1L), n_layers), rep.int(rep.int(origin_names,
          n_origins)[-seq.int(1L, n_origins^2L, n_origins + 1L)], n_layers))]
      } else {
        paths[, c("origins", "destinations") := list(rep.int(rep(1:n_origins, each = n_origins - 1L), n_layers), rep.int(rep.int(1:n_origins,
          n_origins)[-seq.int(1L, n_origins^2L, n_origins + 1L)], n_layers))]
      }
    } else {
      if(origin_nms_specified) {
        paths[, origins := rep.int(rep.int(origin_names, (n_origins - 1L):0), n_layers)]
        if(is.character(origin_names)) {
          paths[, destinations := rep.int(destination_names_character(origin_names), n_layers)]
        } else if(is.integer(origin_names)) {
          paths[, destinations := rep.int(destination_names_integer(origin_names), n_layers)]
        } else {
          paths[, destinations := rep.int(destination_names_numeric(origin_names), n_layers)]
        }
      } else {
        paths[, c("origins", "destinations") := list(rep.int(rep.int(1:n_origins, (n_origins - 1L):0), n_layers), rep.int(destination_names_auto(n_origins),
          n_layers))]
      }
    }
  }
  if(origin_nms_specified) rm(origin_names)
  if(dest_specified && destination_nms_specified) rm(destination_names)
  
  # Layer numbers
  if(upd_rst_specified) {
    paths[, layer := rep(0:(n_layers - 1L), each = nrow(paths) / n_layers)]
  }
  
  # Bidirectional
  if(bidirectional && !dest_specified && !tr_directed) {
    paths <- rbind(paths, data.table::setcolorder(data.table::setnames(data.table::copy(paths), c("origins", "destinations"), c("destinations", "origins")),
      names(paths)))
    if(output_lines) {
      if(output_spatvector) {
        paths_lines <- rbind(paths_lines, paths_lines)
      } else {
        paths_lines <- c(paths_lines, paths_lines)
      }
    }
  }
  
  data.table::setcolorder(paths, c("origins", "destinations"))
  
  # Output
  if(output_lines) {
    if(output_spatvector) {
      terra::values(paths_lines) <- paths
    } else {
      paths_lines <- list(attributes = paths, lines = paths_lines, crs = r_crs)
    }
    return(paths_lines)
  } else {
    if(output_datatable) {
      return(paths[])
    } else {
      data.table::setDF(paths)
      return(paths)
    }
  }
}

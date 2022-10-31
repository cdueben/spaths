#' Shortest Paths between Geographic Locations on Earth
#' 
#' The function computes the shortest paths between locations on planet Earth using Dijstra's algorithm.
#' 
#' @param rst SpatRaster (terra) or RasterLayer (raster) object. RasterLayers are converted to SpatRasters. Pixels with non-NA values in all layers mark 
#' the cells through which the algorithm may pass.
#' @param origins A single SpatVector (terra), sf (sf), or Spatial* (sp) object, or a list of them. Non-SpatVectors are converted to SpatVectors. 
#' Polygons and lines are converted to points using their centroid. Passing a list of spatial objects is a way not to relate all origins to all 
#' destinations. Details are outlined below.
#' @param destinations A single SpatVector (terra), sf (sf), or Spatial* (sp) object, a list of them, or \code{NULL}. It defaults to \code{NULL}, 
#' resulting in the function to compute shortest paths between the \code{origins} points. Non-SpatVectors are converted to SpatVectors. Polygons and 
#' lines are converted to points using their centroid. Passing a list of spatial objects is a way not to relate all origins to all destinations. Details 
#' are outlined below.
#' @param output \code{"lines"} (default) or \code{"distances"}. \code{"lines"} returns the shortest paths as SpatVector lines. \code{"distances"} lists 
#' the total transition costs along the shortest paths. By default, it is the distance between origin and destination in meters. If you pass another 
#' function to \code{tr_fun}, the total transition cost is measured in the units of \code{tr_fun}'s results. \code{"distances"} is faster and requires 
#' less RAM than \code{"lines"}.
#' @param origin_names Character specifying the name of the column in the \code{origins} object used to label the origins in the output object. It 
#' defaults to row numbers, which are more efficient than character labels.
#' @param destination_names Character specifying the name of the column in the \code{destinations} object used to label the origins in the output object. 
#' It defaults to row numbers, which are more efficient than character labels.
#' @param pairwise Logical specifying whether to compute pairwise paths, if \code{origins} and \code{destinations} have equally many rows. If \code{TRUE}, 
#' the function computes the shortest path between the first origin and the first destination, the second origin and the second destination etc. 
#' Otherwise, it derives the shortest paths from all origins to all destinations.
#' @param contiguity "queen" (default) for queen's case contiguity or "rook" for rook's case contiguity. In the latter case, the algorithm is only 
#' allowed to move between horizontally or vertically adjacent cells. In the former case, it is also allowed to travel between diagonally adjacent cells. 
#' "rook" is more efficient than "queen" as it implies fewer edges.
#' @param dist_comp \code{"spaths"} (default) or \code{"terra"}. The former uses spherical (Haversine) distances in case of lonlat data and planar 
#' (Euclidean) distances in case of projected (non-lonlat) data. The functions are optimized based on the fact that many inter-pixel distances are 
#' identical. Modelling the planet as a perfect sphere is in line with e.g. the s2 package, but is of course an oversimplification. With \code{"terra"}, 
#' the function derives distances via \code{terra::distance}. Because this computes all inter-pixel distances separately, it is slower than the 
#' \code{"spaths"} approach. It does take the non-spherical nature of the planet into account though. And as there is no guarantee that the \code{"spaths"} 
#' Euclidean function get the distances right in case of all projections, \code{"terra"} is the safer option for non-lonlat data.
#' @param tr_fun The transition function based on which to compute edge weights, i.e. the travel cost between grid cells. Defaults to the geographic 
#' distance between pixel centroids. Permitted function parameter names are \code{d} (distance between the pixel centroids), \code{x1} (x coordinate or 
#' longitude of the first cell), \code{x2} (x coordinate or longitude of the second cell), \code{y1} (y coordinate or latitude of the first cell), 
#' \code{y2} (y coordinate or latitude of the second cell), \code{v1} (\code{rst} layers' values from the first cell), \code{v2} (\code{rst} layers' 
#' values from the second cell), \code{nc} (number of CPU cores), and \code{cl} (cluster object returned by \code{parallel::makePSOCKcluster(ncores)} or 
#' \code{parallel::makeForkCluster(ncores)}). If \code{lonlat} is \code{TRUE} or if \code{dist_comp = "terra"}, \code{d} is measured in the meters. 
#' Otherwise, it uses the units of the CRS. \code{nc} is meant to be used in C++ functions. A parallel backend at the R level pointed to by \code{cl} is 
#' already registered before \code{tr_fun} is called, in case \code{ncores > 1}. If \code{rst} has one layer, the values are passed to \code{v1} and 
#' \code{v2} as vectors, otherwise they are passed as a data table where the first column refers to the first layer, the second column to the second 
#' layer etc. Note that data tables are also data frames.
#' @param v_matrix Logical specifying whether to pass values to \code{v1} and \code{v2} in \code{tr_fun} as matrices (\code{TRUE}) instead of data tables 
#' in the multi-layer case and vectors in the single-layer case (\code{FALSE}). It defaults to \code{FALSE}. Setting it to \code{TRUE} might e.g. be 
#' useful when defining \code{tr_fun} as a C++ Armadillo function.
#' @param tr_directed Logical specifying whether \code{tr_fun} creates a directed graph. In a directed graph, transition costs can be asymmetric. 
#' Traveling from cells A to B may imply a different cost than traveling from B to A. It defaults to \code{TRUE} and only has an effect when 
#' \code{tr_fun} is not \code{NULL}. The default without \code{tr_fun} constructs an undirected graph.
#' @param update_rst A single SpatVector (terra), sf (sf), or Spatial* (sp) object, a list of them, or \code{NULL}. It defaults to \code{NULL}, 
#' corresponding to \code{rst} not being updated. sf and sp objects are converted to SpatVectors. If \code{update_rst} is a single vector object, the 
#' shortest paths are computed twice. The function first derives the paths based on the specified \code{rst}, then removes all cells intersecting with 
#' the vector object, and derives the paths based on the updated set of cells. The result is a list of length two, where the first element entails the 
#' paths from the unmodified \code{rst} and the second element the paths from the updated \code{rst}. If \code{update_rst} is a list, the function 
#' produces a list that is one element longer than \code{update_rst} is. The first element again consists of shortest paths based on the unupdated 
#' \code{rst}, while the subsequent elements accomodate paths derived from the updated \code{rst}s. The updates always start from the original \code{rst}. 
#' I.e. \code{update_rst[[2]]} does not modify the \code{rst} updated by \code{update_rst[[1]]}, but the \code{rst} unadjusted by any element of 
#' \code{update_rst}. Note that \code{update_rst} only removes cells. Any pixels with NA values in the original \code{rst} remain excluded from the set 
#' of transition cells. Neither does \code{update_rst} alter cell values from any of \code{rst}'s layers.
#' @param touches Logical specifying the \code{touches} argument of \code{terra::extract} used when \code{update_rst} is not \code{NULL}. It defaults to 
#' \code{TRUE}. If \code{FALSE}, the function only removes cells on the line render path or with the center point inside a polygon.
#' @param copy Logical specifying whether to copy paths that are unaffected by \code{update_rst} geometries (\code{TRUE}) rather than re-estimating all 
#' paths on the grids updated by \code{update_rst} (\code{FALSE}). It defaults to \code{TRUE} and only has an effect when \code{update_rst} is not 
#' \code{NULL}. If \code{TRUE}, the function first computes paths on the unmodified \code{rst} and then checks which paths' cells intersect with the 
#' geometries in \code{update_rst}. If you know that all paths are affected by \code{update_rst}, setting \code{copy} to \code{FALSE} can be the more 
#' efficient choice. This especially holds in two scenarios: (i) if you use multiple cores and \code{par_lvl = "update_rst"}, \code{copy = FALSE} directly 
#' parallelizes the function over the different grids, incl. the unmodified \code{rst}, rather than waiting for the paths on the initial grid to be 
#' computed before parallelizing over the remaining grids; (ii) if \code{output = "distances"}, the function only computes distances when 
#' \code{copy = FALSE}, but engages in the more time consuming task of returning lines for the first grid when \code{copy = TRUE}. Thus, the function 
#' tends to be faster with \code{copy = TRUE} than with \code{copy = FALSE}, unless quasi all paths are affected by \code{update_rst}, with both 
#' \code{output = "lines"} and \code{output = "distances"}.
#' @param round_dist Logical specifying whether to round edge weights, i.e. the transition cost or distance between neighboring cells, to integers. It 
#' defaults to \code{FALSE}. Setting it to \code{TRUE} reduces the function's RAM requirements slightly.
#' @param ncores An integer specifying the number of CPU cores to use. It defaults to the number of cores installed on the machine. A value of 1 
#' accordingly induces a single-threaded process.
#' @param par_lvl \code{"points"} (default), \code{"points_lists"}, or \code{"update_rst"}, indicating the level at which to parallelize when using 
#' multiple cores. \code{"points"} parallelizes over the origin (and destination) point combinations. \code{"points_lists"} parallelizes over the list 
#' elements of \code{origins} (and \code{destinations}), if these arguments are lists. \code{"update_rst"} parallelizes over the list of graphs specified 
#' by \code{rst} and \code{update_rst}.
#' @param rst_par_lvl \code{"points"} (default), \code{"points_lists"}, or \code{"none"}, indicating the level at which to parallelize the paths 
#' computations on the unmodified \code{rst} when using multiple cores, \code{update_rst} is not \code{NULL}, \code{par_lvl = "update_rst"}, and 
#' \code{copy = TRUE}. Because the function requires the paths from the unmodified \code{rst} before looping over the updated grids when 
#' \code{copy = TRUE}, the parallel structure requested by \code{par_lvl = "update_rst"} only commences after those initial paths have been identified. 
#' \code{rst_par_lvl} defines the parallel level at which those initial paths are derived. \code{"points"} (default) and \code{"points_list"} are 
#' equivalent to their counterparts in \code{par_lvl}, \code{"none"} induces a serial, i.e. non-parallel, execution.
#' @param cluster \code{"PSOCK"} or \code{"FORK"}, indicating the type of \code{parallel} cluster that the function employs when \code{ncores > 1} or 
#' \code{paths_ncores > 1}. The function defaults to \code{"PSOCK"} on Windows and to \code{"FORK"} otherwise. Windows machines must use \code{"PSOCK"}, 
#' while Mac and Linux can employ either of the two options.
#' @param paths_ncores An integer specifying the number of CPU cores to use in shortest paths computations. It defaults to the value of \code{ncores}. 
#' Thus, only set it, if you want edge weights and shortest paths be computed with differently many cores. The \code{dist_comp = "spaths"} edge weight
#' computations employ efficient C++ level parallelization. The shortest paths sections, in contrast, parallelize on the R level. If you use a \code{PSOCK} 
#' cluster, \code{spaths_earth} copies various objects to the workers before the paths algorithm is applied. This can make the parallel execution slower 
#' than its serial counterpart. Thus, consider setting \code{paths_ncores = 1}, especially when working with \code{PSOCK} clusters.
#' @param verbose Logical specifying whether the function prints messages on progress. It defaults to \code{FALSE}.
#' 
#' @details This function computes shortest paths between locations on Earth, taking custom transition costs into account. Examples are a ship navigating 
#' around land masses or a hiker traversing mountains.
#' 
#' The algorithm links \code{origins} and \code{destinations} by passing through the centroids of \code{rst}'s non-NA cells.
#' 
#' The \href{../doc/spaths_introduction.html}{vignette} provides further details on the function.
#' 
#' @return In the basic cases, i.e. when neither \code{origins} nor \code{destinations} are lists and \code{update_rst} is not specified, the function 
#' returns a SpatVector lines object with \code{output = "lines"} and a data table with \code{output = "distances"}. If \code{origins} or 
#' \code{destinations} are lists or \code{update_rst} is specified, it returns a list or nested list of SpatVectors or data tables. In the nested case, 
#' the outer list refers to the \code{update_rst} level and the inner list to the \code{origins} and \code{destinations} level. Consult the 
#' \href{../doc/spaths_introduction.html}{vignette} for further details.
#' 
#' The distances are the total transition costs along the shortest paths. If \code{lonlat} is \code{TRUE} or if \code{dist_comp = "terra"}, these 
#' distances are measured in the meters. Otherwise, they use the units of the CRS. If you pass a function to \code{tr_fun}, the total transition cost is 
#' measured in the units of \code{tr_fun}'s results.
#' 
#' @seealso \link{spaths_general}
#' 
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(20146)
#' rst <- terra::rast(crs = "epsg:4326", resolution = 1L, vals = sample(c(1L, NA), 64800L, TRUE,
#'   c(0.9, 0.1)))
#' origins <- rnd_locations(5L, output_type = "terra")
#' destinations <- rnd_locations(5L, output_type = "terra")
#' 
#' # Compute shortest paths
#' spaths_earth(rst, origins)
#' spaths_earth(rst, origins, destinations)
#' }
#' 
#' @importFrom data.table :=
#' @importFrom data.table .N
#' @importFrom data.table .SD
#' @importFrom data.table %chin%
#' @importFrom data.table %between%
#' 
#' @useDynLib spaths, .registration = TRUE
#' 
#' @export
spaths_earth <- function(rst, origins, destinations = NULL, output = c("lines", "distances"), origin_names = NULL, destination_names = NULL, pairwise = FALSE,
  contiguity = c("queen", "rook"), dist_comp = c("spaths", "terra"), tr_fun = NULL, v_matrix = FALSE, tr_directed = TRUE, update_rst = NULL, touches = TRUE,
  copy = TRUE, round_dist = FALSE, ncores = NULL, par_lvl = c("points", "points_lists", "update_rst"), rst_par_lvl = c("points", "points_lists", "none"),
  cluster = NULL, paths_ncores = NULL, verbose = FALSE) {
  
  if(length(verbose) != 1L || !is.logical(verbose) || is.na(verbose)) stop("verbose must be logical and of length one")
  if(verbose) message("Checking arguments")
  
  # Convert RasterLayer to SpatRaster
  if(all(class(rst) != "SpatRaster")) rst <- terra::rast(rst)
  
  # Check CRS
  r_crs <- terra::crs(rst)
  if(r_crs == "") stop("rst without specified CRS")
  lonlat <- terra::is.lonlat(rst)
  
  # Check origin_names and destination_names
  origin_nms_specified <- !is.null(origin_names)
  if(origin_nms_specified && (!(is.character(origin_names) && length(origin_names) == 1L) || is.na(origin_names))) {
    stop("origin_names must be NULL or a character object of length one")
  }
  destination_nms_specified <- !is.null(destination_names)
  if(destination_nms_specified && (!(is.character(destination_names) && length(destination_names) == 1L) || is.na(destination_names))) {
    stop("destination_names must be NULL or a character object of length one")
  }
  
  # Convert origins
  origin_list <- any(class(origins) == "list")
  if(origin_list) {
    origins <- lapply(origins, convert_points, rst, r_crs, origin_names, origin_nms_specified)
  } else {
    origins <- convert_points(origins, rst, r_crs, origin_names, origin_nms_specified)
  }
  
  # Convert destinations
  dest_specified <- !is.null(destinations)
  if(dest_specified) {
    dest_list <- any(class(destinations) == "list")
    if(dest_list) {
      destinations <- lapply(destinations, convert_points, rst, r_crs, destination_names, destination_nms_specified, FALSE)
    } else {
      destinations <- convert_points(destinations, rst, r_crs, destination_names, destination_nms_specified, FALSE)
    }
  } else {
    dest_list <- NULL
  }
  
  # Check other arguments
  if(length(pairwise) != 1L || !is.logical(pairwise) || is.na(pairwise)) stop("pairwise must be logical and of length one")
  if(pairwise) {
    if(dest_specified) {
      if(origin_list) {
        if(dest_list) {
          if(length(origins) != length(destinations)) stop("If both origins and destinations are lists, they must be lists of equal length")
          if(!all(vapply(origins, NROW, integer(1L), USE.NAMES = FALSE) == vapply(destinations, NROW, integer(1L), USE.NAMES = FALSE))) {
            stop("origins and destinations must have the same number of rows, if pairwise is TRUE")
          }
        } else if(!all(vapply(origins, NROW, integer(1L), USE.NAMES = FALSE) == NROW(destinations))) {
          stop("origins and destinations must have the same number of rows, if pairwise is TRUE")
        }
      } else if(dest_list) {
        if(!all(NROW(origins) == vapply(destinations, NROW, integer(1L), USE.NAMES = FALSE))) {
          stop("origins and destinations must have the same number of rows, if pairwise is TRUE")
        }
      } else if(NROW(origins) != NROW(destinations)) {
        stop("origins and destinations must have the same number of rows, if pairwise is TRUE")
      }
    } else {
      stop("destinations must not be NULL, if pairwise is TRUE")
    }
  } else {
    if(dest_specified) {
      if(origin_list) {
        if(dest_list) {
          if(length(origins) != length(destinations)) stop("If both origins and destinations are lists, they must be lists of equal length")
        } else {
          stop("If pairwise is FALSE and origins is a list, destinations must also be a list. Pass origins and destinations as single vector objects, ",
            "not as lists, to compute shortest paths from all origins to all destinations.")
        }
      } else if(dest_list) {
        stop("If pairwise is FALSE and destinations is a list, origins must also be a list. Pass origins and destinations as single vector objects, not as ",
          "lists, to compute shortest paths from all origins to all destinations.")
      }
    } else {
      if(origin_list) {
        if(!all(vapply(origins, function(O) NROW(O) > 1L, logical(1L), USE.NAMES = FALSE))) {
          stop("Each origins list element must have more than one row, if destinations are not specified")
        }
      } else if(NROW(origins) <= 1L) {
        stop("origins must have more than one row, if destinations are not specified")
      }
    }
  }
  output <- match.arg(output)
  contiguity <- match.arg(contiguity)
  dist_comp <- match.arg(dist_comp)
  tr_fun_specified <- !is.null(tr_fun)
  if(tr_fun_specified) {
    if(is.function(tr_fun)) {
      tr_fun_v <- names(formals(tr_fun))                                        # Extract function parameter names
      if(!all(tr_fun_v %chin% c("d", "x1", "x2", "y1", "y2", "v1", "v2", "nc", "cl"))) {
        stop("The permitted parameter names in tr_fun are d, x1, x2, y1, y2, v1, v2, nc, and cl")
      }
      args_used <- c("d", "x1", "y1", "x2", "y2", "v1", "v2", "nc", "cl") %chin% tr_fun_v
    } else {
      stop("tr_fun must be either NULL or a function")
    }
    if(length(v_matrix) != 1L || !is.logical(v_matrix) || is.na(v_matrix)) stop("v_matrix must be logical and of length one")
    if(length(tr_directed) != 1L || !is.logical(tr_directed) || is.na(tr_directed)) stop("tr_directed must be logical and of length one")
  }
  upd_rst_specified <- !is.null(update_rst)
  if(upd_rst_specified) {
    update_rst_list <- any(class(update_rst) == "list")
    if(update_rst_list) {
      update_rst <- lapply(update_rst, function(v) {
        if(all(class(v) != "SpatVector")) v <- terra::vect(v)
        if(terra::crs(v) != r_crs) v <- terra::project(v, r_crs)
        return(v)
      })
    } else {
      if(all(class(update_rst) != "SpatVector")) update_rst <- terra::vect(update_rst)
      if(terra::crs(update_rst) != r_crs) update_rst <- terra::project(update_rst, r_crs)
    }
    if(length(touches) != 1L || !is.logical(touches) || is.na(touches)) stop("touches must be logical and of length one")
    if(length(copy) != 1L || !is.logical(copy) || is.na(copy)) stop("copy must be logical and of length one")
    if(copy) rst_par_lvl <- match.arg(rst_par_lvl)
  }
  if(length(round_dist) != 1L || !is.logical(round_dist) || is.na(round_dist)) stop("round_dist must be logical and of length one")
  if(is.null(ncores)) {
    ncores <- parallel::detectCores()
  } else if(length(ncores) != 1L || ncores < 1) {
    stop("ncores must be either NULL or a positive integer")
  }
  paths_ncores_specified <- !is.null(paths_ncores)
  if(paths_ncores_specified && (length(paths_ncores) != 1L || paths_ncores < 1)) {
    stop("paths_ncores must either be NULL or a positive integer")
  }
  ncoresg1 <- ncores > 1
  if(ncoresg1 || (paths_ncores_specified && paths_ncores > 1)) {
    if(is.null(cluster)) {
      if(.Platform$OS.type == "windows") {
        cluster <- "PSOCK"
      } else {
        cluster <- "FORK"
      }
    } else if(length(cluster) != 1L || !(cluster %chin% c("PSOCK", "FORK"))) {
      stop('cluster must be NULL, "PSOCK", or "FORK"')
    }
  }
  # par_lvl is checked before shortest paths computation
    
  if(verbose) message("Constructing graph")
  
  # List non-NA cells as data.table
  if(tr_fun_specified && any(args_used[6:7])) {
    if(any(c("x", "y", "c_n") %chin% names(rst))) {
      p_names <- c("x", "y", "c_n") %chin% names(rst)
      if(sum(p_names) > 1L) {
        if(sum(p_names) > 2L) {
          p_names <- "x, y, and c_n"
        } else {
          p_names <- paste0(c("x", "y", "c_n")[p_names], collapse = " and ")
        }
        p_names <- paste0(p_names, " are not permitted rst layer names")
      } else {
        p_names <- paste0(c("x", "y", "c_n")[p_names], " is not a permitted rst layer name")
      }
      stop(paste0(p_names, ", when using a transition function with parameter v"))
    }
    crd <- terra::values(rst, dataframe = TRUE)                                 # Extract layer values
    data.table::setDT(crd)
    crd[, c_n := 1:.N]
    crd <- stats::na.omit(crd, cols = 1:(NCOL(crd) - 1L))                       # Subset to cells with non-NA values in all layers
    crd[, c("x", "y") := data.table::as.data.table(terra::xyFromCell(rst, crd[["c_n"]]))] # Obtain coordinates of non-NA cells
  } else {
    if(terra::nlyr(rst) > 1) {
      crd <- (1:terra::ncell(rst))[stats::complete.cases(terra::values(rst))]   # Obtain ids of non-NA cells
    } else {
      crd <- (1:terra::ncell(rst))[!is.na(terra::values(rst))]                  # Obtain ids of non-NA cells
    }
    crd <- data.table::as.data.table(terra::xyFromCell(rst, crd))[, c_n := crd] # Obtain coordinates of non-NA cells
  }
  data.table::setkey(crd, "c_n")
  
  # Obtain various rst attributes used in distance computation
  if(dist_comp == "spaths") {
    yr <- terra::yres(rst)
    xr <- terra::xres(rst)
    nr <- terra::nrow(rst)
    ym <- terra::ymin(rst)
  } else {
    yr <- NULL
    xr <- NULL
    nr <- NULL
    ym <- NULL
  }
  
  if(upd_rst_specified) rst_upd <- rst[[1L]]
  crd[, c_n_c := 1:.N]
  
  origins <- update_points(origins, origin_list, crd, origin_nms_specified)
  if(dest_specified) destinations <- update_points(destinations, dest_list, crd, destination_nms_specified)
  
  rst <- data.table::as.data.table(terra::adjacent(rst, crd[["c_n"]], contiguity, TRUE)) # Obtain adjacency matrix
  if(tr_fun_specified && tr_directed) {
    rst <- rst[crd[, c("c_n", "c_n_c")], nomatch = NULL, on = "to==c_n"][, to := NULL] # Subset to non-NA destination cells
  } else {
    rst <- rst[crd[, c("c_n", "c_n_c")], nomatch = NULL, on = "to==c_n"][from < to,][, to := NULL] # Subset to non-NA destination cells and one edge per pair
  }
  data.table::setnames(rst, "c_n_c", "to")
  rst <- rst[crd[, c("c_n", "c_n_c")], nomatch = NULL, on = "from==c_n"][, from := NULL]
  data.table::setnames(rst, "c_n_c", "from")
  if(upd_rst_specified) {
    crd[, c_n_c := NULL]
  } else {
    crd[, c("c_n", "c_n_c") := NULL]
  }
  
  # Register parallel backend
  dist_comp_terra <- dist_comp == "terra"
  if(ncoresg1 && ((tr_fun_specified && ((args_used[1L] && dist_comp_terra) || args_used[9L])) || (!tr_fun_specified && dist_comp_terra))) {
    if(cluster == "PSOCK") {
      cl <- parallel::makePSOCKcluster(ncores)
    } else {
      cl <- parallel::makeForkCluster(ncores)
    }
  } else {
    cl <- NULL
  }
  
  if(tr_fun_specified) {
    # Collect transition function arguments
    tr_fun_args <- list()
    if(args_used[1L]) {
      tr_fun_args$d <- compute_dists(rst, crd[, c("x", "y")], dist_comp_terra, round_dist, contiguity, yr, xr, nr, ym, lonlat, ncoresg1, ncores, cl)
    }
    if(args_used[2L]) tr_fun_args$x1 <- crd[rst[["from"]], "x"][["x"]]
    if(args_used[3L]) tr_fun_args$y1 <- crd[rst[["from"]], "y"][["y"]]
    if(args_used[4L]) tr_fun_args$x2 <- crd[rst[["to"]], "x"][["x"]]
    if(args_used[5L]) tr_fun_args$y2 <- crd[rst[["to"]], "y"][["y"]]
    if(args_used[6L]) {
      if(v_matrix) {
        tr_fun_args$v1 <- as.matrix(crd[rst[["from"]], .SD, .SDcols = setdiff(names(crd), c("x", "y", "c_n"))])
      } else {
        v_vars <- setdiff(names(crd), c("x", "y", "c_n"))
        if(length(v_vars) > 1L) {
          tr_fun_args$v1 <- crd[rst[["from"]], .SD, .SDcols = v_vars]
        } else {
          tr_fun_args$v1 <- crd[rst[["from"]], .SD, .SDcols = v_vars][[1L]]
        }
        rm(v_vars)
      }
    }
    if(args_used[7L]) {
      if(v_matrix) {
        tr_fun_args$v2 <- as.matrix(crd[rst[["to"]], .SD, .SDcols = setdiff(names(crd), c("x", "y", "c_n"))])
      } else {
        v_vars <- setdiff(names(crd), c("x", "y", "c_n"))
        if(length(v_vars) > 1L) {
          tr_fun_args$v2 <- crd[rst[["to"]], .SD, .SDcols = v_vars]
        } else {
          tr_fun_args$v2 <- crd[rst[["to"]], .SD, .SDcols = v_vars][[1L]]
        }
        rm(v_vars)
      }
    }
    if(any(args_used[6:7])) crd[, setdiff(names(crd), c("x", "y", "c_n")) := NULL]
    if(args_used[8L]) tr_fun_args$nc <- ncores
    if(args_used[9L]) tr_fun_args$cl <- cl
    tr_fun_args <- do.call(tr_fun, tr_fun_args[tr_fun_v])
    if(!is.vector(tr_fun_args)) {
      if(!is.null(cl)) parallel::stopCluster(cl)
      stop("tr_fun must return a vector")
    }
    if(any(tr_fun_args < 0)) {
      if(!is.null(cl)) parallel::stopCluster(cl)
      stop("tr_fun must not return negative values")
    }
    if(length(tr_fun_args) != NROW(rst)) stop("The number of values returned by tr_fun must equal the number of edges")
    rst <- igraph::set_edge_attr(igraph::graph_from_edgelist(as.matrix(rst), directed = tr_directed), "weight", value = tr_fun_args) # Construct graph
    rm(tr_fun_args, args_used, tr_fun_v)
  } else {
    rst <- igraph::set_edge_attr(igraph::graph_from_edgelist(as.matrix(rst), directed = FALSE), "weight", value = compute_dists(rst, crd[, c("x", "y")],
      dist_comp_terra, round_dist, contiguity, yr, xr, nr, ym, lonlat, ncoresg1, ncores, cl)) # Construct graph
  }
  
  if(verbose) message("Computing shortest paths")
  
  # Register parallel backend
  if(paths_ncores_specified && ncores != paths_ncores) {
    ncores <- paths_ncores
    ncoresg1 <- ncores > 1
    if(!is.null(cl)) {
      parallel::stopCluster(cl)
      cl <- NULL
    }
  }
  if(ncoresg1) {
    par_lvl <- match.arg(par_lvl)
    if(par_lvl == "update_rst") {
      if(!upd_rst_specified) {
        if(!is.null(cl)) parallel::stopCluster(cl)
        stop('If par_lvl is "update_rst", update_rst must not be NULL')
      }
      if(copy && !update_rst_list) {
        if(!is.null(cl)) parallel::stopCluster(cl)
        stop('If copy is TRUE and update_rst is not a list, par_lvl = "update_rst" does not induce parallelization')
      }
    }
    if(!origin_list && dest_specified && !dest_list && par_lvl == "points_lists") {
      if(!is.null(cl)) parallel::stopCluster(cl)
      stop('If par_lvl is "points_lists", origins or destinations must be a list')
    }
    if(is.null(cl)) {
      if(cluster == "PSOCK") {
        cl <- parallel::makePSOCKcluster(ncores)
      } else {
        cl <- parallel::makeForkCluster(ncores)
      }
    }
  } else {
    par_lvl <- ""
  }
  
  rm(yr, xr, nr, ym, tr_fun_specified, paths_ncores_specified, dist_comp_terra)
  
  # Compute shortest paths
  if(upd_rst_specified) {
    output_lines <- output == "lines"
    if(update_rst_list) {
      update_rst <- lapply(update_rst, function(V) crd[.(terra::extract(rst_upd, V, cells = TRUE, ID = FALSE, touches = touches)$cell), nomatch = NULL,
        which = TRUE, on = "c_n"])
    } else {
      update_rst <- crd[.(terra::extract(rst_upd, update_rst, cells = TRUE, ID = FALSE, touches = touches)$cell), nomatch = NULL, which = TRUE, on = "c_n"]
    }
    rm(rst_upd)
    crd[, c_n := NULL]
    if((update_rst_list && any(lengths(update_rst) > 0L)) || (!update_rst_list && length(update_rst) > 0L)) {
      p_list <- origin_list || (dest_specified && dest_list)
      if(copy) {
        if(rst_par_lvl == "none") {
          p1 <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list, dest_list,
            r_crs, output_lines, pairwise, FALSE, NULL, NULL, NULL, TRUE)
        } else {
          p1 <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list, dest_list,
            r_crs, output_lines, pairwise, ncoresg1, ncores, rst_par_lvl, cl, TRUE)
        }
        if(p_list) {
          p2 <- lapply(p1, `[[`, 2L)
          if(output_lines) p3 <- lapply(p1, `[[`, 3L)
          p1 <- lapply(p1, `[[`, 1L)
        } else {
          p2 <- p1[[2L]]
          if(output_lines) p3 <- p1[[3L]]
          p1 <- p1[[1L]]
        }
        # Lines output
        if(output_lines) {
          if(p_list) {
            paths <- function(u) {
              p_affected <- lapply(p1, function(p) unique(p[.(u), "g", nomatch = NULL, on = "cls"][["g"]]))
              if(any(lengths(p_affected) > 0L)) {
                rst_u <- igraph::delete_vertices(rst, u)
                v <- data.table::data.table(c_n_c = 1:NROW(crd))[-u,]
                p <- function(P_AFFECTED, P1, P2, P3 = NULL) {
                  if(length(P_AFFECTED) > 0L) {
                    origin_c <- P2[P_AFFECTED, "origin_c"]
                    destination_c <- P2[P_AFFECTED, "destination_c"]
                    if(any(origin_c %in% u)) {
                      report_points_ust(origin_c, u, TRUE, cl)
                    }
                    if(any(destination_c %in% u)) {
                      report_points_ust(destination_c, u, FALSE, cl)
                    }
                    origin_c <- v[.(origin_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                    destination_c <- v[.(destination_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                    P <- crd[-u, c("x", "y")]
                    if(ncoresg1 && par_lvl == "points") {
                      P <- rbind(P1[!.(P_AFFECTED), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::clusterMap(cl, function(oc, dc, pa) {
                        return(P[igraph::shortest_paths(rst_u, oc, dc, output = "vpath", algorithm = "dijkstra")$vpath[[1L]],][, g := pa])
                      }, origin_c, destination_c, P_AFFECTED, USE.NAMES = FALSE), use.names = FALSE)[, c("g", "x", "y")])
                      data.table::setorder(P, g)
                      P <- terra::vect(as.matrix(P), type = "line", atts = P3, crs = r_crs)
                    } else {
                      P <- rbind(P1[!.(P_AFFECTED), c("g", "x", "y"), on = "g"], data.table::rbindlist(lapply(1:length(origin_c), function(O) {
                        return(P[igraph::shortest_paths(rst_u, origin_c[O], destination_c[O], output = "vpath",
                          algorithm = "dijkstra")$vpath[[1L]],][, g := P_AFFECTED[O]])
                      }), use.names = FALSE)[, c("g", "x", "y")])
                      data.table::setorder(P, g)
                      if(!ncoresg1) P <- terra::vect(as.matrix(P), type = "line", atts = P3, crs = r_crs)
                    }
                  } else {
                    P <- P1[, c("g", "x", "y")]
                  }
                  return(P)
                }
                if(ncoresg1 && par_lvl == "points_lists") {
                  v <- parallel::clusterMap(cl, p, p_affected, p1, p2, USE.NAMES = FALSE, .scheduling = "dynamic")
                  v <- mapply(function(P, P3) terra::vect(as.matrix(P), type = "line", atts = P3, crs = r_crs), v, p3, SIMPLIFY = FALSE, USE.NAMES = FALSE)
                } else if(ncoresg1 && par_lvl == "update_rst") {
                  v <- mapply(p, p_affected, p1, p2, SIMPLIFY = FALSE, USE.NAMES = FALSE)
                } else {
                  v <- mapply(p, p_affected, p1, p2, p3, SIMPLIFY = FALSE, USE.NAMES = FALSE)
                }
              } else {
                if(ncoresg1 && par_lvl == "update_rst") {
                  v <- lapply(p1, function(p) p[, c("g", "x", "y")])
                } else {
                  v <- mapply(function(P1, P3) terra::vect(as.matrix(P1[, c("g", "x", "y")]), type = "line", atts = P3, crs = r_crs), p1, p3,
                    SIMPLIFY = FALSE, USE.NAMES = FALSE)
                }
              }
              return(v)
            }
            if(update_rst_list) {
              if(ncoresg1 && par_lvl == "update_rst") {
                paths <- c(list(lapply(p1, function(p) p[, c("g", "x", "y")])), parallel::parLapply(cl, update_rst, paths))
                rm(p1, p2, update_rst)
                paths <- lapply(paths, function(p) mapply(function(P, P3) terra::vect(as.matrix(P), type = "line", atts = P3, crs = r_crs), p, p3,
                  SIMPLIFY = FALSE, USE.NAMES = FALSE))
              } else {
                paths <- c(list(mapply(function(P1, P3) terra::vect(as.matrix(P1[, c("g", "x", "y")]), type = "line", atts = P3, crs = r_crs), p1, p3,
                  SIMPLIFY = FALSE, USE.NAMES = FALSE)), lapply(update_rst, paths))
              }
            } else {
              paths <- list(mapply(function(P1, P3) terra::vect(as.matrix(P1[, c("g", "x", "y")]), type = "line", atts = P3, crs = r_crs), p1, p3,
                SIMPLIFY = FALSE, USE.NAMES = FALSE), paths(update_rst))
            }
          } else {
            paths <- function(u) {
              p_affected <- unique(p1[.(u), "g", nomatch = NULL, on = "cls"][["g"]])
              if(length(p_affected) > 0L) {
                origin_c <- p2[p_affected, "origin_c"]
                destination_c <- p2[p_affected, "destination_c"]
                if(any(origin_c %in% u)) {
                  report_points_ust(origin_c, u, TRUE, cl)
                }
                if(any(destination_c %in% u)) {
                  report_points_ust(destination_c, u, FALSE, cl)
                }
                v <- data.table::data.table(c_n_c = 1:NROW(crd))[-u,]
                origin_c <- v[.(origin_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                destination_c <- v[.(destination_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                v <- crd[-u, c("x", "y")]
                rst_u <- igraph::delete_vertices(rst, u)
                if(ncoresg1 && par_lvl == "points") {
                  v <- rbind(p1[!.(p_affected), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::clusterMap(cl, function(oc, dc, pa) {
                    return(v[igraph::shortest_paths(rst_u, oc, dc, output = "vpath", algorithm = "dijkstra")$vpath[[1L]],][, g := pa])
                  }, origin_c, destination_c, p_affected, USE.NAMES = FALSE), use.names = FALSE)[, c("g", "x", "y")])
                } else {
                  v <- rbind(p1[!.(p_affected), c("g", "x", "y"), on = "g"], data.table::rbindlist(lapply(1:length(origin_c), function(O) {
                    return(v[igraph::shortest_paths(rst_u, origin_c[O], destination_c[O], output = "vpath",
                      algorithm = "dijkstra")$vpath[[1L]],][, g := p_affected[O]])
                  }), use.names = FALSE)[, c("g", "x", "y")])
                }
                data.table::setorder(v, g)
              } else {
                v <- p1[, c("g", "x", "y")]
              }
              if(!(ncoresg1 && par_lvl == "update_rst")) v <- terra::vect(as.matrix(v), type = "line", atts = p3, crs = r_crs)
              return(v)
            }
            if(update_rst_list) {
              if(ncoresg1 && par_lvl == "update_rst") {
                paths <- c(list(p1[, c("g", "x", "y")]), parallel::parLapply(cl, update_rst, paths))
                rm(p1, p2, update_rst)
                paths <- lapply(paths, function(p) terra::vect(as.matrix(p), type = "line", atts = p3, crs = r_crs))
              } else {
                paths <- c(list(terra::vect(as.matrix(p1[, c("g", "x", "y")]), type = "line", atts = p3, crs = r_crs)), lapply(update_rst, paths))
              }
            } else {
              paths <- list(terra::vect(as.matrix(p1[, c("g", "x", "y")]), type = "line", atts = p3, crs = r_crs), paths(update_rst))
            }
          }
        # Distances output
        } else {
          if(p_list) {
            paths <- function(u) {
              p_affected <- lapply(p2, function(p) unique(p[.(u), "g", nomatch = NULL, on = "cls"][["g"]]))
              if(any(lengths(p_affected) > 0L)) {
                rst_u <- igraph::delete_vertices(rst, u)
                v <- data.table::data.table(c_n_c = 1:NROW(crd))[-u,]
                p <- function(P_AFFECTED, P1) {
                  if(length(P_AFFECTED) > 0L) {
                    origin_c <- P1[P_AFFECTED, "origin_c"][["origin_c"]]
                    destination_c <- P1[P_AFFECTED, "destination_c"][["destination_c"]]
                    if(any(origin_c %in% u)) {
                      report_points_ust(origin_c, u, TRUE, cl)
                    }
                    if(any(destination_c %in% u)) {
                      report_points_ust(destination_c, u, FALSE, cl)
                    }
                    origin_c <- v[.(origin_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                    destination_c <- v[.(destination_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                    P <- data.table::copy(P1[, c("origin", "destination", "distance")])
                    if(ncoresg1 && par_lvl == "points") {
                      P[P_AFFECTED, distance := parallel::clusterMap(cl, function(oc, dc) igraph::distances(rst_u, oc, dc, mode = "out",
                        algorithm = "dijkstra"), origin_c, destination_c, SIMPLIFY = TRUE, USE.NAMES = FALSE)]
                    } else {
                      P[P_AFFECTED, distance := vapply(1:length(origin_c), function(O) igraph::distances(rst_u, origin_c[O], destination_c[O],
                        mode = "out", algorithm = "dijkstra"), numeric(1L), USE.NAMES = FALSE)]
                    }
                  } else {
                    P <- P1[, c("origin", "destination", "distance")]
                  }
                  return(P)
                }
                if(ncoresg1 && par_lvl == "points_lists") {
                  v <- parallel::clusterMap(cl, p, p_affected, p1, USE.NAMES = FALSE, .scheduling = "dynamic")
                } else {
                  v <- mapply(p, p_affected, p1, SIMPLIFY = FALSE, USE.NAMES = FALSE)
                }
              } else {
                v <- lapply(p1, function(p) p[, c("origin", "destination", "distance")])
              }
              return(v)
            }
            if(update_rst_list) {
              if(ncoresg1 && par_lvl == "update_rst") {
                paths <- c(list(lapply(p1, function(p) p[, c("origin", "destination", "distance")])), parallel::parLapply(cl, update_rst, paths))
              } else {
                paths <- c(list(lapply(p1, function(p) p[, c("origin", "destination", "distance")])), lapply(update_rst, paths))
              }
            } else {
              paths <- list(lapply(p1, function(p) p[, c("origin", "destination", "distance")]), paths(update_rst))
            }
          } else {
            paths <- function(u) {
              p_affected <- unique(p2[.(u), "g", nomatch = NULL, on = "cls"][["g"]])
              if(length(p_affected) > 0L) {
                origin_c <- p1[p_affected, "origin_c"][["origin_c"]]
                destination_c <- p1[p_affected, "destination_c"][["destination_c"]]
                if(any(origin_c %in% u)) {
                  report_points_ust(origin_c, u, TRUE, cl)
                }
                if(any(destination_c %in% u)) {
                  report_points_ust(destination_c, u, FALSE, cl)
                }
                v <- data.table::data.table(c_n_c = 1:NROW(crd))[-u,]
                origin_c <- v[.(origin_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                destination_c <- v[.(destination_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                v <- data.table::copy(p1[, c("origin", "destination", "distance")])
                rst_u <- igraph::delete_vertices(rst, u)
                v[p_affected, distance := vapply(1:length(origin_c), function(O) igraph::distances(rst_u, origin_c[O], destination_c[O], mode = "out",
                  algorithm = "dijkstra"), numeric(1L), USE.NAMES = FALSE)]
              } else {
                v <- p1[, c("origin", "destination", "distance")]
              }
              return(v)
            }
            if(update_rst_list) {
              if(ncoresg1 && par_lvl == "update_rst") {
                paths <- c(list(p1[, c("origin", "destination", "distance")]), parallel::parLapply(cl, update_rst, paths))
              } else {
                paths <- c(list(p1[, c("origin", "destination", "distance")]), lapply(update_rst, paths))
              }
            } else {
              paths <- list(p1[, c("origin", "destination", "distance")], paths(update_rst))
            }
          }
        }
      # Not copy
      } else {
        if(ncoresg1 && par_lvl == "update_rst") {
          paths <- function(V) {
            if(V == 0L) {
              v <- compute_spaths1(rst, crd[, c("x", "y")], origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified,
                origin_list, dest_list, r_crs, output_lines, pairwise, FALSE)
            } else {
              if(update_rst_list) {
                v <- compute_spaths1(igraph::delete_vertices(rst, update_rst[[V]]), crd[-update_rst[[V]], c("x", "y")], origins, destinations, dest_specified,
                  origin_nms_specified, destination_nms_specified, origin_list, dest_list, r_crs, output_lines, pairwise, FALSE)
              } else {
                v <- compute_spaths1(igraph::delete_vertices(rst, update_rst), crd[-update_rst, c("x", "y")], origins, destinations, dest_specified,
                  origin_nms_specified, destination_nms_specified, origin_list, dest_list, r_crs, output_lines, pairwise, FALSE)
              }
            }
            return(v)
          }
          if(update_rst_list) {
            paths <- parallel::parLapply(cl, 0:length(update_rst), paths)
          } else {
            paths <- parallel::parLapply(cl, 0:1, paths)
          }
          if(output_lines) {
            if(p_list) {
              if(origin_list) {
                if(dest_specified) {
                  if(dest_list) {
                    if(pairwise) {
                      if(origin_nms_specified && destination_nms_specified) {
                        paths <- lapply(paths, function(P) mapply(function(p, o, d) terra::vect(as.matrix(p), type = "line",
                          atts = data.frame(origin = o[["nms"]], destination = d[["nms"]]), crs = r_crs), P, origins, destinations, SIMPLIFY = FALSE,
                          USE.NAMES = FALSE))
                      } else if(origin_nms_specified) {
                        rm(destinations)
                        paths <- lapply(paths, function(P) mapply(function(p, o) terra::vect(as.matrix(p), type = "line",
                          atts = data.frame(origin = o[["nms"]], destination = 1:NROW(p)), crs = r_crs), P, origins, SIMPLIFY = FALSE, USE.NAMES = FALSE))
                      } else if(destination_nms_specified) {
                        rm(origins)
                        paths <- lapply(paths, function(P) mapply(function(p, d) terra::vect(as.matrix(p), type = "line",
                          atts = data.frame(origin = 1:NROW(p), destination = d[["nms"]]), crs = r_crs), P, destinations, SIMPLIFY = FALSE,
                          USE.NAMES = FALSE))
                      } else {
                        rm(origins, destinations)
                        paths <- lapply(paths, function(P) lapply(P, function(p) terra::vect(as.matrix(p), type = "line",
                          atts = data.frame(origin = 1:NROW(p), destination = 1:NROW(p)), crs = r_crs)))
                      }
                    } else {
                      if(origin_nms_specified && destination_nms_specified) {
                        paths <- lapply(paths, function(P) mapply(function(p, o, d) terra::vect(as.matrix(p), type = "line",
                          atts = data.frame(origin = rep.int(o, rep.int(NROW(d), NROW(o))), destination = rep.int(d, NROW(o))), crs = r_crs), P, origins,
                          destinations, SIMPLIFY = FALSE, USE.NAMES = FALSE))
                      } else if(origin_nms_specified) {
                        destinations <- lengths(destinations)
                        paths <- lapply(paths, function(P) mapply(function(p, o, d) terra::vect(as.matrix(p), type = "line",
                          atts = data.frame(origin = rep.int(o, rep.int(d, NROW(o))), destination = rep.int(1:d, NROW(o))), crs = r_crs), P, origins,
                          destinations, SIMPLIFY = FALSE, USE.NAMES = FALSE))
                      } else if(destination_nms_specified) {
                        origins <- lengths(origins)
                        paths <- lapply(paths, function(P) mapply(function(p, o, d) terra::vect(as.matrix(p), type = "line",
                          atts = data.frame(origin = rep.int(1:o, rep.int(NROW(d), o)), destination = rep.int(d, o)), crs = r_crs), P, origins,
                          destinations, SIMPLIFY = FALSE, USE.NAMES = FALSE))
                      } else {
                        origins <- lengths(origins)
                        destinations <- lengths(destinations)
                        paths <- lapply(paths, function(P) mapply(function(p, o, d) terra::vect(as.matrix(p), type = "line",
                          atts = data.frame(origin = rep.int(1:o, rep.int(d, o)), destination = rep.int(1:d, o)), crs = r_crs), P, origins, destinations,
                          SIMPLIFY = FALSE, USE.NAMES = FALSE))
                      }
                    }
                  } else {
                    # The function demands pairwise to be TRUE in the origin_list = TRUE, dest_list = FALSE setting
                    if(origin_nms_specified && destination_nms_specified) {
                      destinations[, cls := NULL]
                      paths <- lapply(paths, function(P) mapply(function(p, o) terra::vect(as.matrix(p), type = "line",
                        atts = data.frame(origin = o[["nms"]], destination = destinations), crs = r_crs), P, origins, SIMPLIFY = FALSE,
                        USE.NAMES = FALSE))
                    } else if(origin_nms_specified) {
                      rm(destinations)
                      paths <- lapply(paths, function(P) mapply(function(p, o) terra::vect(as.matrix(p), type = "line",
                        atts = data.frame(origin = o[["nms"]], destination = 1:NROW(p)), crs = r_crs), P, origins, SIMPLIFY = FALSE, USE.NAMES = FALSE))
                    } else if(destination_nms_specified) {
                      rm(origins)
                      destinations <- destinations[["nms"]]
                      paths <- lapply(paths, function(P) lapply(P, function(p) terra::vect(as.matrix(p), type = "line",
                        atts = data.frame(origin = 1:NROW(p), destination = destinations), crs = r_crs)))
                    } else {
                      rm(origins)
                      destinations <- length(destinations)
                      paths <- lapply(paths, function(P) lapply(P, function(p) terra::vect(as.matrix(p), type = "line",
                        atts = data.frame(origin = 1:destinations, destination = 1:destinations), crs = r_crs)))
                    }
                  }
                } else {
                  if(origin_nms_specified) {
                    paths <- lapply(paths, function(P) mapply(function(p, o) terra::vect(as.matrix(p), type = "line",
                      atts = stats::setNames(as.data.frame(t(utils::combn(o[["nms"]], 2L))), c("origin", "destination")), crs = r_crs), P, origins,
                      SIMPLIFY = FALSE, USE.NAMES = FALSE))
                  } else {
                    origins <- lengths(origins)
                    paths <- lapply(paths, function(P) mapply(function(p, o) terra::vect(as.matrix(p), type = "line",
                      atts = stats::setNames(as.data.frame(t(utils::combn(1:o, 2L))), c("origin", "destination")), crs = r_crs), P, origins,
                      SIMPLIFY = FALSE, USE.NAMES = FALSE))
                  }
                }
              } else {
                # The function demands pairwise to be TRUE in the origin_list = FALSE, dest_list = TRUE setting
                if(origin_nms_specified && destination_nms_specified) {
                  origins[, cls := NULL]
                  paths <- lapply(paths, function(P) mapply(function(p, d) terra::vect(as.matrix(p), type = "line",
                    atts = data.frame(origin = origins, destination = d[["nms"]]), crs = r_crs), P, destinations, SIMPLIFY = FALSE, USE.NAMES = FALSE))
                } else if(origin_nms_specified) {
                  origins[, cls := NULL]
                  paths <- lapply(paths, function(P) lapply(P, function(p) terra::vect(as.matrix(p), type = "line", atts = data.frame(origin = origins,
                    destination = 1:NROW(p)), crs = r_crs)))
                } else if(destination_nms_specified) {
                  origins <- length(origins)
                  paths <- lapply(paths, function(P) mapply(function(p, d) terra::vect(as.matrix(p), type = "line", atts = data.frame(origin = 1:origins,
                    destination = d[["nms"]]), crs = r_crs), P, destinations, SIMPLIFY = FALSE, USE.NAMES = FALSE))
                } else {
                  rm(destinations)
                  origins <- length(origins)
                  paths <- lapply(paths, function(P) lapply(P, function(p) terra::vect(as.matrix(p), type = "line", atts = data.frame(origin = 1:origins,
                    destination = 1:origins), crs = r_crs)))
                }
              }
            } else {
              if(dest_specified) {
                if(pairwise) {
                  if(origin_nms_specified && destination_nms_specified) {
                    origins[, c("cls", "destination") := list(NULL, destinations[["nms"]])]
                    rm(destinations)
                    data.table::setnames(origins, "cls", "origin")
                    paths <- lapply(paths, function(P) terra::vect(as.matrix(P), type = "line", atts = origins, crs = r_crs))
                  } else if(origin_nms_specified) {
                    rm(destinations)
                    origins[, c("cls", "destination") := list(NULL, 1:.N)]
                    data.table::setnames(origins, "nms", "origin")
                    paths <- lapply(paths, function(P) terra::vect(as.matrix(P), type = "line", atts = origins, crs = r_crs))
                  } else if(destination_nms_specified) {
                    rm(origins)
                    destinations[, c("cls", "origin") := c(NULL, 1:.N)]
                    data.table::setnames(origins, "nms", "destination")
                    data.table::setcolorder(destinations, c("origin", "destination"))
                    paths <- lapply(paths, function(P) terra::vect(as.matrix(P), type = "line", atts = destinations, crs = r_crs))
                  } else {
                    rm(destinations)
                    origins <- 1:NROW(origins)
                    origins <- data.frame(origin = origins, destination = origins)
                    paths <- lapply(paths, function(P) terra::vect(as.matrix(P), type = "line", atts = origins, crs = r_crs))
                  }
                } else {
                  if(origin_nms_specified && destination_nms_specified) {
                    origins <- data.frame(origin = rep.int(origins[["nms"]], rep.int(NROW(destinations), NROW(origins))),
                      destination = rep.int(destinations[["nms"]], NROW(origins)))
                    rm(destinations)
                    paths <- lapply(paths, function(P) terra::vect(as.matrix(P), type = "line", atts = origins, crs = r_crs))
                  } else if(origin_nms_specified) {
                    destinations <- length(destinations)
                    origins <- data.frame(origin = rep.int(origins[["nms"]], rep.int(destinations, NROW(origins))), destination = rep.int(1:destinations,
                      NROW(origins)))
                    paths <- lapply(paths, function(P) terra::vect(as.matrix(P), type = "line", atts = origins, crs = r_crs))
                  } else if(destination_nms_specified) {
                    origins <- length(origins)
                    destinations <- data.frame(origin = rep.int(1:origins, rep.int(NROW(destinations), origins)), destination = rep.int(destinations[["nms"]],
                      origins))
                    paths <- lapply(paths, function(P) terra::vect(as.matrix(P), type = "line", atts = destinations, crs = r_crs))
                  } else {
                    origins <- length(origins)
                    destinations <- length(destinations)
                    origins <- data.frame(origin = rep.int(1:origins, rep.int(destinations, origins)), destination = rep.int(1:destinations, origins))
                    paths <- lapply(paths, function(P) terra::vect(as.matrix(P), type = "line", atts = origins, crs = r_crs))
                  }
                }
              } else {
                if(origin_nms_specified) {
                  origins <- stats::setNames(as.data.frame(t(utils::combn(origins[["nms"]], 2L))), c("origin", "destination"))
                  paths <- lapply(paths, function(P) terra::vect(as.matrix(P), type = "line", atts = origins, crs = r_crs))
                } else {
                  origins <- stats::setNames(as.data.frame(t(utils::combn(1:length(origins), 2L))), c("origin", "destination"))
                  paths <- lapply(paths, function(P) terra::vect(as.matrix(p), type = "line", atts = origins, crs = r_crs))
                }
              }
            }
          }
        # Serial or not parallel at the update_rst level
        } else {
          if(update_rst_list) {
            paths <- lapply(0:length(update_rst), function(V) {
              if(V == 0L) {
                v <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified,
                  origin_list, dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, cl)
              } else {
                v <- compute_spaths1(igraph::delete_vertices(rst, update_rst[[V]]), crd[-update_rst[[V]],], origins, destinations, dest_specified,
                  origin_nms_specified, destination_nms_specified, origin_list, dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, cl)
              }
              return(v)
            })
          } else {
            paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified,
              origin_list, dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, cl)
            rst <- igraph::delete_vertices(rst, update_rst)
            crd <- crd[-update_rst,]
            paths <- rbind(paths, compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified,
              origin_list, dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, cl))
          }
        }
      }
    # update_rst does not mask any non-NA pixels
    } else {
      warning("The update_rst geometries do not mask any non-NA pixels. Thus, the results are the same as those from the unmodified rst.")
      if(update_rst_list) {
        update_rst <- length(update_rst) + 1L
      } else {
        update_rst <- 2L
      }
      if(!ncoresg1 || (par_lvl == "update_rst" && copy && rst_par_lvl == "none")) {
        paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list,
          dest_list, r_crs, output_lines, pairwise, FALSE)
      } else {
        if(par_lvl == "update_rst" && copy) par_lvl <- rst_par_lvl
        paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list,
          dest_list, r_crs, output_lines, pairwise, TRUE, ncores, par_lvl, cl)
      }
      paths <- replicate(update_rst, paths, FALSE)
    }
  # update_rst not specified
  } else {
    paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list, dest_list, r_crs,
      output == "lines", pairwise, ncoresg1, ncores, par_lvl, cl)
  }
  if(ncoresg1) parallel::stopCluster(cl)
  
  return(paths)
}

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
#' @param round_dist Logical specifying whether to round edge weights, i.e. the transition cost or distance between neighboring cells, to integers. It 
#' defaults to \code{FALSE}. Setting it to \code{TRUE} reduces the function's RAM requirements slightly.
#' @param ncores An integer specifying the number of CPU cores to use. It defaults to the number of cores installed on the machine. A value of 1 
#' accordingly induces a single-threaded process.
#' @param par_lvl \code{"points"} (default), \code{"points_lists"}, or \code{"update_rst"}, indicating the level at which to parallelize when 
#' \code{ncores > 1}. \code{"points"} parallelizes over the origin (and destination) point combinations. \code{"points_lists"} parallelizes over the list 
#' elements of \code{origins} (and \code{destinations}), if these arguments are lists. \code{"update_rst"} parallelizes over the list of graphs specified 
#' by \code{rst} and \code{update_rst}.
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
#' The vignette provides further details on the function.
#' 
#' @return In the basic cases, i.e. when neither \code{origins} nor \code{destinations} are lists and \code{update_rst} is not specified, the function 
#' returns a SpatVector lines object with \code{output = "lines"} and a data table with \code{output = "distances"}. If \code{origins} or 
#' \code{destinations} are lists or \code{update_rst} is specified, it returns a list or nested list of SpatVectors or data tables. In the nested case, 
#' the outer list refers to the \code{update_rst} level and the inner list to the \code{origins} and \code{destinations} level. Consult the vignette for 
#' further details.
#' 
#' The distances are the total transition costs along the shortest paths. If \code{lonlat} is \code{TRUE} or if \code{dist_comp = "terra"}, these 
#' distances are measured in the meters. Otherwise, they use the units of the CRS. If you pass a function to \code{tr_fun}, the total transition cost is 
#' measured in the units of \code{tr_fun}'s results.
#' 
#' @seealso \link{spaths_general}
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
  round_dist = FALSE, ncores = NULL, par_lvl = c("points", "points_lists", "update_rst"), cluster = NULL, paths_ncores = NULL, verbose = FALSE) {
  
  if(length(verbose) != 1L || !is.logical(verbose) || is.na(verbose)) stop("verbose must be logical and of length one")
  if(verbose) message("Checking arguments")
  
  # Convert RasterLayer to SpatRaster
  if(all(class(rst) != "SpatRaster")) rst <- terra::rast(rst)
  
  # Check CRS
  r_crs <- terra::crs(rst)
  if(r_crs == "") stop("rst without specified CRS")
  lonlat <- terra::is.lonlat(rst)
  
  # Check origin_names and destination_names
  origin_nms_null <- is.null(origin_names)
  if(!origin_nms_null && (!(is.character(origin_names) && length(origin_names) == 1L) || is.na(origin_names))) {
    stop("origin_names must be NULL or a character object of length one")
  }
  destination_nms_null <- is.null(destination_names)
  if(!destination_nms_null && (!(is.character(destination_names) && length(destination_names) == 1L) || is.na(destination_names))) {
    stop("destination_names must be NULL or a character object of length one")
  }
  
  # Convert origins
  origin_list <- is.list(origins)
  if(origin_list) {
    origins <- lapply(origins, convert_points, rst, r_crs, origin_names, origin_nms_null)
  } else {
    origins <- convert_points(origins, rst, r_crs, origin_names, origin_nms_null)
  }
  
  # Convert destinations
  dest_specified <- !is.null(destinations)
  if(dest_specified) {
    dest_list <- is.list(destinations)
    if(dest_list) {
      destinations <- lapply(destinations, convert_points, rst, r_crs, destination_names, destination_nms_null, FALSE)
    } else {
      destinations <- convert_points(destinations, rst, r_crs, destination_names, destination_nms_null, FALSE)
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
    if(is.list(update_rst)) {
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
  rst <- data.table::as.data.table(terra::adjacent(rst, crd[["c_n"]], contiguity, TRUE)) # Obtain adjacency matrix
  crd[, c_n_c := 1:.N]
  
  origins <- update_points(origins, origin_list, crd, origin_nms_null)
  if(dest_specified) destinations <- update_points(destinations, dest_list, crd, destination_nms_null)
  
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
    if(!is.vector(tr_fun_args)) stop("tr_fun must return a vector")
    if(any(tr_fun_args < 0)) stop("tr_fun must not return negative values")
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
    if(!upd_rst_specified && par_lvl == "update_rst") stop('If par_lvl is "update_rst", update_rst must not be NULL')
    if(!origin_list && dest_specified && !dest_list && par_lvl == "points_lists") {
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
    if(is.list(update_rst)) {
      if(ncoresg1 && par_lvl == "update_rst") {
        update_rst <- lapply(update_rst, function(V) crd[.(terra::extract(rst_upd, V, cells = TRUE, ID = FALSE, touches = touches)$cell), nomatch = NULL,
          which = TRUE, on = "c_n"])
        paths <- function(V) {
          if(V == 0L) {
            v <- compute_spaths1(rst, crd[, c("x", "y")], origins, destinations, dest_specified, origin_nms_null, destination_nms_null, origin_list,
              dest_list, r_crs, output_lines, pairwise, FALSE)
          } else {
            v <- compute_spaths1(igraph::delete_vertices(rst, update_rst[[V]]), crd[-update_rst[[V]], c("x", "y")], origins, destinations, dest_specified,
              origin_nms_null, destination_nms_null, origin_list, dest_list, r_crs, output_lines, pairwise, FALSE)
          }
          return(v)
        }
        paths <- parallel::parLapply(cl, 0:length(update_rst), paths)
        if(output_lines) {
          if(origin_list || (dest_specified && dest_list)) {
            paths <- lapply(paths, function(P) lapply(P, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]], crs = r_crs)))
          } else {
            paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]], crs = r_crs))
          }
        }
      } else {
        paths <- lapply(0:length(update_rst), function(V) {
          if(V == 0L) {
            v <- compute_spaths1(rst, crd[, c("x", "y")], origins, destinations, dest_specified, origin_nms_null, destination_nms_null, origin_list,
              dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, cl)
          } else {
            v <- crd[.(terra::extract(rst_upd, update_rst[[V]], cells = TRUE, ID = FALSE, touches = touches)$cell), nomatch = NULL, which = TRUE, on = "c_n"]
            v <- compute_spaths1(igraph::delete_vertices(rst, v), crd[-v, c("x", "y")], origins, destinations, dest_specified, origin_nms_null,
              destination_nms_null, origin_list, dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, cl)
          }
          return(v)
        })
      }
    } else {
      rst_upd <- crd[.(terra::extract(rst_upd, update_rst, cells = TRUE, ID = FALSE, touches = touches)$cell), nomatch = NULL, which = TRUE, on = "c_n"]
      rm(update_rst)
      crd[, c_n := NULL]
      if(ncoresg1 && par_lvl == "update_rst") {
        paths <- function(V) {
          if(V == 1L) {
            v <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_null, destination_nms_null, origin_list, dest_list, r_crs,
              output_lines, pairwise, FALSE)
          } else {
            v <- compute_spaths1(igraph::delete_vertices(rst, rst_upd), crd[-rst_upd,], origins, destinations, dest_specified, origin_nms_null,
              destination_nms_null, origin_list, dest_list, r_crs, output_lines, pairwise, FALSE)
          }
          return(v)
        }
        paths <- parallel::parLapply(cl, 1:2, paths)
        if(output_lines) {
          if(origin_list) {
            paths <- lapply(paths, function(P) lapply(P, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]], crs = r_crs)))
          } else {
            paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]], crs = r_crs))
          }
        }
      } else {
        paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_null, destination_nms_null, origin_list, dest_list, r_crs,
          output_lines, pairwise, ncoresg1, ncores, par_lvl, cl)
        rst <- igraph::delete_vertices(rst, rst_upd)
        crd <- crd[-rst_upd,]
        rm(rst_upd)
        paths <- list(paths, compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_null, destination_nms_null, origin_list, dest_list,
          r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, cl))
      }
    }
  } else {
    paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_null, destination_nms_null, origin_list, dest_list, r_crs,
      output == "lines", pairwise, ncoresg1, ncores, par_lvl, cl)
  }
  if(ncoresg1) parallel::stopCluster(cl)
  
  return(paths)
}

# Convert origins and destinations
convert_points <- function(v, rst, r_crs, nms, nms_null, o = TRUE) {
  if(all(class(v) != "SpatVector")) v <- terra::vect(v)
  if(!terra::is.points(v)) v <- terra::centroids(v)
  if(terra::crs(v) != r_crs) v <- terra::project(v, r_crs)
  if(!nms_null) {
    if(!(nms %chin% names(v))) {
      stop(ifelse(o, "origin", "destination"), "_names must either be NULL or the name of a column in ", ifelse(o, "origins", "destinations"))
    }
    nms <- unlist(terra::values(v[, nms]), use.names = FALSE)
  }
  v <- terra::extract(rst, v, cells = TRUE, ID = FALSE)
  if(any(!stats::complete.cases(v[, 1:terra::nlyr(rst)]))) {
    v <- which(!stats::complete.cases(v[, 1:terra::nlyr(rst)]))
    v_length <- length(v)
    if(v_length > 1L) {
      if(v_length > 2L) {
        v <- paste0(paste0(v[1:(v_length - 1L)], collapse = ", "), ", and ", v[v_length])
      } else {
        v <- paste0(v, collapse = " and ")
      }
      v <- paste0(ifelse(o, "Origins", "Destinations"), " ", v, " located on NA cells")
    } else {
      v <- paste0(ifelse(o, "Origin", "Destination"), " ", v, " located on NA cell")
    }
    stop(v)
  }
  if(nms_null) {
    v <- v$cell
  } else {
    v <- data.table::data.table(cls = v$cell, nms = nms)
  }
  return(v)
}

# Update points' cell numbers
update_points <- function(v, v_list, crd, nms_null) {
  if(v_list) {
    if(nms_null) {
      v <- lapply(v, function(V) {
        return(crd[.(V), "c_n_c", on = "c_n"][["c_n_c"]])
      })
    } else {
      v <- lapply(v, function(V) {
        return(data.table::setnames(crd[V, c("c_n_c", "nms"), on = "c_n==cls"], "c_n_c", "cls"))
      })
    }
  } else {
    if(nms_null) {
      v <- crd[.(v), "c_n_c", on = "c_n"][["c_n_c"]]
    } else {
      v <- data.table::setnames(crd[v, c("c_n_c", "nms"), on = "c_n==cls"], "c_n_c", "cls")
    }
  }
  return(v)
}

# Compute distances between neighboring cells
# Capitalized object names differing from the ones in the spaths_earth function indicate that these objects might differ from the spaths_earth
# counterparts, e.g. by being subsets
compute_dists <- function(rst, CRD, dist_comp_terra, round_dist, contiguity, yr, xr, nr, ym, lonlat, ncoresg1, ncores, cl) {
  if(dist_comp_terra) {
    if(ncoresg1) {
      d <- function(RST) {
        return(terra::distance(as.matrix(CRD[RST[["from"]],]), as.matrix(CRD[RST[["to"]],]), lonlat = lonlat, pairwise = TRUE))
      }
      d <- parallel::parLapply(cl, split(rst, cut(1:NROW(rst), ncores, labels = FALSE)), d)
    } else {
      d <- terra::distance(as.matrix(CRD[rst[["from"]],]), as.matrix(CRD[rst[["to"]],]), lonlat = lonlat, pairwise = TRUE)
    }
  } else {
    if(round_dist) {
      if(contiguity == "queen") {
        d <- dists_queen_i(CRD[rst[["from"]], "y"][["y"]], CRD[rst[["from"]], "x"][["x"]], CRD[rst[["to"]], "y"][["y"]],
          CRD[rst[["to"]], "x"][["x"]], yr, xr, nr, ym, lonlat, ncores)
      } else {
        d <- dists_rook_i(CRD[rst[["from"]], "y"][["y"]], CRD[rst[["from"]], "x"][["x"]], CRD[rst[["to"]], "y"][["y"]],
          CRD[rst[["to"]], "x"][["x"]], yr, xr, nr, ym, lonlat, ncores)
      }
    } else {
      if(contiguity == "queen") {
        d <- dists_queen_d(CRD[rst[["from"]], "y"][["y"]], CRD[rst[["from"]], "x"][["x"]], CRD[rst[["to"]], "y"][["y"]],
          CRD[rst[["to"]], "x"][["x"]], yr, xr, nr, ym, lonlat, ncores)
      } else {
        d <- dists_rook_d(CRD[rst[["from"]], "y"][["y"]], CRD[rst[["from"]], "x"][["x"]], CRD[rst[["to"]], "y"][["y"]],
          CRD[rst[["to"]], "x"][["x"]], yr, xr, nr, ym, lonlat, ncores)
      }
    }
  } 
  return(d)
}

# Function calling compute_spaths2
compute_spaths1 <- function(rst, crd, origins, destinations, dest_specified, origin_nms_null, destination_nms_null, origin_list, dest_list, r_crs,
  output_lines, pairwise, NCORESG1, ncores = NULL, par_lvl = NULL, cl = NULL) {
  if(origin_list) {
    if(dest_specified) {
      if(dest_list) {
        if(NCORESG1 && par_lvl == "points_lists") {
          paths <- function(O) {
            return(compute_spaths2(origins[[O]], rst, crd, dest_specified, origin_nms_null, r_crs, output_lines, pairwise, FALSE, NULL, NULL, TRUE,
              destinations[[O]], destination_nms_null))
          }
          paths <- parallel::parLapplyLB(cl, 1:length(origins), paths)
          if(output_lines) paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]], crs = r_crs))
        } else {
          paths <- mapply(compute_spaths2, ORIGINS = origins, DESTINATIONS = destinations, MoreArgs = list(rst = rst, crd = crd,
            dest_specified = dest_specified, origin_nms_null = origin_nms_null, r_crs = r_crs, output_lines = output_lines, pairwise = pairwise,
            NCORESG1 = NCORESG1, ncores = ncores, cl = cl, nvect = is.null(par_lvl), destination_nms_null = destination_nms_null), SIMPLIFY = FALSE,
            USE.NAMES = FALSE)
        }
      } else {
        if(NCORESG1 && par_lvl == "points_lists") {
          paths <- function(O) {
            return(compute_spaths2(O, rst, crd, dest_specified, origin_nms_null, r_crs, output_lines, pairwise, FALSE, NULL, NULL, TRUE, destinations,
              destination_nms_null))
          }
          paths <- parallel::parLapplyLB(cl, origins, paths)
          if(output_lines) paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]], crs = r_crs))
        } else {
          paths <- lapply(origins, compute_spaths2, rst, crd, dest_specified, origin_nms_null, r_crs, output_lines, pairwise, NCORESG1, ncores, cl,
            is.null(par_lvl), destinations, destination_nms_null)
        }
      }
    } else {
      paths <- lapply(origins, compute_spaths2, rst, crd, dest_specified, origin_nms_null, r_crs, output_lines, pairwise, NCORESG1, ncores, cl,
        is.null(par_lvl))
    }
  } else {
    if(dest_specified) {
      if(dest_list) {
        if(NCORESG1 && par_lvl == "points_lists") {
          paths <- function(D) {
            return(compute_spaths2(origins, rst, crd, dest_specified, origin_nms_null, r_crs, output_lines, pairwise, FALSE, NULL, NULL, TRUE,
              destinations[[D]], destination_nms_null))
          }
          paths <- parallel::parLapplyLB(cl, 1:length(destinations), paths)
          if(output_lines) paths <- lapply(paths, function(D) terra::vect(D[[1L]], type = "line", atts = D[[2L]], crs = r_crs))
        } else {
          paths <- lapply(destinations, function(d) compute_spaths2(origins, rst, crd, dest_specified, origin_nms_null, r_crs, output_lines, pairwise,
            NCORESG1, ncores, cl, is.null(par_lvl), d, destination_nms_null))
        }
      } else {
        paths <- compute_spaths2(origins, rst, crd, dest_specified, origin_nms_null, r_crs, output_lines, pairwise, NCORESG1, ncores, cl, is.null(par_lvl),
          destinations, destination_nms_null)
      }
    } else {
      paths <- compute_spaths2(origins, rst, crd, dest_specified, origin_nms_null, r_crs, output_lines, pairwise, NCORESG1, ncores, cl, is.null(par_lvl))
    }
  }
  return(paths)
}

# Compute shortest paths
compute_spaths2 <- function(ORIGINS, rst, crd, dest_specified, origin_nms_null, r_crs, output_lines, pairwise, NCORESG1, ncores, cl, nvect,
  DESTINATIONS = NULL, destination_nms_null = TRUE) {
  os_length <- NROW(ORIGINS)
  if(origin_nms_null) {
    on <- 1:os_length
  } else {
    on <- ORIGINS[["nms"]]
    ORIGINS <- ORIGINS[["cls"]]
  }
  # Shortest paths when destinations are specified
  if(dest_specified) {
    ds_length <- NROW(DESTINATIONS)
    if(destination_nms_null) {
      dn <- 1:ds_length
    } else {
      dn <- DESTINATIONS[["nms"]]
      DESTINATIONS <- DESTINATIONS[["cls"]]
    }
    # Lines output
    if(output_lines) {
      if(pairwise) {
        if(NCORESG1) {
          p <- function(O) {
            return(data.table::rbindlist(lapply(O, function(o) data.table::data.table(cls = igraph::shortest_paths(rst, ORIGINS[o], DESTINATIONS[o],
              output = "vpath", algorithm = "dijkstra")$vpath[[1L]], g = o))))
          }
          p <- terra::vect(as.matrix(data.table::rbindlist(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)),
            p))[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line", atts = data.frame(origin = on, destination = dn), crs = r_crs)
        } else if(nvect) {
          p <- list(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
            return(data.table::rbindlist(lapply(igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS[O], output = "vpath", algorithm = "dijkstra")$vpath,
              function(D) crd[D,][, g := O])))
          }))[, c("g", "x", "y")]), data.frame(origin = on, destination = dn))
        } else {
          p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
            return(data.table::rbindlist(lapply(igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS[O], output = "vpath", algorithm = "dijkstra")$vpath,
              function(D) crd[D,][, g := O])))
          }))[, c("g", "x", "y")]), type = "line", atts = data.frame(origin = on, destination = dn), crs = r_crs)
        }
      } else {
        if(NCORESG1) {
          p <- function(O) {
            return(data.table::rbindlist(lapply(O, function(o) {
              s <- igraph::shortest_paths(rst, ORIGINS[o], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
              o1 <- (o - 1L) * ds_length
              return(data.table::rbindlist(lapply(1:ds_length, function(D) data.table::data.table(cls = s[[D]], g = o1 + D))))
            })))
          }
          p <- terra::vect(as.matrix(data.table::rbindlist(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)),
            p))[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line", atts = data.frame(origin = rep.int(on, rep.int(ds_length, os_length)),
            destination = rep.int(dn, os_length)), crs = r_crs)
        } else if(nvect) {
          p <- list(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
            s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
            o1 <- (O - 1L) * ds_length
            return(data.table::rbindlist(lapply(1:ds_length, function(D) crd[s[[D]],][, g := o1 + D])))
          }))[, c("g", "x", "y")]), data.frame(origin = rep.int(on, rep.int(ds_length, os_length)), destination = rep.int(dn, os_length)))
        } else {
          p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
            s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
            o1 <- (O - 1L) * ds_length
            return(data.table::rbindlist(lapply(1:ds_length, function(D) crd[s[[D]],][, g := o1 + D])))
          }))[, c("g", "x", "y")]), type = "line", atts = data.frame(origin = rep.int(on, rep.int(ds_length, os_length)), destination = rep.int(dn,
            os_length)), crs = r_crs)
        }
      }
    # Distances output
    } else {
      if(pairwise) {
        if(NCORESG1) {
          p <- function(O) {
            return(vapply(O, function(o) igraph::distances(rst, ORIGINS[o], DESTINATIONS[o], mode = "out", algorithm = "dijkstra"), numeric(1L),
              USE.NAMES = FALSE))
          }
          p <- data.table::data.table(origin = on, destination = dn, distance = do.call(c, parallel::parLapply(cl, split(1:os_length, cut(1:os_length,
            ncores, labels = FALSE)), p)))
        } else {
          p <- data.table::data.table(origin = on, destination = dn, distance = vapply(1:os_length, function(O) igraph::distances(rst, ORIGINS[O],
            DESTINATIONS[O], mode = "out", algorithm = "dijkstra"), numeric(1L), USE.NAMES = FALSE))
        }
      } else {
        if(NCORESG1) {
          if(os_length >= ncores || os_length > ds_length) {
            p <- function(O) {
              return(data.table::as.data.table(igraph::distances(rst, O, DESTINATIONS, mode = "out", algorithm = "dijkstra"), na.rm = FALSE))
            }
            p <- data.table::rbindlist(parallel::parLapply(cl, split(ORIGINS, cut(seq_along(ORIGINS), ncores, labels = FALSE)), p))
          } else {
            p <- function(D) {
              return(igraph::distances(rst, ORIGINS, D, mode = "out", algorithm = "dijkstra"))
            }
            p <- data.table::as.data.table(do.call(cbind, parallel::parLapply(cl, split(DESTINATIONS, cut(seq_along(DESTINATIONS), ncores, labels = FALSE)),
              p)), na.rm = FALSE)
          }
        } else {
          p <- data.table::as.data.table(igraph::distances(rst, ORIGINS, DESTINATIONS, mode = "out", algorithm = "dijkstra"), na.rm = FALSE)
        }
        if(is.numeric(dn)) {
          data.table::setnames(p, as.character(dn))
          p[, origin := on]
          if(is.integer(dn)) {
            p <- data.table::melt(p, "origin", variable.name = "destination", value.name = "distance",
              variable.factor = FALSE)[, destination := as.integer(destination)]
          } else {
            p <- data.table::melt(p, "origin", variable.name = "destination", value.name = "distance",
              variable.factor = FALSE)[, destination := as.numeric(destination)]
          }
        } else {
          data.table::setnames(p, dn)
          p[, origin := on]
          p <- data.table::melt(p, "origin", variable.name = "destination", value.name = "distance", variable.factor = FALSE)
        }
      }
    }
  # Shortest paths when destinations are not specified
  } else {
    # Lines output
    if(output_lines) {
      o1 <- os_length - 1L
      if(NCORESG1) {
        p <- function(O) {
          return(data.table::rbindlist(lapply(split(O, by = "V1"), function(o) {
            s <- o[1L, "V1"][["V1"]]
            npO <- os_length - s
            i <- sum(o1:npO) - npO
            s <- igraph::shortest_paths(rst, ORIGINS[s], ORIGINS[o[["V2"]]], output = "vpath", algorithm = "dijkstra")$vpath
            return(data.table::rbindlist(lapply(1:NROW(o), function(D) data.table::data.table(cls = s[[D]], g = i + o[D, "V2"][["V2"]]))))
          })))
        }
        p <- terra::vect(as.matrix(data.table::rbindlist(parallel::parLapply(cl, split(data.table::as.data.table(t(utils::combn(1:os_length, 2L))),
         cut(1:(os_length * o1 / 2L), ncores, labels = FALSE)), p))[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line",
           atts = stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")), crs = r_crs)
      } else if(nvect) {
        p <- list(as.matrix(data.table::rbindlist(lapply(1:o1, function(O) {
          npO <- os_length - O
          i <- sum(o1:npO) - npO
          s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
          return(data.table::rbindlist(lapply(1:npO, function(D) crd[s[[D]],][, g := i + D])))
        }))[, c("g", "x", "y")]), stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")))
      } else {
        p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:o1, function(O) {
          npO <- os_length - O
          i <- sum(o1:npO) - npO
          s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
          return(data.table::rbindlist(lapply(1:npO, function(D) crd[s[[D]],][, g := i + D])))
        }))[, c("g", "x", "y")]), type = "line", atts = stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")), crs = r_crs)
      }
    # Distances output
    } else {
      on_n <- is.numeric(on)
      if(on_n) on_i <- is.integer(on)
      if(NCORESG1) {
        p <- function(O) {
          return(data.table::rbindlist(lapply(split(O, by = "V1"), function(o) {
            g <- o[1L, "V1"][["V1"]]
            s <- data.table::as.data.table(igraph::distances(rst, ORIGINS[g], ORIGINS[o[["V2"]]], mode = "out", algorithm = "dijkstra"), na.rm = FALSE)
            if(on_n) {
              data.table::setnames(s, as.character(on[o[["V2"]]]))
              s[, origin := on[g]]
              if(on_i) {
                s <- data.table::melt(s, "origin", variable.name = "destination", value.name = "distance",
                  variable.factor = FALSE)[, destination := as.integer(destination)]
              } else {
                s <- data.table::melt(s, "origin", variable.name = "destination", value.name = "distance",
                  variable.factor = FALSE)[, destination := as.numeric(destination)]
              }
            } else {
              data.table::setnames(s, on[o[["V2"]]])
              s[, origin := on[g]]
              s <- data.table::melt(s, "origin", variable.name = "destination", value.name = "distance", variable.factor = FALSE)
            }
            return(s)
          })))
        }
        p <- data.table::rbindlist(parallel::parLapply(cl, split(data.table::as.data.table(t(utils::combn(1:os_length, 2L))),
          cut(1:(os_length * (os_length - 1L) / 2L), ncores, labels = FALSE)), p))
      } else {
        p <- data.table::rbindlist(lapply(1:(os_length - 1L), function(O) {
          s <- data.table::as.data.table(igraph::distances(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], mode = "out", algorithm = "dijkstra"),
            na.rm = FALSE)
          if(on_n) {
            data.table::setnames(s, as.character(on[(O + 1L):os_length]))
            s[, origin := on[O]]
            if(on_i) {
              s <- data.table::melt(s, "origin", variable.name = "destination", value.name = "distance",
                variable.factor = FALSE)[, destination := as.integer(destination)]
            } else {
              s <- data.table::melt(s, "origin", variable.name = "destination", value.name = "distance",
                variable.factor = FALSE)[, destination := as.numeric(destination)]
            }
          } else {
            data.table::setnames(s, on[(O + 1L):os_length])
            s[, origin := on[O]]
            s <- data.table::melt(s, "origin", variable.name = "destination", value.name = "distance", variable.factor = FALSE)
          }
          return(s)
        }))
      }
    }
  }
  return(p)
}

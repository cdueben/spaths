#' Shortest Paths between Geographic Locations in General
#' 
#' The function computes the shortest paths between locations in general using Dijstra's algorithm. While \code{spaths_earth} derives shortest paths between 
#' locations on Earth, this function deduces them on any sphere or plane.
#' 
#' @param rst Matrix or list of matrices denoting the grid through via the points are connected. Pixels with non-NA values in all matrices mark the cells 
#' through which the algorithm may pass. A list of matrices is handled as a multi-layered grid, allowing for a cell to have multiple values. Thus, only 
#' use a list instead of a single matrix, if you are utilizing the cell values, i.e. passing a function with the \code{v1} or \code{v2} parameters to 
#' \code{tr_fun}.
#' @param xres Numeric specifying \code{rst}'s resolution along the x dimension, as in \code{terra::xres}.
#' @param yres Numeric specifying \code{rst}'s resolution along the y dimension, as in \code{terra::yres}.
#' @param xmin Numeric specifying \code{rst}'s xmin, as in \code{terra::xmin}.
#' @param ymin Numeric specifying \code{rst}'s ymin, as in \code{terra::ymin}.
#' @param origins A single matrix or data frame, or a list of them. The first column states the x and the second column the y coordinates, irrespective 
#' of what the column names are.
#' @param destinations A single matrix or data frame, a list of them. The first column states the x and the second column the y coordinates, irrespective 
#' of what the column names are. If not specified, the function to computes shortest paths between the \code{origins} points. Passing a list of data 
#' frames or matrices is a way not to relate all origins to all destinations. Details are outlined below.
#' @param lonlat Logical specifying whether the data is in lonlat (\code{TRUE}), i.e. spherical, or in planar (\code{FALSE}) format. It defaults to 
#' \code{TRUE}.
#' @param radius The radius of the sphere, if \code{lonlat} is \code{TRUE}. It is ignored, if \code{lonlat} is \code{FALSE} of if the function passed to 
#' \code{tr_fun} does not use the parameter \code{d}.
#' @param output \code{"lines"} (default) or \code{"distances"}. \code{"lines"} returns the shortest paths as SpatVector lines. \code{"distances"} lists 
#' the total transition costs along the shortest paths, which is by default the distance between origin and destination. If \code{lonlat} is \code{TRUE}, 
#' these distances are measured in the same units as \code{radius}. If \code{lonlat} is \code{FALSE}, the units are the same as those of \code{yres}, 
#' \code{yres}, \code{xmin}, \code{ymin}, \code{origins}, and \code{destinations}. If you pass a function to \code{tr_fun}, the total transition cost is 
#' measured in the units of \code{tr_fun}'s results. \code{"distances"} is faster and requires less RAM than \code{"lines"}.
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
#' @param tr_fun The transition function based on which to compute edge weights, i.e. the travel cost between grid cells. Defaults to the geographic 
#' distance between pixel centroids. Permitted function parameter names are \code{d} (distance between the pixel centroids), \code{x1} (x 
#' coordinate or longitude of the first cell), \code{x2} (x coordinate or longitude of the second cell), \code{y1} (y coordinate or latitude of the first 
#' cell), \code{y2} (y coordinate or latitude of the second cell), \code{v1} (\code{rst}'s value from the first cell), \code{v2} (\code{rst}'s 
#' value from the second cell), \code{nc} (number of CPU cores), and \code{cl} (cluster object returned by \code{parallel::makePSOCKcluster(ncores)} or 
#' \code{parallel::makeForkCluster(ncores)}). If \code{lonlat} is \code{TRUE}, \code{d} is measured in the same units as \code{radius}. If \code{lonlat} 
#' is \code{FALSE}, the units are the same as those of \code{yres}, \code{yres}, \code{xmin}, \code{ymin}, \code{origins}, and \code{destinations}. 
#' \code{nc} is meant to be used in C++ functions. A parallel backend at the R level pointed to by \code{cl} is already registered before \code{tr_fun} 
#' is called, in case \code{ncores > 1}. If \code{rst} is a single matrix, the values are passed to \code{v1} and \code{v2} as vectors, otherwise they 
#' are passed as a data table where the first column refers to the first matrix, the second column to the second matrix etc. Note that data tables are 
#' also data frames.
#' @param v_matrix Logical specifying whether to pass values to \code{v1} and \code{v2} in \code{tr_fun} as matrices (\code{TRUE}) instead of data tables 
#' in the list case and vectors in the single matrix case (\code{FALSE}). It defaults to \code{FALSE}. Setting it to \code{TRUE} might e.g. be useful when 
#' defining \code{tr_fun} as a C++ Armadillo function.
#' @param tr_directed Logical specifying whether \code{tr_fun} creates a directed graph. In a directed graph, transition costs can be asymmetric. 
#' Traveling from cell A to cell B may imply a different cost than traveling from B to A. It defaults to \code{TRUE} and only has an effect when 
#' \code{tr_fun} is not \code{NULL}. The default without \code{tr_fun} constructs an undirected graph.
#' @param round_dist Logical specifying whether to round edge weights, i.e. the transition cost between neighboring cells, to integers. It defaults to 
#' \code{FALSE}. Setting it to \code{TRUE} reduces the function's RAM requirements slightly.
#' @param ncores An integer specifying the number of CPU cores to use. It defaults to the number of cores installed on the machine. A value of 1 
#' accordingly induces a single-threaded process.
#' @param par_lvl \code{"points"} (default) or \code{"points_lists"}, indicating the level at which to parallelize when \code{ncores > 1}. \code{"points"} 
#' parallelizes over the origin (and destination) point combinations. \code{"points_lists"} parallelizes over the list elements of \code{origins} (and 
#' \code{destinations}), if these arguments are lists.
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
#' @details This function computes shortest paths between points on a sphere or a plane, taking custom transition costs into account. The sphere could be 
#' the moon or an orange. The plane could be the projected surface of a planet or a table top. The function computes shortest paths, irrespective of 
#' whether it represents a rover navigating the Marsian topography or an ant traveling between a coffee mug and a laptop on your desk.
#' 
#' The algorithm links \code{origins} and \code{destinations} by passing through the centroids of \code{rst}'s non-NA cells.
#' 
#' The \href{../doc/spaths_introduction.html}{vignette} provides further details on the function.
#' 
#' @return In the basic cases, i.e. when neither \code{origins} nor \code{destinations} are lists, the function returns a SpatVector lines object with 
#' \code{output = "lines"} and a data table with \code{output = "distances"}. If \code{origins} or \code{destinations} are lists, it returns a list of 
#' SpatVectors or data tables. Consult the \href{../doc/spaths_introduction.html}{vignette} for further details.
#' 
#' The distances are the total transition costs along the shortest paths. If \code{lonlat} is \code{TRUE}, these distances are measured in the same units 
#' as \code{radius}. If \code{lonlat} is \code{FALSE}, the units are the same as those of \code{yres}, \code{yres}, \code{xmin}, \code{ymin}, 
#' \code{origins}, and \code{destinations}. If you pass a function to \code{tr_fun}, the total transition cost is measured in the units of \code{tr_fun}'s 
#' results.
#' 
#' @seealso \link{spaths_earth}
#' 
#' @export
spaths_general <- function(rst, xres, yres, xmin, ymin, origins, destinations = NULL, lonlat = TRUE, radius = NULL, output = c("lines", "distances"),
  origin_names = NULL, destination_names = NULL, pairwise = FALSE, contiguity = c("queen", "rook"), tr_fun = NULL, v_matrix = FALSE, tr_directed = TRUE,
  round_dist = FALSE, ncores = NULL, par_lvl = c("points", "points_lists"), cluster = NULL, paths_ncores = NULL, verbose = FALSE) {
  
  if(length(verbose) != 1L || !is.logical(verbose) || is.na(verbose)) stop("verbose must be logical and of length one")
  if(verbose) message("Checking arguments")
  
  # Check rst
  rst_list <- is.list(rst)
  if(rst_list) {
    if(!all(vapply(rst, function(r) any(class(rst) == "matrix"), logical(1L), USE.NAMES = FALSE))) stop("All of rst's list elements must be matrices")
    if(!do.call(identical, lapply(rst, dim))) stop("The matrices must have the same dimensions")
    nr <- NROW(rst[[1L]])
    nc <- NCOL(rst[[1L]])
  } else {
    if(!any(class(rst) == "matrix")) stop("rst must be a matrix or a list of matrices")
    nr <- NROW(rst)
    nc <- NCOL(rst)
  }
  
  # Check origin_names and destination_names
  origin_nms_null <- is.null(origin_names)
  if(!origin_nms_null && (!(is.character(origin_names) && length(origin_names) == 1L) || is.na(origin_names))) {
    stop("origin_names must be NULL or a character object of length one")
  }
  destination_nms_null <- is.null(destination_names)
  if(!destination_nms_null && (!(is.character(destination_names) && length(destination_names) == 1L) || is.na(destination_names))) {
    stop("destination_names must be NULL or a character object of length one")
  }
  
  # Check rst attributes
  if(length(xres) != 1L || xres <= 0) stop("xres must be a positive number")
  if(length(yres) != 1L || yres <= 0) stop("yres must be a positive number")
  if(length(xmin) != 1L || !is.finite(xmin)) stop("xmin must be a number")
  if(length(ymin) != 1L || !is.finite(ymin)) stop("ymin must be a number")
  if(length(lonlat) != 1L || !is.logical(lonlat) || is.na(lonlat)) stop("lonlat must be logical and of length one")
  if(lonlat) {
    if(xmin < -180 || xmin + xres * nc > 180) {
      stop("If lonlat is TRUE, xmin and xmax, i.e. xmin + xres * NCOL(rst), have to be in the [-180, 180] interval.")
    }
    if(ymin < -90 || ymin + yres * nr > 90) {
      stop("If lonlat is TRUE, ymin and ymax, i.e. ymin + yres * NROW(rst), have to be in the [-90, 90] interval.")
    }
  }
  
  # Convert origins
  origin_list <- any(class(origins) == "list")
  if(origin_list) {
    origins <- lapply(origins, convert_points_g, rst, rst_list, nr, nc, xres, yres, xmin, ymin, origin_names, origin_nms_null)
  } else {
    origins <- convert_points_g(origins, rst, rst_list, nr, nc, xres, yres, xmin, ymin, origin_names, origin_nms_null)
  }
  
  # Convert destinations
  dest_specified <- !is.null(destinations)
  if(dest_specified) {
    dest_list <- is.list(destinations)
    if(dest_list) {
      destinations <- lapply(destinations, convert_points_g, rst, rst_list, nr, nc, xres, yres, xmin, ymin, destination_names, destination_nms_null, FALSE)
    } else {
      destinations <- convert_points_g(destinations, rst, rst_list, nr, nc, xres, yres, xmin, ymin, destination_names, destination_nms_null, FALSE)
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
        stop("If pairwise is FALSE and destinations is a list, origins must also be a list. Pass origins and destinations as single vector objects, not ",
          "as lists, to compute shortest paths from all origins to all destinations.")
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
  if(rst_list) {
    crd <- lapply(rst, function(r) as.vector(t(r)))                             # Extract matrix values
    data.table::setDT(crd)
  } else {
    crd <- data.table::data.table(V1 = as.vector(t(rst)))                       # Extract matrix values
  }
  crd[, c("x", "y", "c_n") := list(rep.int(seq.int(xmin + 0.5 * xres, by = xres, length.out = nc), nr), rep(seq.int(ymin + yres * (nr - 0.5),
    by = -yres, length.out = nr), each = nc), 1:.N)]                            # Obtain cell numbers and coordinates
  crd <- stats::na.omit(crd, cols = 1:(NCOL(crd) - 3L))                         # Subset to cells with non-NA values in all matrices
  if(!(tr_fun_specified && any(args_used[6:7]))) crd[, setdiff(names(crd), c("x", "y", "c_n")) := NULL]
  data.table::setkey(crd, "c_n")
  
  # Obtain adjacency matrix
  if(lonlat) {
    if(rst_list) {
      rst <- data.table::as.data.table(terra::adjacent(terra::rast(rst[[1L]], crs = "epsg:4326", extent = terra::ext(c(xmin, xmin + xres * nc, ymin,
        ymin + yres * nr))), crd[["c_n"]], contiguity, TRUE))
    } else {
      rst <- data.table::as.data.table(terra::adjacent(terra::rast(rst, crs = "epsg:4326", extent = terra::ext(c(xmin, xmin + xres * nc, ymin,
        ymin + yres * nr))), crd[["c_n"]], contiguity, TRUE))                   
    }
  } else {
    if(rst_list) {
      rst <- data.table::as.data.table(terra::adjacent(terra::rast(rst[[1L]]), crd[["c_n"]], contiguity, TRUE))
    } else {
      rst <- data.table::as.data.table(terra::adjacent(terra::rast(rst), crd[["c_n"]], contiguity, TRUE))
    }
  }
  
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
  crd[, c("c_n", "c_n_c") := NULL]
  
  # Register parallel backend
  if(ncoresg1 && tr_fun_specified && args_used[9L]) {
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
    if(args_used[1L]) tr_fun_args$d <- compute_dists_g(rst, crd[, c("x", "y")], round_dist, contiguity, yres, xres, nr, ymin, lonlat, radius, ncores)
    if(args_used[2L]) tr_fun_args$x1 <- crd[rst[["from"]], "x"][["x"]]
    if(args_used[3L]) tr_fun_args$y1 <- crd[rst[["from"]], "y"][["y"]]
    if(args_used[4L]) tr_fun_args$x2 <- crd[rst[["to"]], "x"][["x"]]
    if(args_used[5L]) tr_fun_args$y2 <- crd[rst[["to"]], "y"][["y"]]
    if(args_used[6L]) {
      if(v_matrix) {
        tr_fun_args$v1 <- as.matrix(crd[rst[["from"]], .SD, .SDcols = setdiff(names(crd), c("x", "y"))])
      } else if(rst_list) {
        tr_fun_args$v1 <- crd[rst[["from"]], .SD, .SDcols = setdiff(names(crd), c("x", "y"))]
      } else {
        tr_fun_args$v1 <- crd[rst[["from"]], .SD, .SDcols = setdiff(names(crd), c("x", "y"))][[1L]]
      }
    }
    if(args_used[7L]) {
      if(v_matrix) {
        tr_fun_args$v2 <- as.matrix(crd[rst[["to"]], .SD, .SDcols = setdiff(names(crd), c("x", "y"))])
      } else if(rst_list) {
        tr_fun_args$v2 <- crd[rst[["to"]], .SD, .SDcols = setdiff(names(crd), c("x", "y"))]
      } else {
        tr_fun_args$v2 <- crd[rst[["to"]], .SD, .SDcols = setdiff(names(crd), c("x", "y"))][[1L]]
      }
    }
    if(any(args_used[6:7])) crd[, setdiff(names(crd), c("x", "y")) := NULL]
    if(args_used[8L]) tr_fun_args$nc <- ncores
    if(args_used[9L]) tr_fun_args$cl <- cl
    tr_fun_args <- do.call(tr_fun, tr_fun_args[tr_fun_v])
    if(!is.vector(tr_fun_args)) stop("tr_fun must return a vector")
    if(any(tr_fun_args < 0)) stop("tr_fun must not return negative values")
    if(length(tr_fun_args) != NROW(rst)) stop("The number of values returned by tr_fun must equal the number of edges")
    rst <- igraph::set_edge_attr(igraph::graph_from_edgelist(as.matrix(rst), directed = tr_directed), "weight", value = tr_fun_args) # Construct graph
    rm(tr_fun_args, args_used, tr_fun_v)
  } else {
    rst <- igraph::set_edge_attr(igraph::graph_from_edgelist(as.matrix(rst), directed = FALSE), "weight", value = compute_dists_g(rst, crd[, c("x", "y")],
      round_dist, contiguity, yres, xres, nr, ymin, lonlat, radius, ncores))    # Construct graph
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
  }
  
  rm(tr_fun_specified, paths_ncores_specified)
  
  # Compute shortest paths
  output_lines <- output == "lines"
  if(origin_list) {
    if(dest_specified) {
      if(dest_list) {
        if(ncoresg1 && par_lvl == "points_lists") {
          paths <- function(O) {
            return(compute_spaths_g(origins[[O]], rst, crd, dest_specified, origin_nms_null, output_lines, pairwise, FALSE, NULL, NULL, TRUE,
              destinations[[O]], destination_nms_null))
          }
          paths <- parallel::parLapplyLB(cl, 1:length(origins), paths)
          if(output_lines) paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]]))
        } else {
          paths <- mapply(compute_spaths_g, ORIGINS = origins, DESTINATIONS = destinations, MoreArgs = list(rst = rst, crd = crd,
            dest_specified = dest_specified, origin_nms_null = origin_nms_null, output_lines = output_lines, pairwise = pairwise, NCORESG1 = ncoresg1,
            ncores = ncores, cl = cl, nvect = FALSE, destination_nms_null = destination_nms_null), SIMPLIFY = FALSE, USE.NAMES = FALSE)
        }
      } else {
        if(ncoresg1 && par_lvl == "points_lists") {
          paths <- function(O) {
            return(compute_spaths_g(O, rst, crd, dest_specified, origin_nms_null, output_lines, pairwise, FALSE, NULL, NULL, TRUE, destinations,
              destination_nms_null))
          }
          paths <- parallel::parLapplyLB(cl, origins, paths)
          if(output_lines) paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]]))
        } else {
          paths <- lapply(origins, compute_spaths_g, rst, crd, dest_specified, origin_nms_null, output_lines, pairwise, ncoresg1, ncores, cl, FALSE,
            destinations, destination_nms_null)
        }
      }
    } else {
      paths <- lapply(origins, compute_spaths_g, rst, crd, dest_specified, origin_nms_null, output_lines, pairwise, ncoresg1, ncores, cl, FALSE)
    }
  } else {
    if(dest_specified) {
      if(dest_list) {
        if(ncoresg1 && par_lvl == "points_lists") {
          paths <- function(D) {
            return(compute_spaths_g(origins, rst, crd, dest_specified, origin_nms_null, output_lines, pairwise, FALSE, NULL, NULL, TRUE, destinations[[D]],
              destination_nms_null))
          }
          paths <- parallel::parLapplyLB(cl, 1:length(destinations), paths)
          if(output_lines) paths <- lapply(paths, function(D) terra::vect(D[[1L]], type = "line", atts = D[[2L]]))
        } else {
          paths <- lapply(destinations, function(d) compute_spaths_g(origins, rst, crd, dest_specified, origin_nms_null, output_lines, pairwise, ncoresg1,
            ncores, cl, FALSE, d, destination_nms_null))
        }
      } else {
        paths <- compute_spaths_g(origins, rst, crd, dest_specified, origin_nms_null, output_lines, pairwise, ncoresg1, ncores, cl, FALSE, destinations,
          destination_nms_null)
      }
    } else {
      paths <- compute_spaths_g(origins, rst, crd, dest_specified, origin_nms_null, output_lines, pairwise, ncoresg1, ncores, cl, FALSE)
    }
  }
  if(ncoresg1) parallel::stopCluster(cl)
  
  return(paths)
}

# Check points' location and list the corresponding grid cell numbers
convert_points_g <- function(v, rst, rst_list, nr, nc, xres, yres, xmin, ymin, nms, nms_null, o = TRUE) {
  if(!any(class(v) %chin% c("matrix", "data.frame"))) stop(ifelse(o, "origins", "destinations"), " must be a matrix or data frame, or a list of them")
  if(!nms_null) {
    if(!(nms %chin% names(v))) {
      stop(ifelse(o, "origin", "destination"), "_names must either be NULL or the name of a column in ", ifelse(o, "origins", "destinations"))
    }
    nms <- unlist(v[, nms], use.names = FALSE)
  }
  if(data.table::is.data.table(v)) {
    v <- nr * (check_position(v[, 2L][[1L]], ymin, yres, nr, o) + 1L) - check_position(v[, 1L][[1L]], xmin, xres, nc, o)
  } else {
    v <- nr * (check_position(v[, 2L], ymin, yres, nr, o) + 1L) - check_position(v[, 1L], xmin, xres, nc, o)
  }
  if(rst_list) {
    if(any(vapply(rst, function(r) any(is.na(r[v])), logical(1L), USE.NAMES = FALSE))) {
      report_points(unique(do.call(c, lapply(rst, function(r) which(is.na(r[v]))))), o)
    }
  } else if(any(is.na(rst[v]))) {
    report_points(which(is.na(rst[v])), o)
  }
  if(!nms_null) v <- data.table::data.table(cls = v, nms = nms)
  return(v)
}

# Check points' rows and columns in grid
check_position <- function(v, rmin, rres, nd, o) {
  if(!all(v %between% c(rmin, rmin + rres * nd))) {
    v <- which(!(v %between% c(rmin, rmin + rres * nd)))
    v_length <- length(v)
    if(v_length > 1L) {
      if(v_length > 2L) {
        v <- paste0(paste0(v[1:(v_length - 1L)], collapse = ", "), ", and ", v[v_length])
      } else {
        v <- paste0(v, collapse = " and ")
      }
      v <- paste0(ifelse(o, "Origins", "Destinations"), " ", v, " located outside extent of rst")
    } else {
      v <- paste0(ifelse(o, "Origin", "Destination"), " ", v, " located outside extent of rst")
    }
    stop(v)
  }
  return(floor((v - rmin) / rres))
}

# Report points on NA cells
report_points <- function(v, o) {
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

# Compute distances between neighboring cells
# Capitalized object names differing from the ones in the spaths_earth function indicate that these objects might differ from the spaths_earth
# counterparts, e.g. by being subsets
compute_dists_g <- function(rst, CRD, round_dist, contiguity, yr, xr, nr, ym, lonlat, radius, ncores) {
  if(lonlat) {
    if(is.null(radius) || radius <= 0) stop("If lonlat is TRUE, radius must be a number larger than zero")
  } else {
    radius <- 0
  }
  if(round_dist) {
    if(contiguity == "queen") {
      d <- dists_queen_i(CRD[rst[["from"]], "y"][["y"]], CRD[rst[["from"]], "x"][["x"]], CRD[rst[["to"]], "y"][["y"]], CRD[rst[["to"]], "x"][["x"]], yr,
        xr, nr, ym, lonlat, ncores, radius * 2)
    } else {
      d <- dists_rook_i(CRD[rst[["from"]], "y"][["y"]], CRD[rst[["from"]], "x"][["x"]], CRD[rst[["to"]], "y"][["y"]], CRD[rst[["to"]], "x"][["x"]], yr,
        xr, nr, ym, lonlat, ncores, radius * 2)
    }
  } else {
    if(contiguity == "queen") {
      d <- dists_queen_d(CRD[rst[["from"]], "y"][["y"]], CRD[rst[["from"]], "x"][["x"]], CRD[rst[["to"]], "y"][["y"]], CRD[rst[["to"]], "x"][["x"]], yr,
        xr, nr, ym, lonlat, ncores, radius * 2)
    } else {
      d <- dists_rook_d(CRD[rst[["from"]], "y"][["y"]], CRD[rst[["from"]], "x"][["x"]], CRD[rst[["to"]], "y"][["y"]], CRD[rst[["to"]], "x"][["x"]], yr,
        xr, nr, ym, lonlat, ncores, radius * 2)
    }
  }
  return(d)
}

# Compute shortest paths
compute_spaths_g <- function(ORIGINS, rst, crd, dest_specified, origin_nms_null, output_lines, pairwise, NCORESG1, ncores, cl, nvect, DESTINATIONS = NULL,
  destination_nms_null = TRUE) {
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
            p))[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line", atts = data.frame(origin = on, destination = dn))
        } else if(nvect) {
          p <- list(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
            return(data.table::rbindlist(lapply(igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS[O], output = "vpath", algorithm = "dijkstra")$vpath,
              function(D) crd[D,][, g := O])))
          }))[, c("g", "x", "y")]), data.frame(origin = on, destination = dn))
        } else {
          p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
            return(data.table::rbindlist(lapply(igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS[O], output = "vpath", algorithm = "dijkstra")$vpath,
              function(D) crd[D,][, g := O])))
          }))[, c("g", "x", "y")]), type = "line", atts = data.frame(origin = on, destination = dn))
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
            destination = rep.int(dn, os_length)))
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
            os_length)))
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
            p <- data.table::as.data.table(do.call(cbind, parallel::parLapply(cl, split(DESTINATIONS, cut(seq_along(DESTINATIONS), ncores,
              labels = FALSE)), p)), na.rm = FALSE)
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
          atts = stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")))
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
        }))[, c("g", "x", "y")]), type = "line", atts = stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")))
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

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
#' of what the column names are. It defaults to \code{NULL}, resulting in the function to compute shortest paths between the \code{origins} points. 
#' Passing a list of data frames or matrices is a way not to relate all origins to all destinations. Details are outlined below.
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
#' @param cluster \code{"PSOCK"}, \code{"FORK"}, or \code{"MPI"}, indicating the type of \code{parallel} cluster that the function employs when 
#' \code{ncores > 1} or \code{paths_ncores > 1}. The function defaults to \code{"PSOCK"} on Windows and to \code{"FORK"} otherwise. \code{"FORK"} is not 
#' available on Windows machines. \code{"MPI"} requires the Rmpi package to be installed. There are various ways of parallelizing R with MPI. This 
#' function utilizes the variant implemented in the \code{snow} package.
#' @param paths_ncores An integer specifying the number of CPU cores to use in shortest paths computations. It defaults to the value of \code{ncores}. 
#' Thus, only set it, if you want edge weights and shortest paths be computed with differently many cores. The \code{dist_comp = "spaths"} edge weight
#' computations employ efficient C++ level parallelization. The shortest paths sections, in contrast, parallelize on the R level. If you use a \code{PSOCK} 
#' cluster, \code{spaths_earth} copies various objects to the workers before the paths algorithm is applied. This can make the parallel execution slower 
#' than its serial counterpart. Thus, consider setting \code{paths_ncores = 1}, especially when working with \code{PSOCK} clusters.
#' @param write_dir Directory to which write output. If \code{NULL} (default), the output is not written to disk but returned by the function. When 
#' specifying a directory, \code{spaths_general} writes the output to it and returns \code{NULL}. Keeping the results in RAM is commonly faster than 
#' writing them to disk. Thus, the recommendation is to keep the default unless your machine has insufficient RAM. In cases in which the function would 
#' return a list, i.e. if \code{origins} or \code{destinations} are lists, \code{spaths_general} writes one file per list element. File names state which 
#' element the file represents: the number following \code{p} points to the element of \code{origins} or \code{destinations}. Values begin at 1. In 
#' generating file names of equal length, the function employs leading zeros. So, if the function computes 100 SpatVectors, the file corresponding to the 
#' first list element is not \code{p1} but \code{p001}. In the basic case, in which the result is a single data.table or SpatVector object, the output 
#' file is named \code{results}. The directory must not contain any file named as one of the output files - irrespective of file type - if 
#' \code{output = "lines"}.
#' @param file_type The output file type when \code{write_dir} is specified. It defaults to \code{"shp"} in case \code{output = "lines"}, but can also 
#' be set to \code{"kml"}, \code{"json"}, or any other vector format that \code{terra::writeVector} can write. \code{terra::gdal(drivers = TRUE)} lists 
#' drivers. \code{output = "distances"} defaults to \code{"rds"} as file type, but can alternatively be set to \code{"csv"}.
#' @param unconnected_error Logical specifying whether the function throws an error when trying to compute the distance between locations unconnected by 
#' the \code{rst} grid. If \code{TRUE} (default), the function throws and error and aborts. If \code{FALSE}, unconnected places have an \code{Inf} 
#' distance. If \code{output = "lines"}, attempting to link unconnected places always throws an error.
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
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(20146)
#' rst <- matrix(sample(c(1L, NA), 64800L, TRUE, c(0.9, 0.1)), 180L, 360L)
#' origins <- rnd_locations(5L)
#' destinations <- rnd_locations(5L)
#' 
#' # Compute shortest paths
#' spaths_general(rst, 1, 1, -180, -90, origins, radius = 6000)
#' spaths_general(rst, 1, 1, -180, -90, origins, destinations, radius = 6000)
#' }
#' 
#' @export
spaths_general <- function(rst, xres, yres, xmin, ymin, origins, destinations = NULL, lonlat = TRUE, radius = NULL, output = c("lines", "distances"),
  origin_names = NULL, destination_names = NULL, pairwise = FALSE, contiguity = c("queen", "rook"), tr_fun = NULL, v_matrix = FALSE, tr_directed = TRUE,
  round_dist = FALSE, ncores = NULL, par_lvl = c("points", "points_lists"), cluster = NULL, paths_ncores = NULL, write_dir = NULL, file_type = NULL,
  unconnected_error = TRUE, verbose = FALSE) {
  
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
  origin_nms_specified <- !is.null(origin_names)
  if(origin_nms_specified && (!(is.character(origin_names) && length(origin_names) == 1L) || is.na(origin_names))) {
    stop("origin_names must be NULL or a character object of length one")
  }
  destination_nms_specified <- !is.null(destination_names)
  if(destination_nms_specified && (!(is.character(destination_names) && length(destination_names) == 1L) || is.na(destination_names))) {
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
    origins <- lapply(origins, convert_points_g, rst, rst_list, nr, nc, xres, yres, xmin, ymin, origin_names, origin_nms_specified)
  } else {
    origins <- convert_points_g(origins, rst, rst_list, nr, nc, xres, yres, xmin, ymin, origin_names, origin_nms_specified)
  }
  
  # Convert destinations
  dest_specified <- !is.null(destinations)
  if(dest_specified) {
    dest_list <- any(class(destinations) == "list")
    if(dest_list) {
      destinations <- lapply(destinations, convert_points_g, rst, rst_list, nr, nc, xres, yres, xmin, ymin, destination_names, destination_nms_specified,
        FALSE)
    } else {
      destinations <- convert_points_g(destinations, rst, rst_list, nr, nc, xres, yres, xmin, ymin, destination_names, destination_nms_specified, FALSE)
    }
  } else {
    dest_list <- NULL
  }
  p_list <- origin_list || (dest_specified && dest_list)
  
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
  output_lines <- output == "lines"
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
        nfork <- TRUE
      } else {
        cluster <- "FORK"
        nfork <- FALSE
      }
    } else if(length(cluster) != 1L || !(cluster %chin% c("PSOCK", "FORK", "MPI"))) {
      stop('cluster must be NULL, "PSOCK", "FORK", or "MPI"')
    } else if(cluster == "MPI" && length(find.package("Rmpi", quiet = TRUE)) == 0L) {
      stop('Install the Rmpi package to use cluster = "MPI"')
    } else {
      nfork <- cluster != "FORK"
    }
  } else {
    nfork <- NULL
  }
  # par_lvl is checked before shortest paths computation
  write_disk <- !is.null(write_dir)
  if(write_disk) {
    if(length(write_dir) != 1L || !is.character(write_dir) || !dir.exists(write_dir)) stop("write_dir must be NULL or the path of an existent directory")
    if(is.null(file_type)) {
      if(output_lines) {
        file_type <- "shp"
        file_type_rds <- NULL
      } else {
        file_type_rds <- TRUE
      }
    } else {
      if(length(file_type) != 1L || !is.character(file_type)) stop("file_type must be NULL or a character string")
      if(output_lines) {
        file_type_rds <- NULL
      } else {
        if(!(file_type %chin% c("rds", "csv"))) stop('file_type must be NULL, "rds", or "csv" if output = "distances"')
        file_type_rds <- file_type == "rds"
      }
    }
    if(p_list) {
      if(origin_list) {
        wp <- nchar(length(origins))
      } else {
        wp <- nchar(length(destinations))
      }
    }
    if(output_lines && length(list.files(write_dir, "^results[.]")) > 0L) {
      stop("results already exists in ", write_dir)
    }
  }
  if(length(unconnected_error) != 1L || !is.logical(unconnected_error) || is.na(unconnected_error)) {
    stop("unconnected_error must be logical and of length one")
  }
  
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
  data.table::setkey(crd, c_n)
  
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
  
  origins <- update_points(origins, origin_list, crd, origin_nms_specified)
  if(dest_specified) destinations <- update_points(destinations, dest_list, crd, destination_nms_specified)
  
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
  if(ncoresg1 && tr_fun_specified && args_used[9L] && nfork) {
    if(cluster == "PSOCK") {
      cl <- parallel::makePSOCKcluster(ncores, useXDR = FALSE)
    } else {
      cl <- parallel::makeCluster(type = "MPI")
    }
    parallel::clusterEvalQ(cl, igraph::igraph_options(return.vs.es = FALSE))
    on.exit({
      if(!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
    }, add = TRUE)
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
    if(!is.vector(tr_fun_args)) {
      if(!is.null(cl)) parallel::stopCluster(cl)
      stop("tr_fun must return a vector")
    }
    if(any(tr_fun_args < 0)) {
      if(!is.null(cl)) parallel::stopCluster(cl)
      stop("tr_fun must not return negative values")
    }
    if(length(tr_fun_args) != NROW(rst)) {
      if(!is.null(cl)) parallel::stopCluster(cl)
      stop("The number of values returned by tr_fun must equal the number of edges")
    }
    rst <- igraph::set_edge_attr(igraph::graph_from_edgelist(as.matrix(rst), directed = tr_directed), "weight", value = tr_fun_args) # Construct graph
    rm(tr_fun_args, args_used, tr_fun_v)
  } else {
    rst <- igraph::set_edge_attr(igraph::graph_from_edgelist(as.matrix(rst), directed = FALSE), "weight", value = compute_dists_g(rst, crd[, c("x", "y")],
      round_dist, contiguity, yres, xres, nr, ymin, lonlat, radius, ncores))    # Construct graph
    tr_directed <- FALSE
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
    par_lvl <- match(match.arg(par_lvl), c("points", "points_lists"))
    if(!origin_list && dest_specified && !dest_list && par_lvl == 2L) {
      if(!is.null(cl)) parallel::stopCluster(cl)
      stop('If par_lvl is "points_lists", origins or destinations must be a list')
    }
    if(is.null(cl) && nfork) {
      if(cluster == "PSOCK") {
        cl <- parallel::makePSOCKcluster(ncores, useXDR = FALSE)
      } else {
        cl <- parallel::makeCluster(type = "MPI")
      }
      parallel::clusterEvalQ(cl, igraph::igraph_options(return.vs.es = FALSE))
      on.exit({
        if(!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
      }, add = TRUE)
    }
  } else {
    par_lvl <- 0L
  }
  
  rm(tr_fun_specified, paths_ncores_specified)
  
  # Compute shortest paths
  if(write_disk && p_list) {
    write_file <- function(o, P) {
      if(output_lines) {
        terra::writeVector(o, paste0(write_dir, "/", formatC(P, flag = "0", width = wp), ".", file_type))
      } else if(file_type_rds) {
        saveRDS(o, paste0(write_dir, "/", formatC(P, flag = "0", width = wp), ".rds"))
      } else {
        data.table::fwrite(o, paste0(write_dir, "/", formatC(P, flag = "0", width = wp), ".csv"))
      }
      return(NULL)
    }
  }
  if(origin_list) {
    if(dest_specified) {
      if(dest_list) {
        if(ncoresg1 && par_lvl == 2L) {
          if(write_disk) {
            paths <- function(O, D, P) {
              return(write_file(compute_spaths_g(O, rst, crd, TRUE, origin_nms_specified, output_lines, pairwise, FALSE, NULL, NULL, NULL, FALSE,
                unconnected_error, tr_directed, D, destination_nms_specified), P))
            }
            if(nfork) {
              parallel::clusterMap(cl, paths, origins, destinations, 1:length(origins), USE.NAMES = FALSE, .scheduling = "dynamic")
            } else {
              parallel::mcmapply(paths, origins, destinations, 1:length(origins), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE,
                mc.silent = TRUE, mc.cores = ncores)
            }
          } else {
            paths <- function(O, D) {
              return(compute_spaths_g(O, rst, crd, TRUE, origin_nms_specified, output_lines, pairwise, FALSE, NULL, FALSE, NULL, TRUE,
                unconnected_error, tr_directed, D, destination_nms_specified))
            }
            if(nfork) {
              paths <- parallel::clusterMap(cl, paths, origins, destinations, USE.NAMES = FALSE, .scheduling = "dynamic")
            } else {
              paths <- parallel::mcmapply(paths, origins, destinations, SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE, mc.silent = TRUE,
                mc.cores = ncores)
            }
            if(output_lines) paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]]))
          }
        } else {
          if(write_disk) {
            mapply(function(O, D, P) write_file(compute_spaths_g(O, rst, crd, TRUE, origin_nms_specified, output_lines, pairwise, ncoresg1, ncores, nfork,
              cl, FALSE, unconnected_error, tr_directed, D, destination_nms_specified), P), origins, destinations, 1:length(origins), SIMPLIFY = FALSE,
              USE.NAMES = FALSE)
          } else {
            paths <- mapply(compute_spaths_g, ORIGINS = origins, DESTINATIONS = destinations, MoreArgs = list(rst = rst, crd = crd,
              dest_specified = dest_specified, origin_nms_specified = origin_nms_specified, output_lines = output_lines, pairwise = pairwise,
              NCORESG1 = ncoresg1, ncores = ncores, nfork = nfork, cl = cl, nvect = FALSE, unconnected_error = unconnected_error,
              tr_directed = tr_directed, destination_nms_specified = destination_nms_specified), SIMPLIFY = FALSE, USE.NAMES = FALSE)
          }
        }
      } else {
        if(ncoresg1 && par_lvl == 2L) {
          if(write_disk) {
            paths <- function(O, P) {
              return(write_file(compute_spaths_g(O, rst, crd, TRUE, origin_nms_specified, output_lines, pairwise, FALSE, NULL, NULL, NULL, FALSE,
                unconnected_error, tr_directed, destinations, destination_nms_specified), P))
            }
            if(nfork) {
              parallel::clusterMap(cl, paths, origins, 1:length(origins), USE.NAMES = FALSE, .scheduling = "dynamic")
            } else {
              parallel::mcmapply(paths, origins, 1:length(origins), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE, mc.silent = TRUE,
                mc.cores = ncores)
            }
          } else {
            paths <- function(O) {
              return(compute_spaths_g(O, rst, crd, TRUE, origin_nms_specified, output_lines, pairwise, FALSE, NULL, NULL, NULL, TRUE,
                unconnected_error, tr_directed, destinations, destination_nms_specified))
            }
            if(nfork) {
              paths <- parallel::parLapplyLB(cl, origins, paths)
            } else {
              paths <- parallel::mclapply(origins, paths, mc.preschedule = FALSE, mc.silent = TRUE, mc.cores = ncores)
            }
            if(output_lines) paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]]))
          }
        } else {
          if(write_disk) {
            mapply(function(O, P) write_file(compute_spaths_g(O, rst, crd, TRUE, origin_nms_specified, output_lines, pairwise, ncoresg1,
              ncores, nfork, cl, FALSE, unconnected_error, tr_directed, destinations, destination_nms_specified), P), origins, 1:length(origins),
              SIMPLIFY = FALSE, USE.NAMES = FALSE)
          } else {
            paths <- lapply(origins, compute_spaths_g, rst, crd, TRUE, origin_nms_specified, output_lines, pairwise, ncoresg1, ncores, nfork,
              cl, FALSE, unconnected_error, tr_directed, destinations, destination_nms_specified)
          }
        }
      }
    } else {
      if(write_disk) {
        mapply(function(O, P) write_file(compute_spaths_g(O, rst, crd, FALSE, origin_nms_specified, output_lines, pairwise, ncoresg1, ncores, nfork, cl,
          FALSE, unconnected_error, tr_directed), P), origins, 1:length(origins), SIMPLIFY = FALSE, USE.NAMES = FALSE)
      } else {
        paths <- lapply(origins, compute_spaths_g, rst, crd, FALSE, origin_nms_specified, output_lines, pairwise, ncoresg1, ncores, nfork, cl, FALSE,
          unconnected_error, tr_directed)
      }
    }
  } else {
    if(dest_specified) {
      if(dest_list) {
        if(ncoresg1 && par_lvl == 2L) {
          if(write_disk) {
            paths <- function(D, P) {
              return(write_file(compute_spaths_g(origins, rst, crd, TRUE, origin_nms_specified, output_lines, pairwise, FALSE, NULL, NULL, NULL, FALSE,
                unconnected_error, tr_directed, D, destination_nms_specified), P))
            }
            if(nfork) {
              parallel::clusterMap(cl, paths, destinations, 1:length(destinations), USE.NAMES = FALSE, .scheduling = "dynamic")
            } else {
              parallel::mcmapply(paths, destinations, 1:length(destinations), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE,
                mc.silent = TRUE, mc.cores = ncores)
            }
          } else {
            paths <- function(D) {
              return(compute_spaths_g(origins, rst, crd, TRUE, origin_nms_specified, output_lines, pairwise, FALSE, NULL, NULL, NULL, TRUE,
                unconnected_error, tr_directed, D, destination_nms_specified))
            }
            if(nfork) {
              paths <- parallel::parLapplyLB(cl, destinations, paths)
            } else {
              paths <- parallel::mclapply(destinations, paths, mc.preschedule = FALSE, mc.silent = TRUE, mc.cores = ncores)
            }
            if(output_lines) paths <- lapply(paths, function(D) terra::vect(D[[1L]], type = "line", atts = D[[2L]]))
          }
        } else {
          if(write_disk) {
            mapply(function(d, P) write_file(compute_spaths_g(origins, rst, crd, dest_specified, origin_nms_specified, output_lines, pairwise, ncoresg1,
              ncores, fork, cl, FALSE, unconnected_error, tr_directed, d, destination_nms_specified), P), destinations, 1:length(destinations),
              SIMPLIFY = FALSE, USE.NAMES = FALSE)
          } else {
            paths <- lapply(destinations, function(d) compute_spaths_g(origins, rst, crd, dest_specified, origin_nms_specified, output_lines, pairwise,
              ncoresg1, ncores, fork, cl, FALSE, unconnected_error, tr_directed, d, destination_nms_specified))
          }
        }
      } else {
        paths <- compute_spaths_g(origins, rst, crd, dest_specified, origin_nms_specified, output_lines, pairwise, ncoresg1, ncores, nfork, cl, FALSE,
          unconnected_error, tr_directed, destinations, destination_nms_specified)
        if(write_disk) {
          if(output_lines) {
            terra::writeVector(paths, paste0(write_dir, "/results.", file_type))
          } else if(file_type_rds) {
            saveRDS(paths, paste0(write_dir, "/results.rds"))
          } else {
            data.table::fwrite(paths, paste0(write_dir, "/results.csv"))
          }
        }
      }
    } else {
      paths <- compute_spaths_g(origins, rst, crd, FALSE, origin_nms_specified, output_lines, pairwise, ncoresg1, ncores, nfork, cl, FALSE,
        unconnected_error, tr_directed)
      if(write_disk) {
        if(output_lines) {
          terra::writeVector(paths, paste0(write_dir, "/results.", file_type))
        } else if(file_type_rds) {
          saveRDS(paths, paste0(write_dir, "/results.rds"))
        } else {
          data.table::fwrite(paths, paste0(write_dir, "/results.csv"))
        }
      }
    }
  }
  if(ncoresg1 && nfork && !is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- NULL
  }
  if(write_disk) paths <- NULL
  return(paths)
}

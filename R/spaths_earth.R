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
#' \code{parallel::makeCluster(type = "MPI")}). If \code{lonlat} is \code{TRUE} or if \code{dist_comp = "terra"}, \code{d} is measured in the meters. 
#' Otherwise, it uses the units of the CRS. \code{nc} is meant to be used in C++ functions or in \code{FORK}-based parallelization. A parallel backend at 
#' the R level pointed to by \code{cl} is already registered before \code{tr_fun} is called, in case \code{ncores > 1} and \code{cluster} is not 
#' \code{"FORK"}. Use \code{multicore}-based functions like \code{parallel::mclapply} rather than their less efficient \code{snow}-based counterparts 
#' like \code{parallel::parLapply} in the \code{"FORK"} case. If \code{rst} has one layer, the values are passed to \code{v1} and \code{v2} as vectors, 
#' otherwise they are passed as a data table where the first column refers to the first layer, the second column to the second layer etc. Note that data 
#' tables are also data frames.
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
#' @param cluster \code{"PSOCK"}, \code{"FORK"}, or \code{"MPI"}, indicating the type of \code{parallel} cluster that the function employs when 
#' \code{ncores > 1} or \code{paths_ncores > 1}. The function defaults to \code{"PSOCK"} on Windows and to \code{"FORK"} otherwise. \code{"FORK"} is not 
#' available on Windows machines. \code{"MPI"} requires the Rmpi package to be installed. There are various ways of parallelizing R with MPI. This 
#' function utilizes the variant implemented in the \code{snow} package.
#' @param paths_ncores Integer specifying the number of CPU cores to use in shortest paths computations. It defaults to the value of \code{ncores}. 
#' Thus, only set it, if you want edge weights and shortest paths be computed with differently many cores. The \code{dist_comp = "spaths"} edge weight
#' computations employ efficient C++ level parallelization. The shortest paths sections, in contrast, parallelize on the R level. If you use a \code{PSOCK} 
#' cluster, \code{spaths_earth} copies various objects to the workers before the paths algorithm is applied. This can make the parallel execution slower 
#' than its serial counterpart. Thus, consider setting \code{paths_ncores = 1}, especially when working with \code{PSOCK} clusters.
#' @param write_dir Directory to which write output. If \code{NULL} (default), the output is not written to disk but returned by the function. When 
#' specifying a directory, \code{spaths_earth} writes the output to it and returns \code{NULL}. Keeping the results in RAM is commonly faster than 
#' writing them to disk. Thus, the recommendation is to keep the default unless your machine has insufficient RAM. In cases in which the function would 
#' return a list, i.e. if \code{update_rst} is specified or if \code{origins} or \code{destinations} are lists, \code{spaths_earth} writes one file per 
#' list element. File names state which element the file represents: the number following \code{u} points to the element of \code{update_rst} and the 
#' number following \code{p} to the element of \code{origins} or \code{destinations}. \code{u} values start at 0 where 0 features results based on 
#' \code{rst} not updated by \code{update_rst}. \code{p} values begin at 1. Hence, when all \code{update_rst}, \code{origins}, and \code{destinations} 
#' are lists, the first element (\code{[[1]][[1]]}), which represents the paths between the points in the first origins list element and the points in 
#' the first destinations list element based on the not updated graph, is written to \code{u0_p1}. In generating file names of equal length, the function 
#' employs leading zeros. So, if \code{update_rst} is a list of 100 SpatVectors, the first file is not \code{u0} but \code{u000}. In the basic case, in 
#' which the result is a single data.table or SpatVector object, the output file is named \code{results}. The directory must not contain any file named 
#' as one of the output files - irrespective of file type - if \code{output = "lines"}.
#' @param file_type The output file type when \code{write_dir} is specified. It defaults to \code{"shp"} in case \code{output = "lines"}, but can also 
#' be set to \code{"kml"}, \code{"json"}, or any other vector format that \code{terra::writeVector} can write. \code{terra::gdal(drivers = TRUE)} lists 
#' drivers. \code{output = "distances"} defaults to \code{"rds"} as file type, but can alternatively be set to \code{"csv"}.
#' @param unconnected_error Logical specifying whether the function throws an error when trying to compute the distance between locations unconnected by 
#' the \code{rst} grid. If \code{TRUE} (default), the function throws and error and aborts. If \code{FALSE}, unconnected places have an \code{Inf} 
#' distance. The argument only has an effect when \code{update_rst = NULL} or \code{copy = FALSE}. Otherwise, i.e. if \code{output = "lines"} or 
#' \code{update_rst} is not \code{NULL} and \code{copy = TRUE}, attempting to link unconnected places always throws an error.
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
#' \href{../doc/spaths_introduction.html}{vignette} for further details, such as the object's number of rows.
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
spaths_earth <- function(rst, origins, destinations = NULL, output = c("lines", "distances"), origin_names = NULL, destination_names = NULL,
  pairwise = FALSE, contiguity = c("queen", "rook"), dist_comp = c("spaths", "terra"), tr_fun = NULL, v_matrix = FALSE, tr_directed = TRUE,
  update_rst = NULL, touches = TRUE, copy = TRUE, round_dist = FALSE, ncores = NULL, par_lvl = c("points", "points_lists", "update_rst"),
  rst_par_lvl = c("points", "points_lists", "none"), cluster = NULL, paths_ncores = NULL, write_dir = NULL, file_type = NULL, unconnected_error = TRUE,
  verbose = FALSE) {
  
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
        if(any(vapply(origins, function(O) NROW(O) < 2L, logical(1L), USE.NAMES = FALSE))) {
          stop("Each origins list element must have more than one row, if destinations are not specified")
        }
      } else if(NROW(origins) < 2L) {
        stop("origins must have more than one row, if destinations are not specified")
      }
    }
  }
  output <- match.arg(output)
  output_lines <- output == "lines"
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
    if(copy) rst_par_lvl <- match(match.arg(rst_par_lvl), c("points", "points_lists", "none"))
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
    if(upd_rst_specified) {
      if(update_rst_list) {
        wu <- nchar(length(update_rst))
      } else {
        wu <- 1L
      }
      u0_n <- paste0("u", formatC(0L, flag = "0", width = wu))
      if(output_lines) {
        if(p_list) {
          if(length(list.files(write_dir, paste0("^", u0_n, "_p", formatC(1L, flag = "0", width = wp), "[.]"))) > 0L) {
            stop(paste0(u0_n, "_p", formatC(1L, flag = "0", width = wp)), " already exists in ", write_dir)
          }
        } else if(length(list.files(write_dir, paste0("^", u0_n, "[.]"))) > 0L) {
          stop(u0_n, " already exists in ", write_dir)
        }
      }
    } else if(output_lines && length(list.files(write_dir, "^results[.]")) > 0L) {
      stop("results already exists in ", write_dir)
    }
  }
  if(length(unconnected_error) != 1L || !is.logical(unconnected_error) || is.na(unconnected_error)) {
    stop("unconnected_error must be logical and of length one")
  }
  
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
  data.table::setkey(crd, c_n)
  
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
  if(ncoresg1 && ((tr_fun_specified && ((args_used[1L] && dist_comp_terra) || args_used[9L])) || (!tr_fun_specified && dist_comp_terra)) && nfork) {
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
    if(length(tr_fun_args) != NROW(rst)) {
      if(!is.null(cl)) parallel::stopCluster(cl)
      stop("The number of values returned by tr_fun must equal the number of edges")
    }
    rst <- igraph::set_edge_attr(igraph::graph_from_edgelist(as.matrix(rst), directed = tr_directed), "weight", value = tr_fun_args) # Construct graph
    rm(tr_fun_args, args_used, tr_fun_v)
  } else {
    rst <- igraph::set_edge_attr(igraph::graph_from_edgelist(as.matrix(rst), directed = FALSE), "weight", value = compute_dists(rst, crd[, c("x", "y")],
      dist_comp_terra, round_dist, contiguity, yr, xr, nr, ym, lonlat, ncoresg1, ncores, cl)) # Construct graph
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
    par_lvl <- match(match.arg(par_lvl), c("points", "points_lists", "update_rst"))
    if(par_lvl == 3L) {
      if(!upd_rst_specified) {
        if(!is.null(cl)) parallel::stopCluster(cl)
        stop('If par_lvl is "update_rst", update_rst must not be NULL')
      }
      if(copy && !update_rst_list) {
        if(!is.null(cl)) parallel::stopCluster(cl)
        stop('If copy is TRUE and update_rst is not a list, par_lvl = "update_rst" does not induce parallelization')
      }
    }
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
  
  rm(yr, xr, nr, ym, tr_fun_specified, paths_ncores_specified, dist_comp_terra)
  
  # Compute shortest paths
  if(upd_rst_specified) {
    if(update_rst_list) {
      update_rst <- lapply(update_rst, function(V) crd[.(terra::extract(rst_upd, V, cells = TRUE, ID = FALSE, touches = touches)$cell), nomatch = NULL,
        which = TRUE, on = "c_n"])
    } else {
      update_rst <- crd[.(terra::extract(rst_upd, update_rst, cells = TRUE, ID = FALSE, touches = touches)$cell), nomatch = NULL, which = TRUE, on = "c_n"]
    }
    rm(rst_upd)
    crd[, c_n := NULL]
    if((update_rst_list && any(lengths(update_rst) > 0L)) || (!update_rst_list && length(update_rst) > 0L)) {
      if(copy) {
        if(rst_par_lvl == 3L) {
          p1 <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list, dest_list,
            r_crs, output_lines, pairwise, FALSE, copy = TRUE, tr_directed = tr_directed)
        } else {
          p1 <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list, dest_list,
            r_crs, output_lines, pairwise, ncoresg1, ncores, rst_par_lvl, nfork, cl, TRUE, tr_directed = tr_directed)
        }
        if(p_list) {
          p2 <- lapply(p1, `[[`, 2L)
          if(output_lines) p3 <- lapply(p1, `[[`, 3L)
          p1 <- lapply(p1, `[[`, 1L)
          if(write_disk) p_l <- length(p1)
        } else {
          p2 <- p1[[2L]]
          if(output_lines) p3 <- p1[[3L]]
          p1 <- p1[[1L]]
        }
        if(ncoresg1 && nfork) parallel::clusterEvalQ(cl, gc(FALSE))
        # Lines output
        if(output_lines) {
          if(p_list) {
            if(write_disk) {
              mapply(function(P1, P3, P) terra::writeVector(terra::vect(as.matrix(P1[, c("g", "x", "y")]), type = "line", atts = P3, crs = r_crs),
                paste0(write_dir, "/", u0_n, "_p", formatC(P, flag = "0", width = wp), ".", file_type)), p1, p3, 1:p_l, SIMPLIFY = FALSE,
                USE.NAMES = FALSE)
              wu2 <- wu + 2L
            }
            paths <- function(u, U = NULL) {
              p_affected <- lapply(p1, function(p) unique(p[.(u), "g", nomatch = NULL, on = "cls"][["g"]]))
              if(write_disk) u_n <- formatC(U, flag = "0", width = wu)
              if(any(lengths(p_affected) > 0L)) {
                rst_u <- igraph::delete_vertices(rst, u)
                v <- data.table::data.table(c_n_c = 1:NROW(crd))[-u,]
                data.table::setkey(v, c_n_c)
                p <- function(P_AFFECTED, P1, P2, P3 = NULL, I = NULL) {
                  if(length(P_AFFECTED) > 0L) {
                    p_c <- data.table::data.table(origin_c = P2[P_AFFECTED, "origin_c"], destination_c = P2[P_AFFECTED, "destination_c"], pa = P_AFFECTED)
                    if(any(p_c[["origin_c"]] %in% u)) {
                      report_points_ust(p_c[["origin_c"]], u, TRUE)
                    }
                    if(any(p_c[["destination_c"]] %in% u)) {
                      report_points_ust(p_c[["destination_c"]], u, FALSE)
                    }
                    p_c[, c("origin_c", "destination_c") := list(v[p_c[, "origin_c"], nomatch = NULL, which = TRUE, on = c(c_n_c = "origin_c")],
                      v[p_c[, "destination_c"], nomatch = NULL, which = TRUE, on = c(c_n_c = "destination_c")])]
                    origin_n <- is.integer(p_c[["origin_c"]])
                    p_c <- split(p_c, by = "origin_c", keep.by = FALSE)
                    P <- crd[-u, c("x", "y")]
                    if(ncoresg1 && par_lvl == 1L) {
                      V <- function(dp, oc) {
                        ps <- igraph::shortest_paths(rst_u, oc, dp[["destination_c"]], output = "vpath", algorithm = "dijkstra")$vpath
                        if(min(lengths(ps), na.rm = TRUE) == 0L) {
                          if(is.null(U)) {
                            stop("Not all points are connected when updating rst with update_rst")
                          } else {
                            stop("Not all points are connected when updating rst with element ", U, " of update_rst")
                          }
                        }
                        return(data.table::rbindlist(mapply(function(pd, pa) P[pd,][, g := pa], ps, dp[["pa"]], SIMPLIFY = FALSE, USE.NAMES = FALSE),
                          use.names = FALSE))
                      }
                      l_p_c <- length(p_c)
                      if(l_p_c < ncores) {
                        if(origin_n) {
                          origin_n <- as.integer(names(p_c))
                        } else {
                          origin_n <- as.numeric(names(p_c))
                        }
                        if(l_p_c > 1L) {
                          if(length(P_AFFECTED) - l_p_c + 1 > .Machine$integer.max) {
                            origin_m <- which.max(sapply(p_c, NROW, USE.NAMES = FALSE))
                          } else {
                            origin_m <- which.max(vapply(p_c, NROW, integer(1L), USE.NAMES = FALSE))
                          }
                          n_m <- NROW(p_c[[origin_m]])
                          l_p_c <- min(c(ncores - l_p_c + 1, n_m))
                          p_c <- c(p_c[-origin_m], split(p_c[[origin_m]], cut(1:n_m, l_p_c, labels = FALSE)))
                          origin_n <- c(origin_n[-origin_m], rep.int(origin_n[origin_m], l_p_c))
                          rm(l_p_c, origin_m, n_m)
                          if(nfork) {
                            P <- rbind(P1[!.(P_AFFECTED), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::clusterMap(cl, V, p_c, origin_n,
                              USE.NAMES = FALSE), use.names = FALSE)[, c("g", "x", "y")])
                          } else {
                            P <- rbind(P1[!.(P_AFFECTED), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::mcmapply(V, p_c, origin_n,
                              SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores), use.names = FALSE)[, c("g", "x", "y")])
                          }
                        } else {
                          rm(l_p_c)
                          p_c <- split(p_c[[1L]], cut(1:NROW(p_c[[1L]]), ncores, labels = FALSE))
                          if(nfork) {
                            P <- rbind(P1[!.(P_AFFECTED), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::parLapply(cl, p_c, V, origin_n),
                              use.names = FALSE)[, c("g", "x", "y")])
                          } else {
                            P <- rbind(P1[!.(P_AFFECTED), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::mclapply(p_c, V, origin_n,
                              mc.silent = TRUE, mc.cores = ncores), use.names = FALSE)[, c("g", "x", "y")])
                          }
                        }
                      } else {
                        if(nfork) {
                          P <- rbind(P1[!.(P_AFFECTED), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::clusterMap(cl, V, p_c,
                            ifelse(origin_n, as.integer(names(p_c)), as.numeric(names(p_c))), USE.NAMES = FALSE, .scheduling = ifelse(l_p_c == ncores,
                            "static", "dynamic")), use.names = FALSE)[, c("g", "x", "y")])
                        } else {
                          P <- rbind(P1[!.(P_AFFECTED), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::mcmapply(V, p_c, ifelse(origin_n,
                            as.integer(names(p_c)), as.numeric(names(p_c))), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = l_p_c == ncores,
                            mc.silent = TRUE, mc.cores = ncores), use.names = FALSE)[, c("g", "x", "y")])
                        }
                      }
                      data.table::setorder(P, g)
                      P <- terra::vect(as.matrix(P), type = "line", atts = P3, crs = r_crs)
                    } else {
                      P <- rbind(P1[!.(P_AFFECTED), c("g", "x", "y"), on = "g"], data.table::rbindlist(mapply(function(dp, oc) {
                        ps <- igraph::shortest_paths(rst_u, oc, dp[["destination_c"]], output = "vpath", algorithm = "dijkstra")$vpath
                        if(min(lengths(ps), na.rm = TRUE) == 0L) {
                          if(is.null(U)) {
                            stop("Not all points are connected when updating rst with update_rst")
                          } else {
                            stop("Not all points are connected when updating rst with element ", U, " of update_rst")
                          }
                        }
                        return(data.table::rbindlist(mapply(function(pd, pa) P[pd,][, g := pa], ps, dp[["pa"]], SIMPLIFY = FALSE, USE.NAMES = FALSE),
                          use.names = FALSE))
                      }, p_c, ifelse(origin_n, as.integer(names(p_c)), as.numeric(names(p_c))), SIMPLIFY = FALSE, USE.NAMES = FALSE),
                        use.names = FALSE)[, c("g", "x", "y")])
                      data.table::setorder(P, g)
                      if(!ncoresg1 || write_disk) P <- terra::vect(as.matrix(P), type = "line", atts = P3, crs = r_crs)
                    }
                    if(write_disk) {
                      terra::writeVector(P, paste0(write_dir, "/u", u_n, "_p", formatC(I, flag = "0", width = wp), ".", file_type))
                      P <- NULL
                    }
                  } else {
                    if(write_disk) {
                      P <- substring(list.files(write_dir, paste0("^", u0_n, "_p", formatC(I, flag = "0", width = wp), "[.]")), wu2)
                      file.copy(paste0(write_dir, "/", u0_n, P), paste0(write_dir, "/u", u_n, P))
                      P <- NULL
                    } else {
                      P <- P1[, c("g", "x", "y")]
                    }
                  }
                  return(P)
                }
                if(write_disk) {
                  if(ncoresg1 && par_lvl == 2L) {
                    if(nfork) {
                      parallel::clusterMap(cl, p, p_affected, p1, p2, p3, 1:p_l, USE.NAMES = FALSE, .scheduling = "dynamic")
                    } else {
                      parallel::mcmapply(p, p_affected, p1, p2, p3, 1:p_l, SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE, mc.silent = TRUE,
                        mc.cores = ncores)
                    }
                  } else {
                    mapply(p, p_affected, p1, p2, p3, 1:p_l, SIMPLIFY = FALSE, USE.NAMES = FALSE)
                  }
                  v <- NULL
                } else {
                  if(ncoresg1 && par_lvl == 2L) {
                    if(nfork) {
                      v <- parallel::clusterMap(cl, p, p_affected, p1, p2, USE.NAMES = FALSE, .scheduling = "dynamic")
                    } else {
                      v <- parallel::mcmapply(p, p_affected, p1, p2, SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE, mc.silent = TRUE,
                        mc.cores = ncores)
                    }
                    v <- mapply(function(P, P3) terra::vect(as.matrix(P), type = "line", atts = P3, crs = r_crs), v, p3, SIMPLIFY = FALSE,
                      USE.NAMES = FALSE)
                  } else if(ncoresg1 && par_lvl == 3L) {
                    v <- mapply(p, p_affected, p1, p2, SIMPLIFY = FALSE, USE.NAMES = FALSE)
                  } else {
                    v <- mapply(p, p_affected, p1, p2, p3, SIMPLIFY = FALSE, USE.NAMES = FALSE)
                  }
                }
              } else {
                if(write_disk) {
                  v <- substring(list.files(write_dir, paste0("^", u0_n, "_p\\d{", wp, "}[.]")), wu2)
                  if(ncoresg1 && par_lvl < 3L) {
                    p <- function(V) {
                      file.copy(paste0(write_dir, "/", u0_n, V), paste0(write_dir, "/", u_n, V))
                      return(NULL)
                    }
                    if(nfork) {
                      parallel::parLapply(cl, v, p)
                    } else {
                      parallel::mclapply(v, p, mc.silent = TRUE, mc.cores = ncores)
                    }
                  } else {
                    file.copy(paste0(write_dir, "/", u0_n, v), paste0(write_dir, "/", u_n, v))
                  }
                  v <- NULL
                } else {
                  if(ncoresg1 && par_lvl == 3L) {
                    v <- lapply(p1, function(p) p[, c("g", "x", "y")])
                  } else {
                    v <- mapply(function(P1, P3) terra::vect(as.matrix(P1[, c("g", "x", "y")]), type = "line", atts = P3, crs = r_crs), p1, p3,
                      SIMPLIFY = FALSE, USE.NAMES = FALSE)
                  }
                }
              }
              return(v)
            }
            if(update_rst_list) {
              if(ncoresg1 && par_lvl == 3L) {
                if(write_disk) {
                  if(nfork) {
                    parallel::clusterMap(cl, paths, update_rst, 1:length(update_rst), USE.NAMES = FALSE)
                  } else {
                    parallel::mcmapply(paths, update_rst, 1:length(update_rst), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores)
                  }
                  paths <- NULL
                } else {
                  if(nfork) {
                    paths <- c(list(lapply(p1, function(p) p[, c("g", "x", "y")])), parallel::parLapply(cl, update_rst, paths))
                    parallel::stopCluster(cl)
                    cl <- NULL
                  } else {
                    paths <- c(list(lapply(p1, function(p) p[, c("g", "x", "y")])), parallel::mclapply(update_rst, paths, mc.silent = TRUE,
                      mc.cores = ncores))
                  }
                  rm(p1, p2, update_rst)
                  paths <- lapply(paths, function(p) mapply(function(P, P3) terra::vect(as.matrix(P), type = "line", atts = P3, crs = r_crs), p, p3,
                    SIMPLIFY = FALSE, USE.NAMES = FALSE))
                }
              } else {
                if(write_disk) {
                  mapply(paths, update_rst, 1:length(update_rst), SIMPLIFY = FALSE, USE.NAMES = FALSE)
                  paths <- NULL
                } else {
                  paths <- c(list(mapply(function(P1, P3) terra::vect(as.matrix(P1[, c("g", "x", "y")]), type = "line", atts = P3, crs = r_crs), p1, p3,
                    SIMPLIFY = FALSE, USE.NAMES = FALSE)), lapply(update_rst, paths))
                }
              }
            } else {
              if(write_disk) {
                paths <- paths(update_rst, 1L)
              } else {
                paths <- list(mapply(function(P1, P3) terra::vect(as.matrix(P1[, c("g", "x", "y")]), type = "line", atts = P3, crs = r_crs), p1, p3,
                  SIMPLIFY = FALSE, USE.NAMES = FALSE), paths(update_rst))
              }
            }
          } else {
            if(write_disk) {
              terra::writeVector(terra::vect(as.matrix(p1[, c("g", "x", "y")]), type = "line", atts = p3, crs = r_crs), paste0(write_dir, "/", u0_n, ".",
                file_type))
              u0_e <- substring(list.files(write_dir, paste0("^", u0_n, "[.]")), wu + 3L)
              u0_n <- paste0(write_dir, "/", u0_n, ".")
            }
            paths <- function(u, U = NULL) {
              p_affected <- unique(p1[.(u), "g", nomatch = NULL, on = "cls"][["g"]])
              if(length(p_affected) > 0L) {
                p_c <- data.table::data.table(origin_c = p2[p_affected, "origin_c"], destination_c = p2[p_affected, "destination_c"], pa = p_affected)
                if(any(p_c[["origin_c"]] %in% u)) {
                  report_points_ust(p_c[["origin_c"]], u, TRUE)
                }
                if(any(p_c[["destination_c"]] %in% u)) {
                  report_points_ust(p_c[["destination_c"]], u, FALSE)
                }
                v <- data.table::data.table(c_n_c = 1:NROW(crd))[-u,]
                data.table::setkey(v, c_n_c)
                p_c[, c("origin_c", "destination_c") := list(v[p_c[, "origin_c"], nomatch = NULL, which = TRUE, on = c(c_n_c = "origin_c")],
                  v[p_c[, "destination_c"], nomatch = NULL, which = TRUE, on = c(c_n_c = "destination_c")])]
                origin_n <- is.integer(p_c[["origin_c"]])
                p_c <- split(p_c, by = "origin_c", keep.by = FALSE)
                v <- crd[-u, c("x", "y")]
                rst_u <- igraph::delete_vertices(rst, u)
                if(ncoresg1 && par_lvl == 1L) {
                  V <- function(dp, oc) {
                    ps <- igraph::shortest_paths(rst_u, oc, dp[["destination_c"]], output = "vpath", algorithm = "dijkstra")$vpath
                    if(min(lengths(ps), na.rm = TRUE) == 0L) {
                      if(is.null(U)) {
                        stop("Not all points are connected when updating rst with update_rst")
                      } else {
                        stop("Not all points are connected when updating rst with element ", U, " of update_rst")
                      }
                    }
                    return(data.table::rbindlist(mapply(function(pd, pa) v[pd,][, g := pa], ps, dp[["pa"]], SIMPLIFY = FALSE, USE.NAMES = FALSE),
                      use.names = FALSE))
                  }
                  l_p_c <- length(p_c)
                  if(l_p_c < ncores) {
                    if(origin_n) {
                      origin_n <- as.integer(names(p_c))
                    } else {
                      origin_n <- as.numeric(names(p_c))
                    }
                    if(l_p_c > 1L) {
                      if(length(p_affected) - l_p_c + 1 > .Machine$integer.max) {
                        origin_m <- which.max(sapply(p_c, NROW, USE.NAMES = FALSE))
                      } else {
                        origin_m <- which.max(vapply(p_c, NROW, integer(1L), USE.NAMES = FALSE))
                      }
                      n_m <- NROW(p_c[[origin_m]])
                      l_p_c <- min(c(ncores - l_p_c + 1, n_m))
                      p_c <- c(p_c[-origin_m], split(p_c[[origin_m]], cut(1:n_m, l_p_c, labels = FALSE)))
                      origin_n <- c(origin_n[-origin_m], rep.int(origin_n[origin_m], l_p_c))
                      rm(l_p_c, origin_m, n_m)
                      if(nfork) {
                        v <- rbind(p1[!.(p_affected), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::clusterMap(cl, V, p_c, origin_n,
                          USE.NAMES = FALSE), use.names = FALSE)[, c("g", "x", "y")])
                      } else {
                        v <- rbind(p1[!.(p_affected), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::mcmapply(V, p_c, origin_n,
                          SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores), use.names = FALSE)[, c("g", "x", "y")])
                      }
                    } else {
                      rm(l_p_c)
                      p_c <- split(p_c[[1L]], cut(1:NROW(p_c[[1L]]), ncores, labels = FALSE))
                      if(nfork) {
                        v <- rbind(p1[!.(p_affected), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::parLapply(cl, p_c, V, origin_n),
                          use.names = FALSE)[, c("g", "x", "y")])
                      } else {
                        v <- rbind(p1[!.(p_affected), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::mclapply(p_c, V, origin_n,
                          mc.silent = TRUE, mc.cores = ncores), use.names = FALSE)[, c("g", "x", "y")])
                      }
                    }
                  } else {
                    if(nfork) {
                      v <- rbind(p1[!.(p_affected), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::clusterMap(cl, V, p_c,
                        ifelse(origin_n, as.integer(names(p_c)), as.numeric(names(p_c))), USE.NAMES = FALSE, .scheduling = ifelse(l_p_c == ncores,
                        "static", "dynamic")), use.names = FALSE)[, c("g", "x", "y")])
                    } else {
                      v <- rbind(p1[!.(p_affected), c("g", "x", "y"), on = "g"], data.table::rbindlist(parallel::mcmapply(V, p_c, ifelse(origin_n,
                        as.integer(names(p_c)), as.numeric(names(p_c))), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = l_p_c == ncores,
                        mc.silent = TRUE, mc.cores = ncores), use.names = FALSE)[, c("g", "x", "y")])
                    }
                  }
                } else {
                  v <- rbind(p1[!.(p_affected), c("g", "x", "y"), on = "g"], data.table::rbindlist(mapply(function(dp, oc) {
                    ps <- igraph::shortest_paths(rst_u, oc, dp[["destination_c"]], output = "vpath", algorithm = "dijkstra")$vpath
                    if(min(lengths(ps), na.rm = TRUE) == 0L) {
                      if(is.null(U)) {
                        stop("Not all points are connected when updating rst with update_rst")
                      } else {
                        stop("Not all points are connected when updating rst with element ", U, " of update_rst")
                      }
                    }
                    return(data.table::rbindlist(mapply(function(pd, pa) P[pd,][, g := pa], ps, dp[["pa"]], SIMPLIFY = FALSE, USE.NAMES = FALSE),
                      use.names = FALSE))
                  }, p_c, ifelse(origin_n, as.integer(names(p_c)), as.numeric(names(p_c))), SIMPLIFY = FALSE, USE.NAMES = FALSE),
                    use.names = FALSE)[, c("g", "x", "y")])
                }
                data.table::setorder(v, g)
                if(write_disk) {
                  terra::writeVector(terra::vect(as.matrix(v), type = "line", atts = p3, crs = r_crs), paste0(write_dir, "/u", formatC(U, flag = "0",
                    width = wu), ".", file_type))
                }
              } else {
                if(write_disk) {
                  file.copy(paste0(u0_n, u0_e), paste0(write_dir, "/u", formatC(U, flag = "0", width = wu), ".", u0_e))
                } else {
                  v <- p1[, c("g", "x", "y")]
                }
              }
              if(write_disk) {
                v <- NULL
              } else if(!(ncoresg1 && par_lvl == 3L)) {
                v <- terra::vect(as.matrix(v), type = "line", atts = p3, crs = r_crs)
              }
              return(v)
            }
            if(update_rst_list) {
              if(ncoresg1 && par_lvl == 3L) {
                if(write_disk) {
                  if(nfork) {
                    parallel::clusterMap(cl, paths, update_rst, 1:length(update_rst), USE.NAMES = FALSE)
                  } else {
                    parallel::mcmapply(paths, update_rst, 1:length(update_rst), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores)
                  }
                  paths <- NULL
                } else {
                  if(nfork) {
                    paths <- c(list(p1[, c("g", "x", "y")]), parallel::parLapply(cl, update_rst, paths))
                    parallel::stopCluster(cl)
                    cl <- NULL
                  } else {
                    paths <- c(list(p1[, c("g", "x", "y")]), parallel::mclapply(update_rst, paths, mc.silent = TRUE, mc.cores = ncores))
                  }
                  rm(p1, p2, update_rst)
                  paths <- lapply(paths, function(p) terra::vect(as.matrix(p), type = "line", atts = p3, crs = r_crs))
                }
              } else {
                if(write_disk) {
                  mapply(paths, update_rst, 1:length(update_rst), SIMPLIFY = FALSE, USE.NAMES = FALSE)
                  paths <- NULL
                } else {
                  paths <- c(list(terra::vect(as.matrix(p1[, c("g", "x", "y")]), type = "line", atts = p3, crs = r_crs)), lapply(update_rst, paths))
                }
              }
            } else {
              if(write_disk) {
                paths <- paths(update_rst, 1L)
              } else {
                paths <- list(terra::vect(as.matrix(p1[, c("g", "x", "y")]), type = "line", atts = p3, crs = r_crs), paths(update_rst))
              }
            }
          }
        # Distances output
        } else {
          if(p_list) {
            if(write_disk) {
              if(file_type_rds) {
                mapply(function(P1, P) saveRDS(P1[, c("origin", "destination", "distance")], paste0(write_dir, "/", u0_n, "_p", formatC(P, flag = "0",
                  width = wp), ".rds")), p1, 1:p_l, SIMPLIFY = FALSE, USE.NAMES = FALSE)
              } else {
                mapply(function(P1, P) data.table::fwrite(P1[, c("origin", "destination", "distance")], paste0(write_dir, "/", u0_n, "_p", formatC(P,
                  flag = "0", width = wp), ".csv")), p1, 1:p_l, SIMPLIFY = FALSE, USE.NAMES = FALSE)
              }
            }
            crd <- NROW(crd)
            paths <- function(u, U = NULL) {
              p_affected <- lapply(p2, function(p) unique(p[.(u), "g", nomatch = NULL, on = "cls"][["g"]]))
              if(write_disk) u_n <- paste0("/u", formatC(U, flag = "0", width = wu))
              if(any(lengths(p_affected) > 0L)) {
                rst_u <- igraph::delete_vertices(rst, u)
                v <- data.table::data.table(c_n_c = 1:crd)[-u,]
                p <- function(P_AFFECTED, P1, I = NULL) {
                  if(length(P_AFFECTED) > 0L) {
                    origin_c <- P1[P_AFFECTED, "origin_c"][["origin_c"]]
                    destination_c <- P1[P_AFFECTED, "destination_c"][["destination_c"]]
                    if(any(origin_c %in% u)) {
                      report_points_ust(origin_c, u, TRUE)
                    }
                    if(any(destination_c %in% u)) {
                      report_points_ust(destination_c, u, FALSE)
                    }
                    origin_c <- v[.(origin_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                    destination_c <- v[.(destination_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                    P <- data.table::copy(P1[, c("origin", "destination", "distance")])
                    if(ncoresg1 && par_lvl == 1L) {
                      p_c <- data.table::data.table(origin_c = P1[P_AFFECTED, "origin_c"][["origin_c"]], destination_c = P1[P_AFFECTED,
                        "destination_c"][["destination_c"]])
                      if(any(p_c[["origin_c"]] %in% u)) {
                        report_points_ust(p_c[["origin_c"]], u, TRUE)
                      }
                      if(any(p_c[["destination_c"]] %in% u)) {
                        report_points_ust(p_c[["destination_c"]], u, FALSE)
                      }
                      p_c[, c("origin_c", "destination_c", "p_affected") := list(v[p_c[, "origin_c"], nomatch = NULL, which = TRUE,
                        on = c(c_n_c = "origin_c")], v[p_c[, "destination_c"], nomatch = NULL, which = TRUE, on = c(c_n_c = "destination_c")],
                        P_AFFECTED)]
                      P <- data.table::copy(p1[, c("origin", "destination", "distance")])
                      n_o <- data.table::uniqueN(p_c[["origin_c"]])
                      n_d <- data.table::uniqueN(p_c[["destination_c"]])
                      l_o <- n_o > n_d
                      if(l_o) {
                        mv <- "origin_c"
                        sv <- "destination_c"
                      } else {
                        mv <- "destination_c"
                        sv <- "origin_c"
                        n_o <- n_d
                      }
                      rm(n_d)
                      data.table::setkeyv(p_c, mv)
                      P_AFFECTED <- p_c[["p_affected"]]
                      p_c[, p_affected := NULL]
                      V <- function(origin_c, destination_c, eq = FALSE) {
                        destination_u <- unique(destination_c)
                        if(eq) {
                          origin_u <- origin_c
                          origin_c <- match(destination_c, destination_u)
                        } else {
                          origin_u <- unique(origin_c)
                          origin_c <- (match(destination_c, destination_u) - 1L) * length(origin_u) + match(origin_c, origin_u)
                        }
                        return(igraph::distances(rst_u, origin_u, destination_u, mode = "out", algorithm = "dijkstra")[origin_c])
                      }
                      if(n_o == ncores) {
                        o_i <- is.integer(p_c[[mv]])
                        p_c <- lapply(split(p_c, by = mv, keep.by = FALSE), `[[`, sv)
                        rm(mv, sv, n_o)
                        if(nfork) {
                          if(l_o) {
                            p_c <- do.call(c, parallel::clusterMap(cl, V, ifelse(o_i, as.integer(names(p_c)), as.numeric(names(p_c))), p_c,
                              MoreArgs = list(eq = TRUE), USE.NAMES = FALSE))
                          } else {
                            p_c <- do.call(c, parallel::clusterMap(cl, V, p_c, ifelse(o_i, as.integer(names(p_c)), as.numeric(names(p_c))),
                              MoreArgs = list(eq = TRUE), USE.NAMES = FALSE))
                          }
                        } else {
                          if(l_o) {
                            p_c <- do.call(c, parallel::mcmapply(V, ifelse(o_i, as.integer(names(p_c)), as.numeric(names(p_c))), p_c,
                              MoreArgs = list(eq = TRUE), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores))
                          } else {
                            p_c <- do.call(c, parallel::mcmapply(V, p_c, ifelse(o_i, as.integer(names(p_c)), as.numeric(names(p_c))),
                              MoreArgs = list(eq = TRUE), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores))
                          }
                        }
                        rm(o_i, l_o)
                      } else {
                        if(n_o < ncores) {
                          l_o <- which.max(p_c[, .(n_r = .N), by = mv][["n_r"]])
                          p_c <- split(p_c, by = mv)
                          p_c <- append(p_c[-l_o], split(p_c[[l_o]], cut(1:NROW(p_c[[l_o]]), ncores - n_o + 1, labels = FALSE)), l_o - 1)
                        } else {
                          p_c <- split(p_c, rep.int(1:ncores, define_ranges(p_c[, .(n_r = .N), by = mv][["n_r"]],
                            floor(length(P_AFFECTED) / ncores), ncores)))
                        }
                        rm(mv, sv, n_o, l_o)
                        if(nfork) {
                          p_c <- do.call(c, parallel::clusterMap(cl, V, lapply(p_c, `[[`, "origin_c"), lapply(p_c, `[[`, "destination_c"),
                            USE.NAMES = FALSE))
                        } else {
                          p_c <- do.call(c, parallel::mcmapply(V, lapply(p_c, `[[`, "origin_c"), lapply(p_c, `[[`, "destination_c"), SIMPLIFY = FALSE,
                            USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores))
                        }
                      }
                      P[P_AFFECTED, distance := p_c]
                      if(unconnected_error && is.infinite(max(p_c, na.rm = TRUE))) {
                        P <- P[P[P_AFFECTED[which.max(p_c)], "origin"], nomatch = NULL, on = "origin"]
                        report_points_unc(P[1L, "origin"][["origin"]], as.integer(is.finite(P[["distance"]])), dest_specified = dest_specified,
                          d = P[["destination"]], u = U)
                      }
                      rm(V, p_c)
                    } else {
                      origin_c <- P1[P_AFFECTED, "origin_c"][["origin_c"]]
                      destination_c <- P1[P_AFFECTED, "destination_c"][["destination_c"]]
                      if(any(origin_c %in% u)) {
                        report_points_ust(origin_c, u, TRUE)
                      }
                      if(any(destination_c %in% u)) {
                        report_points_ust(destination_c, u, FALSE)
                      }
                      origin_c <- v[.(origin_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                      destination_c <- v[.(destination_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                      P <- data.table::copy(P1[, c("origin", "destination", "distance")])
                      origin_u <- unique(origin_c)
                      destination_u <- unique(destination_c)
                      origin_c <- (match(destination_c, destination_u) - 1L) * length(origin_u) + match(origin_c, origin_u)
                      rm(destination_c)
                      origin_c <- igraph::distances(igraph::delete_vertices(rst, u), origin_u, destination_u, mode = "out",
                        algorithm = "dijkstra")[origin_c]
                      P[P_AFFECTED, distance := origin_c]
                      if(unconnected_error && is.infinite(max(origin_c, na.rm = TRUE))) {
                        P <- P[P[P_AFFECTED[which.max(origin_c)], "origin"], nomatch = NULL, on = "origin"]
                        report_points_unc(P[1L, "origin"][["origin"]], as.integer(is.finite(P[["distance"]])), dest_specified = dest_specified,
                          d = P[["destination"]], u = U)
                      }
                      rm(origin_u, destination_u, origin_c)
                    }
                    if(write_disk) {
                      if(file_type_rds) {
                        saveRDS(P, paste0(write_dir, u_n, "_p", formatC(I, flag = "0", width = wp), ".rds"))
                      } else {
                        data.table::fwrite(P, paste0(write_dir, u_n, "_p", formatC(I, flag = "0", width = wp), ".csv"))
                      }
                      P <- NULL
                    }
                  } else {
                    if(write_disk) {
                      if(file_type_rds) {
                        P <- "rds"
                      } else {
                        P <- "csv"
                      }
                      P <- paste0("_p", formatC(I, flag = "0", width = wp), ".", P)
                      file.copy(paste0(write_dir, "/", u0_n, P), paste0(write_dir, u_n, P))
                      P <- NULL
                    } else {
                      P <- P1[, c("origin", "destination", "distance")]
                    }
                  }
                  return(P)
                }
                if(ncoresg1 && par_lvl == 2L) {
                  if(write_disk) {
                    if(nfork) {
                      parallel::clusterMap(cl, p, p_affected, p1, 1:p_l, USE.NAMES = FALSE, .scheduling = "dynamic")
                    } else {
                      parallel::mcmapply(p, p_affected, p1, 1:p_l, SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE, mc.silent = TRUE,
                        mc.cores = ncores)
                    }
                  } else {
                    if(nfork) {
                      v <- parallel::clusterMap(cl, p, p_affected, p1, USE.NAMES = FALSE, .scheduling = "dynamic")
                    } else {
                      v <- parallel::mcmapply(p, p_affected, p1, SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE, mc.silent = TRUE,
                        mc.cores = ncores)
                    }
                  }
                } else {
                  if(write_disk) {
                    mapply(p, p_affected, p1, 1:p_l, SIMPLIFY = FALSE, USE.NAMES = FALSE)
                  } else {
                    v <- mapply(p, p_affected, p1, SIMPLIFY = FALSE, USE.NAMES = FALSE)
                  }
                }
              } else {
                if(write_disk) {
                  if(file_type_rds) {
                    v <- "rds"
                  } else {
                    v <- "csv"
                  }
                  v <- paste0("_p", formatC(1:p_l, flag = "0", width = wp), ".", v)
                  if(ncoresg1 && par_lvl < 3L) {
                    p <- function(V) {
                      file.copy(paste0(write_dir, "/", u0_n, V), paste0(write_dir, u_n, V))
                      return(NULL)
                    }
                    if(nfork) {
                      parallel::parLapply(cl, v, p)
                    } else {
                      parallel::mclapply(v, p, mc.silent = TRUE, mc.cores = ncores)
                    }
                  } else {
                    file.copy(paste0(write_dir, "/", u0_n, v), paste0(write_dir, u_n, v))
                  }
                } else {
                  v <- lapply(p1, function(p) p[, c("origin", "destination", "distance")])
                }
              }
              if(write_disk) v <- NULL
              return(v)
            }
            if(update_rst_list) {
              if(ncoresg1 && par_lvl == 3L) {
                if(write_disk) {
                  if(nfork) {
                    parallel::clusterMap(cl, paths, update_rst, 1:length(update_rst), USE.NAMES = FALSE)
                  } else {
                    parallel::mcmapply(paths, update_rst, 1:length(update_rst), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores)
                  }
                } else {
                  if(nfork) {
                    paths <- c(list(lapply(p1, function(p) p[, c("origin", "destination", "distance")])), parallel::parLapply(cl, update_rst, paths))
                  } else {
                    paths <- c(list(lapply(p1, function(p) p[, c("origin", "destination", "distance")])), parallel::mclapply(update_rst, paths,
                      mc.silent = TRUE, mc.cores = ncores))
                  }
                }
              } else {
                if(write_disk) {
                  mapply(paths, update_rst, 1:length(update_rst), SIMPLIFY = FALSE, USE.NAMES = FALSE)
                } else {
                  paths <- c(list(lapply(p1, function(p) p[, c("origin", "destination", "distance")])), lapply(update_rst, paths))
                }
              }
              if(write_disk) paths <- NULL
            } else {
              if(write_disk) {
                paths <- paths(update_rst, 1L)
              } else {
                paths <- list(lapply(p1, function(p) p[, c("origin", "destination", "distance")]), paths(update_rst))
              }
            }
          } else {
            if(write_disk) {
              if(file_type_rds) {
                saveRDS(p1[, c("origin", "destination", "distance")], paste0(write_dir, "/", u0_n, ".rds"))
              } else {
                data.table::fwrite(p1[, c("origin", "destination", "distance")], paste0(write_dir, "/", u0_n, ".csv"))
              }
            }
            crd <- NROW(crd)
            paths <- function(u, U = NULL) {
              p_affected <- unique(p2[.(u), "g", nomatch = NULL, on = "cls"][["g"]])
              if(length(p_affected) > 0L) {
                if(ncoresg1 && par_lvl == 1L) {
                  p_c <- data.table::data.table(origin_c = p1[p_affected, "origin_c"][["origin_c"]], destination_c = p1[p_affected,
                    "destination_c"][["destination_c"]])
                  if(any(p_c[["origin_c"]] %in% u)) {
                    report_points_ust(p_c[["origin_c"]], u, TRUE)
                  }
                  if(any(p_c[["destination_c"]] %in% u)) {
                    report_points_ust(p_c[["destination_c"]], u, FALSE)
                  }
                  v <- data.table::data.table(c_n_c = 1:crd)[-u,]
                  p_c[, c("origin_c", "destination_c", "p_affected") := list(v[p_c[, "origin_c"], nomatch = NULL, which = TRUE,
                    on = c(c_n_c = "origin_c")], v[p_c[, "destination_c"], nomatch = NULL, which = TRUE, on = c(c_n_c = "destination_c")], p_affected)]
                  v <- data.table::copy(p1[, c("origin", "destination", "distance")])
                  n_o <- data.table::uniqueN(p_c[["origin_c"]])
                  n_d <- data.table::uniqueN(p_c[["destination_c"]])
                  l_o <- n_o > n_d
                  if(l_o) {
                    mv <- "origin_c"
                    sv <- "destination_c"
                  } else {
                    mv <- "destination_c"
                    sv <- "origin_c"
                    n_o <- n_d
                  }
                  rm(n_d)
                  data.table::setkeyv(p_c, mv)
                  p_affected <- p_c[["p_affected"]]
                  p_c[, p_affected := NULL]
                  rst_u <- igraph::delete_vertices(rst, u)
                  p <- function(origin_c, destination_c, eq = FALSE) {
                    destination_u <- unique(destination_c)
                    if(eq) {
                      origin_u <- origin_c
                      origin_c <- match(destination_c, destination_u)
                    } else {
                      origin_u <- unique(origin_c)
                      origin_c <- (match(destination_c, destination_u) - 1L) * length(origin_u) + match(origin_c, origin_u)
                    }
                    return(igraph::distances(rst_u, origin_u, destination_u, mode = "out", algorithm = "dijkstra")[origin_c])
                  }
                  if(n_o == ncores) {
                    o_i <- is.integer(p_c[[mv]])
                    p_c <- lapply(split(p_c, by = mv, keep.by = FALSE), `[[`, sv)
                    rm(mv, sv, n_o)
                    if(nfork) {
                      if(l_o) {
                        p_c <- do.call(c, parallel::clusterMap(cl, p, ifelse(o_i, as.integer(names(p_c)), as.numeric(names(p_c))), p_c,
                          MoreArgs = list(eq = TRUE), USE.NAMES = FALSE))
                      } else {
                        p_c <- do.call(c, parallel::clusterMap(cl, p, p_c, ifelse(o_i, as.integer(names(p_c)), as.numeric(names(p_c))),
                          MoreArgs = list(eq = TRUE), USE.NAMES = FALSE))
                      }
                    } else {
                      if(l_o) {
                        p_c <- do.call(c, parallel::mcmapply(p, ifelse(o_i, as.integer(names(p_c)), as.numeric(names(p_c))), p_c,
                          MoreArgs = list(eq = TRUE), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores))
                      } else {
                        p_c <- do.call(c, parallel::mcmapply(p, p_c, ifelse(o_i, as.integer(names(p_c)), as.numeric(names(p_c))),
                          MoreArgs = list(eq = TRUE), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores))
                      }
                    }
                    rm(o_i, l_o)
                  } else {
                    if(n_o < ncores) {
                      l_o <- which.max(p_c[, .(n_r = .N), by = mv][["n_r"]])
                      p_c <- split(p_c, by = mv)
                      p_c <- append(p_c[-l_o], split(p_c[[l_o]], cut(1:NROW(p_c[[l_o]]), ncores - n_o + 1, labels = FALSE)), l_o - 1)
                    } else {
                      p_c <- split(p_c, rep.int(1:ncores, define_ranges(p_c[, .(n_r = .N), by = mv][["n_r"]],
                        floor(length(p_affected) / ncores), ncores)))
                    }
                    rm(mv, sv, n_o, l_o)
                    if(nfork) {
                      p_c <- do.call(c, parallel::clusterMap(cl, p, lapply(p_c, `[[`, "origin_c"), lapply(p_c, `[[`, "destination_c"),
                        USE.NAMES = FALSE))
                    } else {
                      p_c <- do.call(c, parallel::mcmapply(p, lapply(p_c, `[[`, "origin_c"), lapply(p_c, `[[`, "destination_c"), SIMPLIFY = FALSE,
                        USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores))
                    }
                  }
                  rm(rst_u, p)
                  v[p_affected, distance := p_c]
                  if(unconnected_error && is.infinite(max(p_c, na.rm = TRUE))) {
                    v <- v[v[p_affected[which.max(p_c)], "origin"], nomatch = NULL, on = "origin"]
                    report_points_unc(v[1L, "origin"][["origin"]], as.integer(is.finite(v[["distance"]])), dest_specified = dest_specified,
                      d = v[["destination"]], u = U)
                  }
                  rm(p_affected, p_c)
                } else {
                  origin_c <- p1[p_affected, "origin_c"][["origin_c"]]
                  destination_c <- p1[p_affected, "destination_c"][["destination_c"]]
                  if(any(origin_c %in% u)) {
                    report_points_ust(origin_c, u, TRUE)
                  }
                  if(any(destination_c %in% u)) {
                    report_points_ust(destination_c, u, FALSE)
                  }
                  v <- data.table::data.table(c_n_c = 1:crd)[-u,]
                  origin_c <- v[.(origin_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                  destination_c <- v[.(destination_c), nomatch = NULL, which = TRUE, on = "c_n_c"]
                  v <- data.table::copy(p1[, c("origin", "destination", "distance")])
                  origin_u <- unique(origin_c)
                  destination_u <- unique(destination_c)
                  origin_c <- (match(destination_c, destination_u) - 1L) * length(origin_u) + match(origin_c, origin_u)
                  rm(destination_c)
                  origin_c <- igraph::distances(igraph::delete_vertices(rst, u), origin_u, destination_u, mode = "out", algorithm = "dijkstra")[origin_c]
                  v[p_affected, distance := origin_c]
                  if(unconnected_error && is.infinite(max(origin_c, na.rm = TRUE))) {
                    v <- v[v[p_affected[which.max(origin_c)], "origin"], nomatch = NULL, on = "origin"]
                    report_points_unc(v[1L, "origin"][["origin"]], as.integer(is.finite(v[["distance"]])), dest_specified = dest_specified,
                      d = v[["destination"]], u = U)
                  }
                  rm(p_affected, origin_u, destination_u, origin_c)
                }
                if(write_disk) {
                  if(file_type_rds) {
                    saveRDS(v, paste0(write_dir, "/u", formatC(U, flag = "0", width = wu), ".rds"))
                  } else {
                    data.table::fwrite(v, paste0(write_dir, "/u", formatC(U, flag = "0", width = wu), ".csv"))
                  }
                  v <- NULL
                }
              } else {
                if(write_disk) {
                  if(file_type_rds) {
                    v <- "rds"
                  } else {
                    v <- "csv"
                  }
                  file.copy(paste0(write_dir, "/", u0_n, ".", v), paste0(write_dir, "/u", formatC(U, flag = "0", width = wu), ".", v))
                  v <- NULL
                } else {
                  v <- p1[, c("origin", "destination", "distance")]
                }
              }
              return(v)
            }
            if(update_rst_list) {
              if(ncoresg1 && par_lvl == 3L) {
                if(write_disk) {
                  if(nfork) {
                    parallel::clusterMap(cl, paths, update_rst, 1:length(update_rst), USE.NAMES = FALSE)
                  } else {
                    parallel::mcmapply(paths, update_rst, 1:length(update_rst), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores)
                  }
                } else {
                  if(nfork) {
                    paths <- c(list(p1[, c("origin", "destination", "distance")]), parallel::parLapply(cl, update_rst, paths))
                  } else {
                    paths <- c(list(p1[, c("origin", "destination", "distance")]), parallel::mclapply(paths, update_rst, mc.silent = TRUE,
                      mc.cores = ncores))
                  }
                }
              } else {
                if(write_disk) {
                  mapply(paths, update_rst, 1:length(update_rst), SIMPLIFY = FALSE, USE.NAMES = FALSE)
                } else {
                  paths <- c(list(p1[, c("origin", "destination", "distance")]), lapply(update_rst, paths))
                }
              }
              if(write_disk) paths <- NULL
            } else {
              if(write_disk) {
                paths <- paths(update_rst, 1L)
              } else {
                paths <- list(p1[, c("origin", "destination", "distance")], paths(update_rst))
              }
            }
          }
        }
      # Not copy
      } else {
        if(write_disk) {
          if(p_list) {
            p <- "_p"
          } else {
            p <- ""
            wp <- NULL
          }
        }
        if(ncoresg1 && par_lvl == 3L) {
          if(update_rst_list) {
            update_rst <- c(list(NULL), update_rst)
            if(write_disk) {
              paths <- function(V, u) {
                if(V == 0L) {
                  compute_spaths1(rst, crd[, c("x", "y")], origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified,
                    origin_list, dest_list, r_crs, output_lines, pairwise, FALSE, WRITE_DISK = TRUE, file_type = file_type,
                    file_type_rds = file_type_rds, nm = paste0(write_dir, "/u", formatC(0L, flag = "0", width = wu), p), wp = wp,
                    unconnected_error = unconnected_error, tr_directed = tr_directed)
                } else {
                  compute_spaths1(igraph::delete_vertices(rst, u), crd[-u, c("x", "y")], origins, destinations, dest_specified, origin_nms_specified,
                    destination_nms_specified, origin_list, dest_list, r_crs, output_lines, pairwise, FALSE, WRITE_DISK = TRUE, file_type = file_type,
                    file_type_rds = file_type_rds, nm = paste0(write_dir, "/u", formatC(V, flag = "0", width = wu), p), wp = wp,
                    unconnected_error = unconnected_error, tr_directed = tr_directed)
                }
                return(NULL)
              }
              if(nfork) {
                parallel::clusterMap(cl, paths, 0:length(update_rst), update_rst, USE.NAMES = FALSE)
              } else {
                parallel::mcmapply(paths, 0:length(update_rst), update_rst, SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE, mc.cores = ncores)
              }
              paths <- NULL
            } else {
              paths <- function(V, u) {
                if(V == 0L) {
                  v <- compute_spaths1(rst, crd[, c("x", "y")], origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified,
                    origin_list, dest_list, r_crs, output_lines, pairwise, FALSE, unconnected_error = unconnected_error, tr_directed = tr_directed)
                } else {
                  v <- compute_spaths1(igraph::delete_vertices(rst, u), crd[-u, c("x", "y")], origins, destinations, dest_specified, origin_nms_specified,
                    destination_nms_specified, origin_list, dest_list, r_crs, output_lines, pairwise, FALSE, unconnected_error = unconnected_error,
                    tr_directed = tr_directed)
                }
                return(v)
              }
              if(nfork) {
                paths <- parallel::clusterMap(cl, paths, 0:length(update_rst), update_rst, USE.NAMES = FALSE)
                parallel::stopCluster(cl)
                cl <- NULL
              } else {
                paths <- parallel::mcmapply(paths, 0:length(update_rst), update_rst, SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.silent = TRUE,
                  mc.cores = ncores)
              }
            }
          } else {
            if(write_disk) {
              paths <- function(V) {
                if(V == 0L) {
                  compute_spaths1(rst, crd[, c("x", "y")], origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified,
                    origin_list, dest_list, r_crs, output_lines, pairwise, FALSE, WRITE_DISK = TRUE, file_type = file_type, file_type_rds = file_type_rds,
                    nm = paste0(write_dir, "/u", formatC(0L, flag = "0", width = wu), p), wp = wp, unconnected_error = unconnected_error,
                    tr_directed = tr_directed)
                } else {
                  compute_spaths1(igraph::delete_vertices(rst, update_rst), crd[-update_rst, c("x", "y")], origins, destinations, dest_specified,
                    origin_nms_specified, destination_nms_specified, origin_list, dest_list, r_crs, output_lines, pairwise, FALSE, WRITE_DISK = TRUE,
                    file_type = file_type, file_type_rds = file_type_rds, nm = paste0(write_dir, "/u", formatC(1L, flag = "0", width = wu), p), wp = wp,
                    unconnected_error = unconnected_error, tr_directed = tr_directed)
                }
                return(NULL)
              }
              if(nfork) {
                parallel::parLapply(cl, 0:1, paths)
              } else {
                parallel::mclapply(0:1, paths, mc.silent = TRUE, mc.cores = ncores)
              }
              paths <- NULL
            } else {
              paths <- function(V) {
                if(V == 0L) {
                  v <- compute_spaths1(rst, crd[, c("x", "y")], origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified,
                    origin_list, dest_list, r_crs, output_lines, pairwise, FALSE, unconnected_error = unconnected_error, tr_directed = tr_directed)
                } else {
                  v <- compute_spaths1(igraph::delete_vertices(rst, update_rst), crd[-update_rst, c("x", "y")], origins, destinations, dest_specified,
                    origin_nms_specified, destination_nms_specified, origin_list, dest_list, r_crs, output_lines, pairwise, FALSE,
                    unconnected_error = unconnected_error, tr_directed = tr_directed)
                }
                return(v)
              }
              if(nfork) {
                paths <- parallel::parLapply(cl, 0:1, paths)
                parallel::stopCluster(cl)
                cl <- NULL
              } else {
                paths <- parallel::mclapply(0:1, paths, mc.silent = TRUE, mc.cores = ncores)
              }
            }
          }
          rm(rst, crd, update_rst)
          if(output_lines && !write_disk) {
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
            if(write_disk) {
              lapply(0:length(update_rst), function(V) {
                if(V == 0L) {
                  compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list,
                    dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, nfork, cl, FALSE, TRUE, file_type, file_type_rds,
                    paste0(write_dir, "/u", formatC(0L, flag = "0", width = wu), p), wp, unconnected_error, tr_directed)
                } else {
                  compute_spaths1(igraph::delete_vertices(rst, update_rst[[V]]), crd[-update_rst[[V]],], origins, destinations, dest_specified,
                    origin_nms_specified, destination_nms_specified, origin_list, dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl,
                    nfork, cl, FALSE, TRUE, file_type, file_type_rds, paste0(write_dir, "/u", formatC(V, flag = "0", width = wu), p), wp,
                    unconnected_error, tr_directed)
                }
                return(NULL)
              })
              paths <- NULL
            } else {
              paths <- lapply(0:length(update_rst), function(V) {
                if(V == 0L) {
                  v <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified,
                    origin_list, dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, nfork, cl, unconnected_error = unconnected_error,
                    tr_directed = tr_directed)
                } else {
                  v <- compute_spaths1(igraph::delete_vertices(rst, update_rst[[V]]), crd[-update_rst[[V]],], origins, destinations, dest_specified,
                    origin_nms_specified, destination_nms_specified, origin_list, dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl,
                    nfork, cl, unconnected_error = unconnected_error, tr_directed = tr_directed)
                }
                return(v)
              })
            }
          } else {
            if(write_disk) {
              compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list, dest_list,
                r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, nfork, cl, FALSE, TRUE, file_type, file_type_rds, paste0(write_dir, "/u",
                formatC(0L, flag = "0", width = wu), p), wp, unconnected_error, tr_directed)
            } else {
              paths <- list(compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified,
                origin_list, dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, nfork, cl, unconnected_error = unconnected_error,
                tr_directed = tr_directed))
            }
            rst <- igraph::delete_vertices(rst, update_rst)
            rm(update_rst)
            crd <- crd[-update_rst,]
            if(write_disk) {
              paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list,
                dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, nfork, cl, FALSE, TRUE, file_type, file_type_rds, paste0(write_dir,
                "/u", formatC(1L, flag = "0", width = wu), p), wp, unconnected_error, tr_directed)
            } else {
              paths <- c(paths, compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified,
                origin_list, dest_list, r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, nfork, cl, unconnected_error = unconnected_error,
                tr_directed = tr_directed))
            }
          }
        }
      }
    # update_rst does not mask any non-NA pixels
    } else {
      warning("The update_rst geometries do not mask any non-NA pixels. Thus, the results are the same as those from the unmodified rst.")
      if(write_disk) {
        if(update_rst_list) {
          update_rst <- length(update_rst)
        } else {
          update_rst <- 1L
        }
        if(p_list) {
          p <- "_p"
        } else {
          p <- ""
          wp <- NULL
        }
      } else {
        if(update_rst_list) {
          update_rst <- length(update_rst) + 1L
        } else {
          update_rst <- 2L
        }
      }
      if(!ncoresg1 || (par_lvl == 3L && copy && rst_par_lvl == 3L)) {
        if(write_disk) {
          paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list,
            dest_list, r_crs, output_lines, pairwise, FALSE, WRITE_DISK = TRUE, file_type = file_type, file_type_rds = file_type_rds,
            nm = paste0(write_dir, "/", u0_n, p), wp = wp, unconnected_error = unconnected_error, tr_directed = tr_directed)
        } else {
          paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list,
            dest_list, r_crs, output_lines, pairwise, FALSE, unconnected_error = unconnected_error, tr_directed = tr_directed)
        }
      } else {
        if(par_lvl == 3L && copy) par_lvl <- rst_par_lvl %% 3L
        if(write_disk) {
          paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list,
            dest_list, r_crs, output_lines, pairwise, TRUE, ncores, par_lvl, nfork, cl, WRITE_DISK = TRUE, file_type = file_type,
            file_type_rds = file_type_rds, nm = paste0(write_dir, "/", u0_n, p), wp = wp, unconnected_error = unconnected_error,
            tr_directed = tr_directed)
        } else {
          paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list,
            dest_list, r_crs, output_lines, pairwise, TRUE, ncores, par_lvl, nfork, cl, unconnected_error = unconnected_error, tr_directed = tr_directed)
        }
      }
      rm(rst, crd, origins, destinations)
      if(write_disk) {
        if(output_lines) {
          p <- substring(list.files(write_dir, paste0("^", u0_n, p)), wu + 2L)
        } else {
          p <- substring(list.files(write_dir, paste0("^", u0_n, p, "\\d*[.]", file_type, "$")), wu + 2L)
        }
        if(ncoresg1) {
          P <- function(V) {
            file.copy(paste0(write_dir, "/", u0_n, p), paste0(write_dir, "/u", formatC(V, flag = "0", width = wu), p))
            return(NULL)
          }
          if(nfork) {
            parallel::parLapply(cl, 1:update_rst, P)
          } else {
            parallel::mclapply(1:update_rst, P, mc.silent = TRUE, mc.cores = ncores)
          }
        } else {
          lapply(1:update_rst, function(V) file.copy(paste0(write_dir, "/", u0_n, p), paste0(write_dir, "/u", formatC(V, flag = "0", width = wu), p)))
        }
      } else {
        paths <- replicate(update_rst, paths, FALSE)
      }
    }
  # update_rst not specified
  } else {
    if(write_disk) {
      if(p_list) {
        paths <- ""
      } else {
        paths <- "results"
        wp <- NULL
      }
      paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list, dest_list,
        r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, nfork, cl, WRITE_DISK = TRUE, file_type = file_type, file_type_rds = file_type_rds,
        nm = paths, wp = wp, unconnected_error = unconnected_error, tr_directed = tr_directed)
    } else {
      paths <- compute_spaths1(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list, dest_list,
        r_crs, output_lines, pairwise, ncoresg1, ncores, par_lvl, nfork, cl, unconnected_error = unconnected_error, tr_directed = tr_directed)
    }
  }
  if(ncoresg1 && nfork && !is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- NULL
  }
  return(paths)
}

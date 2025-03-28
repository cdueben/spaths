# Functions called by both shortest_paths

# Convert origins and destinations
convert_points <- function(v, rst, r_crs, nms_column, nms_specified, pairwise, rst_terra, rst_list, n_grids, rst_xmin, rst_xmax, rst_ymin, rst_ymax,
  rst_xres, rst_yres, rst_ncol, o = TRUE) {
  if(rst_terra) {
    if(all(class(v) != "SpatVector")) v <- terra::vect(v)
    if(!terra::is.points(v)) v <- terra::centroids(v)
    if(terra::crs(v) != r_crs) v <- terra::project(v, r_crs)
    if(nms_specified) nms <- names(v)
  } else {
    if(is.matrix(v)) {
      v_df <- FALSE
      nms <- colnames(v)
    } else if(is.data.frame(v)) {
      v_df <- TRUE
      nms <- names(v)
    } else {
      stop(data.table::fifelse(o, "origins", "destinations"), " must be a SpatVector (terra), sf (sf), Spatial* (sp), vector, matrix, or data.frame object")
    }
    if(!all(c("x", "y") %chin% nms)) stop(data.table::fifelse(o, "origins", "destinations"), " must contain columns named x and y")
    if(v_df) {
      v_xmin <- min(v[["x"]])
      v_xmax <- max(v[["x"]])
      v_ymin <- min(v[["y"]])
      v_ymax <- max(v[["y"]])
    } else {
      v_xmin <- min(v[, "x"])
      v_xmax <- max(v[, "x"])
      v_ymin <- min(v[, "y"])
      v_ymax <- max(v[, "y"])
    }
    if(!is.finite(v_xmin) || !is.finite(v_xmax) || !is.finite(v_ymin) || !is.finite(v_ymax)) {
      stop(data.table::fifelse(o, "origins", "destinations"), " must only contain finite, non-NA values")
    }
    if(v_xmin < rst_xmin || v_xmax > rst_xmax || v_ymin < rst_ymin || v_ymax > rst_ymax) {
      stop(data.table::fifelse(o, "origins", "destinations"), " fall outside rst")
    }
    if(nms_specified && nms_column %chin% c("x", "y")) stop(data.table::fifelse(o, "origin", "destination"), "_names must not be x or y") 
  }
  if(nms_specified) {
    if(!(nms_column %chin% nms)) {
      stop(data.table::fifelse(o, "origin", "destination"), "_names must either be NULL or the name of a column in ", data.table::fifelse(o, "origins",
        "destinations"))
    }
    if(rst_terra) {
      nms <- unlist(terra::values(v[, nms_column]), use.names = FALSE)
    } else if(v_df) {
      nms <- v[[nms_column]]
    } else {
      nms <- v[, nms_column]
    }
  }
  if(rst_terra) {
    v_cells <- terra::extract(rst, v, cells = TRUE, ID = FALSE)
    na_points <- !all(stats::complete.cases(v_cells[, 1:n_grids]))
  } else {
    if(v_df) {
      v_cells <- rst_ncol * as.integer((rst_ymax - v[["y"]]) / rst_yres) + as.integer((v[["x"]] - rst_xmin) / rst_xres) + 1L
    } else {
      v_cells <- rst_ncol * as.integer((rst_ymax - v[, "y"]) / rst_yres) + as.integer((v[, "x"] - rst_xmin) / rst_xres) + 1L
    }
    if(rst_list) {
      na_points <- !all(stats::complete.cases(rst[v_cells,]))
    } else {
      na_points <- anyNA(rst[v_cells])
    }
  }
  if(na_points) {
    if(rst_terra) {
      v <- which(!stats::complete.cases(v_cells[, 1:n_grids]))
    } else if(rst_list) {
      v <- which(!stats::complete.cases(rst[v_cells,]))
    } else {
      v <- which(is.na(rst[v_cells]))
    }
    v_length <- length(v)
    if(v_length > 1L) {
      if(v_length > 2L) {
        v <- paste0(paste0(v[1:(v_length - 1L)], collapse = ", "), ", and ", v[v_length])
      } else {
        v <- paste0(v, collapse = " and ")
      }
      v <- paste0(data.table::fifelse(o, "Origins", "Destinations"), " ", v, " located on NA cells")
    } else {
      v <- paste0(data.table::fifelse(o, "Origin", "Destination"), " ", v, " located on NA cell")
    }
    stop(v)
  }
  if(rst_terra) {
    v_cells <- v_cells$cell
  }
  if(!pairwise && anyDuplicated(v_cells) != 0L) {
    w <- which(duplicated(v_cells))
    w_length <- length(w)
    if(w_length > 1L) {
      if(w_length > 2L) {
        w <- paste0(paste0(w[1:(w_length - 1L)], collapse = ", "), ", and ", w[w_length])
      } else {
        w <- paste0(w, collapse = " and ")
      }
      w <- paste0(data.table::fifelse(o, "Origins", "Destinations"), " ", w, " are duplicates, in that they fall")
    } else {
      w <- paste0(data.table::fifelse(o, "Origin", "Destination"), " ", w, " is a duplicate, in that it falls")
    }
    w <- paste0(w, " into the same cell as another ", data.table::fifelse(o, "origin", "destination"))
    warning(w, ". This causes the function to run longer than necessary.")
  }
  if(nms_specified) {
    v <- list(cls = v_cells, nms = nms)
  } else {
    v <- list(cls = v_cells)
  }
  return(v)
}

# Get the maximum number of neighbors stored in an adjacency list in R
get_max_neighbors <- function(n_cells, n_cells_na, queen, rst_ncol, rst_nrow, global) {
  # Remove cross-boundary edges
  max_neighbors <- n_cells * data.table::fifelse(queen, 8, 4) - rst_ncol * data.table::fifelse(queen, 6, 2)
  if(!global) {
    max_neighbors <- max_neighbors - rst_nrow * data.table::fifelse(queen, 6, 2)
  }
  
  # Remove NA cells
  if(queen) {
    if(global) {
      if(n_cells_na > rst_ncol) {
        max_neighbors <- max_neighbors - 5 * rst_ncol - 8 * (n_cells_na - rst_ncol)
      } else {
        max_neighbors <- max_neighbors - 5 * n_cells_na
      }
    } else {
      if(n_cells_na > 4L) {
        max_neighbors <- max_neighbors - 12
        n_cells_na <- n_cells_na - 4L                                          # corner cells
        n_fringe_cells <- (rst_ncol + rst_nrow - 4L) * 2L
        if(n_cells_na > n_fringe_cells) {
          max_neighbors <- max_neighbors - n_fringe_cells * 5
          n_cells_na <- n_cells_na - n_fringe_cells                            # fringe cells
          max_neighbors <- max_neighbors - n_cells_na * 8
        } else {
          max_neighbors <- max_neighbors - n_cells_na * 5
        }
      } else {
        max_neighbors <- max_neighbors - 3 * n_cells_na
      }
    }
  } else {
    if(n_cells_na > 4L) {
      max_neighbors <- max_neighbors - 8
      n_cells_na <- n_cells_na - 4L                                            # corner cells
      n_fringe_cells <- (rst_ncol + rst_nrow - 4L) * 2L
      if(n_cells_na > n_fringe_cells) {
        max_neighbors <- max_neighbors - n_fringe_cells * 3
        n_cells_na <- n_cells_na - n_fringe_cells                              # fringe cells
        max_neighbors <- max_neighbors - n_cells_na * 4
      } else {
        max_neighbors <- max_neighbors - n_cells_na * 3
      }
    } else {
      max_neighbors <- max_neighbors - 2 * n_cells_na
    }
  }
  
  return(max_neighbors)
}

# Check whether the C++ functions implement a conversion for the layer classes
layer_class_check <- function(v_vars_classes) {
  if(!all(v_vars_classes %chin% c("integer", "numeric", "character", "logical"))) {
    v_vars_classes <- paste0(unique(v_vars_classes[!(v_vars_classes %chin% c("integer", "numeric", "character", "logical"))]), collapse = ", ")
    stop("For graphs with potentially more than ", .Machine$integer.max, " edges, tr_fun's v parameters only accept integer, numeric, character, and ",
      "logical, not ", v_vars_classes, ", layers")
  }
}

# Avoid R CMD check note
utils::globalVariables(c(".", "cell_numbers", "cell_numbers_shifted", "connected", "layer", "path", "starts", "from", "to", "xres", "yres"))

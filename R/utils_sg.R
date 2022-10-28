# Functions called by spaths_general
# Capitalized object names differing from the ones in the spaths_general function indicate that these objects might differ from the spaths_general
# counterparts, e.g. by being subsets

# Check points' location and list the corresponding grid cell numbers
convert_points_g <- function(v, rst, rst_list, nr, nc, xres, yres, xmin, ymin, nms, nms_specified, o = TRUE) {
  if(!any(class(v) %chin% c("matrix", "data.frame"))) stop(ifelse(o, "origins", "destinations"), " must be a matrix or data frame, or a list of them")
  if(nms_specified) {
    if(!(nms %chin% names(v))) {
      stop(ifelse(o, "origin", "destination"), "_names must either be NULL or the name of a column in ", ifelse(o, "origins", "destinations"))
    }
    nms <- unlist(v[, nms], use.names = FALSE)
  }
  if(data.table::is.data.table(v)) {
    v_col <- check_position(v[, 1L][[1L]], xmin, xres, nc, o)
    v_row <- check_position(v[, 2L][[1L]], ymin, yres, nr, o)
  } else {
    v_col <- check_position(v[, 1L], xmin, xres, nc, o)
    v_row <- check_position(v[, 2L], ymin, yres, nr, o)
  }
  v <- nr * v_col + v_row + 1L
  if(rst_list) {
    if(any(vapply(rst, function(r) anyNA(r[v]), logical(1L), USE.NAMES = FALSE))) {
      report_points(unique(do.call(c, lapply(rst, function(r) which(is.na(r[v]))))), o)
    }
  } else if(anyNA(rst[v])) {
    report_points(which(is.na(rst[v])), o)
  }
  v <- nc * v_row + v_col + 1L
  if(nms_specified) v <- data.table::data.table(cls = v, nms = nms)
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
compute_spaths_g <- function(ORIGINS, rst, crd, dest_specified, origin_nms_specified, output_lines, pairwise, NCORESG1, ncores, cl, nvect,
  DESTINATIONS = NULL, destination_nms_specified = TRUE) {
  os_length <- NROW(ORIGINS)
  if(origin_nms_specified) {
    on <- ORIGINS[["nms"]]
    ORIGINS <- ORIGINS[["cls"]]
  } else {
    on <- 1:os_length
  }
  # Shortest paths when destinations are specified
  if(dest_specified) {
    ds_length <- NROW(DESTINATIONS)
    if(destination_nms_specified) {
      dn <- DESTINATIONS[["nms"]]
      DESTINATIONS <- DESTINATIONS[["cls"]]
    } else {
      dn <- 1:ds_length
    }
    # Lines output
    if(output_lines) {
      if(pairwise) {
        if(NCORESG1) {
          p <- function(O) {
            s <- lapply(O, function(o) igraph::shortest_paths(rst, ORIGINS[o], DESTINATIONS[o], output = "vpath", algorithm = "dijkstra")$vpath[[1L]])
            return(data.table::data.table(cls = do.call(c, s), g = rep.int(O, lengths(s))))
          }
          p <- terra::vect(as.matrix(data.table::rbindlist(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)),
            p), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line", atts = data.frame(origin = on, destination = dn))
        } else if(nvect) {
          p <- list(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
            return(data.table::rbindlist(lapply(igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS[O], output = "vpath", algorithm = "dijkstra")$vpath,
              function(D) crd[D,][, g := O]), use.names = FALSE))
          }), use.names = FALSE)[, c("g", "x", "y")]), data.frame(origin = on, destination = dn))
        } else {
          p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
            return(data.table::rbindlist(lapply(igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS[O], output = "vpath", algorithm = "dijkstra")$vpath,
              function(D) crd[D,][, g := O]), use.names = FALSE))
          }), use.names = FALSE)[, c("g", "x", "y")]), type = "line", atts = data.frame(origin = on, destination = dn))
        }
      } else {
        if(NCORESG1) {
          p <- function(O) {
            return(data.table::rbindlist(lapply(O, function(o) {
              s <- igraph::shortest_paths(rst, ORIGINS[o], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
              o1 <- (o - 1L) * ds_length
              return(data.table::data.table(cls = do.call(c, s), g = o1 + rep.int(1:ds_length, lengths(s))))
            }), use.names = FALSE))
          }
          p <- terra::vect(as.matrix(data.table::rbindlist(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)),
            p), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line", atts = data.frame(origin = rep.int(on, rep.int(ds_length,
            os_length)), destination = rep.int(dn, os_length)))
        } else if(nvect) {
          p <- list(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
            s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
            o1 <- (O - 1L) * ds_length
            return(data.table::rbindlist(lapply(1:ds_length, function(D) crd[s[[D]],][, g := o1 + D]), use.names = FALSE))
          }), use.names = FALSE)[, c("g", "x", "y")]), data.frame(origin = rep.int(on, rep.int(ds_length, os_length)), destination = rep.int(dn, os_length)))
        } else {
          p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
            s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
            o1 <- (O - 1L) * ds_length
            return(data.table::rbindlist(lapply(1:ds_length, function(D) crd[s[[D]],][, g := o1 + D]), use.names = FALSE))
          }), use.names = FALSE)[, c("g", "x", "y")]), type = "line", atts = data.frame(origin = rep.int(on, rep.int(ds_length, os_length)),
            destination = rep.int(dn, os_length)))
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
            p <- data.table::rbindlist(parallel::parLapply(cl, split(ORIGINS, cut(seq_along(ORIGINS), ncores, labels = FALSE)), p), use.names = FALSE)
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
            return(data.table::data.table(cls = do.call(c, s), g = i + rep.int(o[["V2"]], lengths(s))))
          }), use.names = FALSE))
        }
        p <- terra::vect(as.matrix(data.table::rbindlist(parallel::parLapply(cl, split(data.table::as.data.table(t(utils::combn(1:os_length, 2L))),
          cut(1:(os_length * o1 / 2L), ncores, labels = FALSE)), p), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line",
          atts = stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")))
      } else if(nvect) {
        p <- list(as.matrix(data.table::rbindlist(lapply(1:o1, function(O) {
          npO <- os_length - O
          i <- sum(o1:npO) - npO
          s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
          return(data.table::rbindlist(lapply(1:npO, function(D) crd[s[[D]],][, g := i + D])))
        }), use.names = FALSE)[, c("g", "x", "y")]), stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")))
      } else {
        p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:o1, function(O) {
          npO <- os_length - O
          i <- sum(o1:npO) - npO
          s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
          return(data.table::rbindlist(lapply(1:npO, function(D) crd[s[[D]],][, g := i + D])))
        }), use.names = FALSE)[, c("g", "x", "y")]), type = "line", atts = stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin",
          "destination")))
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
          }), use.names = FALSE))
        }
        p <- data.table::rbindlist(parallel::parLapply(cl, split(data.table::as.data.table(t(utils::combn(1:os_length, 2L))),
          cut(1:(os_length * (os_length - 1L) / 2L), ncores, labels = FALSE)), p), use.names = FALSE)
      } else {
        p <- data.table::rbindlist(lapply(1:(os_length - 1L), function(O) {
          s <- data.table::as.data.table(igraph::distances(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], mode = "out", algorithm = "dijkstra"), na.rm = FALSE)
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
        }), use.names = FALSE)
      }
    }
  }
  return(p)
}

# Functions called by spaths_earth
# Capitalized object names differing from the ones in the spaths_earth function indicate that these objects might differ from the spaths_earth
# counterparts, e.g. by being subsets

# Convert origins and destinations
convert_points <- function(v, rst, r_crs, nms, nms_specified, o = TRUE) {
  if(all(class(v) != "SpatVector")) v <- terra::vect(v)
  if(!terra::is.points(v)) v <- terra::centroids(v)
  if(terra::crs(v) != r_crs) v <- terra::project(v, r_crs)
  if(nms_specified) {
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
  if(nms_specified) {
    v <- data.table::data.table(cls = v$cell, nms = nms)
  } else {
    v <- v$cell
  }
  return(v)
}

# Compute distances between neighboring cells
compute_dists <- function(rst, CRD, dist_comp_terra, round_dist, contiguity, yr, xr, nr, ym, lonlat, ncoresg1, ncores, cl) {
  if(dist_comp_terra) {
    if(ncoresg1) {
      d <- function(RST) {
        return(terra::distance(as.matrix(CRD[RST[["from"]],]), as.matrix(CRD[RST[["to"]],]), lonlat = lonlat, pairwise = TRUE))
      }
      d <- do.call(c, parallel::parLapply(cl, split(rst, cut(1:NROW(rst), ncores, labels = FALSE)), d))
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
compute_spaths1 <- function(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list, dest_list, r_crs,
  output_lines, pairwise, NCORESG1, ncores = NULL, par_lvl = NULL, cl = NULL, copy = FALSE) {
  if(origin_list) {
    if(dest_specified) {
      if(dest_list) {
        if(NCORESG1 && par_lvl == "points_lists") {
          paths <- function(O) {
            return(compute_spaths2(origins[[O]], rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, FALSE, NULL, NULL, TRUE, copy,
              destinations[[O]], destination_nms_specified))
          }
          paths <- parallel::parLapplyLB(cl, 1:length(origins), paths)
          if(output_lines && !copy) paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]], crs = r_crs))
        } else {
          paths <- mapply(compute_spaths2, ORIGINS = origins, DESTINATIONS = destinations, MoreArgs = list(rst = rst, crd = crd, dest_specified = TRUE,
            origin_nms_specified = origin_nms_specified, r_crs = r_crs, output_lines = output_lines, pairwise = pairwise, NCORESG1 = NCORESG1,
            ncores = ncores, cl = cl, nvect = is.null(par_lvl), copy = copy, destination_nms_specified = destination_nms_specified), SIMPLIFY = FALSE,
            USE.NAMES = FALSE)
        }
      } else {
        if(NCORESG1 && par_lvl == "points_lists") {
          paths <- function(O) {
            return(compute_spaths2(O, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, FALSE, NULL, NULL, TRUE, copy, destinations,
              destination_nms_specified))
          }
          paths <- parallel::parLapplyLB(cl, origins, paths)
          if(output_lines && !copy) paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]], crs = r_crs))
        } else {
          paths <- lapply(origins, compute_spaths2, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, cl,
            is.null(par_lvl), copy, destinations, destination_nms_specified)
        }
      }
    } else {
      paths <- lapply(origins, compute_spaths2, rst, crd, FALSE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, cl, is.null(par_lvl),
        copy)
    }
  } else {
    if(dest_specified) {
      if(dest_list) {
        if(NCORESG1 && par_lvl == "points_lists") {
          paths <- function(D) {
            return(compute_spaths2(origins, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, FALSE, NULL, NULL, TRUE, copy,
              destinations[[D]], destination_nms_specified))
          }
          paths <- parallel::parLapplyLB(cl, 1:length(destinations), paths)
          if(output_lines && !copy) paths <- lapply(paths, function(D) terra::vect(D[[1L]], type = "line", atts = D[[2L]], crs = r_crs))
        } else {
          paths <- lapply(destinations, function(d) compute_spaths2(origins, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1,
            ncores, cl, is.null(par_lvl), copy, d, destination_nms_specified))
        }
      } else {
        paths <- compute_spaths2(origins, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, cl, is.null(par_lvl), copy,
          destinations, destination_nms_specified)
      }
    } else {
      paths <- compute_spaths2(origins, rst, crd, FALSE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, cl, is.null(par_lvl), copy)
    }
  }
  return(paths)
}

# Compute shortest paths
compute_spaths2 <- function(ORIGINS, rst, crd, dest_specified, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, cl, nvect, copy,
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
        if(copy) {
          if(NCORESG1) {
            p <- function(O) {
              s <- lapply(O, function(o) igraph::shortest_paths(rst, ORIGINS[o], DESTINATIONS[o], output = "vpath", algorithm = "dijkstra")$vpath[[1L]])
              pl <- lengths(s)
              if(any(pl == 0L)) report_points_unc(O, pl, TRUE)
              return(data.table::data.table(g = rep.int(O, pl), cls = do.call(c, s)))
            }
            p <- list(data.table::rbindlist(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)), p),
              use.names = FALSE)[, c("x", "y") := crd[cls,]][, c("g", "x", "y", "cls")], data.frame(origin_c = ORIGINS, destination_c = DESTINATIONS),
              data.frame(origin = on, destination = dn))
          } else {
            p <- list(data.table::rbindlist(lapply(1:os_length, function(O) {
              D <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS[O], output = "vpath", algorithm = "dijkstra")$vpath[[1L]]
              if(length(D) == 0L) stop("Origin ", O, " is not connected to its respective destination")
              return(crd[D,][, c("g", "cls") := list(O, D)])
            }), use.names = FALSE)[, c("g", "x", "y", "cls")], data.frame(origin_c = ORIGINS, destination_c = DESTINATIONS), data.frame(origin = on,
              destination = dn))
          }
        } else {
          if(NCORESG1) {
            p <- function(O) {
              s <- lapply(O, function(o) igraph::shortest_paths(rst, ORIGINS[o], DESTINATIONS[o], output = "vpath", algorithm = "dijkstra")$vpath[[1L]])
              pl <- lengths(s)
              if(any(pl == 0L)) report_points_unc(O, pl, TRUE)
              return(data.table::data.table(g = rep.int(O, pl), cls = do.call(c, s)))
            }
            p <- terra::vect(as.matrix(data.table::rbindlist(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)), p),
              use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line", atts = data.frame(origin = on, destination = dn), crs = r_crs)
          } else if(nvect) {
            p <- list(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
              D <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS[O], output = "vpath", algorithm = "dijkstra")$vpath[[1L]]
              if(length(D) == 0L) stop("Origin ", O, " is not connected to its respective destination")
              return(crd[D,][, g := O])
            }), use.names = FALSE)[, c("g", "x", "y")]), data.frame(origin = on, destination = dn))
          } else {
            p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
              D <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS[O], output = "vpath", algorithm = "dijkstra")$vpath[[1L]]
              if(length(D) == 0L) stop("Origin ", O, " is not connected to its respective destination")
              return(crd[D,][, g := O])
            }), use.names = FALSE)[, c("g", "x", "y")]), type = "line", atts = data.frame(origin = on, destination = dn), crs = r_crs)
          }
        }
      } else {
        if(copy) {
          if(NCORESG1) {
            p <- function(O) {
              return(data.table::rbindlist(lapply(O, function(o) {
                s <- igraph::shortest_paths(rst, ORIGINS[o], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
                o1 <- (o - 1L) * ds_length
                pl <- lengths(s)
                if(any(pl == 0L)) report_points_unc(o, pl)
                return(data.table::data.table(g = o1 + rep.int(1:ds_length, pl), cls = do.call(c, s)))
              }), use.names = FALSE))
            }
            p <- list(data.table::rbindlist(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)), p),
              use.names = FALSE)[, c("x", "y") := crd[cls,]][, c("g", "x", "y", "cls")], data.frame(origin_c = rep.int(ORIGINS, rep.int(ds_length,
              os_length)), destination_c = rep.int(DESTINATIONS, os_length)), data.frame(origin = rep.int(on, rep.int(ds_length, os_length)),
              destination = rep.int(dn, os_length)))
          } else {
              p <- list(data.table::rbindlist(lapply(1:os_length, function(O) {
                s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
                o1 <- (O - 1L) * ds_length
                if(any(lengths(s) == 0L)) report_points_unc(O, lengths(s))
                return(data.table::rbindlist(lapply(1:ds_length, function(D) crd[s[[D]],][, c("g", "cls") := list(o1 + D, s[[D]])]), use.names = FALSE))
              }), use.names = FALSE)[, c("g", "x", "y", "cls")], data.frame(origin_c = rep.int(ORIGINS, rep.int(ds_length, os_length)),
                destination_c = rep.int(DESTINATIONS, os_length)), data.frame(origin = rep.int(on, rep.int(ds_length, os_length)),
                destination = rep.int(dn, os_length)))
          }
        } else {
          if(NCORESG1) {
            p <- function(O) {
              return(data.table::rbindlist(lapply(O, function(o) {
                s <- igraph::shortest_paths(rst, ORIGINS[o], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
                o1 <- (o - 1L) * ds_length
                pl <- lengths(s)
                if(any(pl == 0L)) report_points_unc(o, pl)
                return(data.table::data.table(g = o1 + rep.int(1:ds_length, pl), cls = do.call(c, s)))
              }), use.names = FALSE))
            }
            p <- terra::vect(as.matrix(data.table::rbindlist(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)),
              p), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line", atts = data.frame(origin = rep.int(on, rep.int(ds_length,
              os_length)), destination = rep.int(dn, os_length)), crs = r_crs)
          } else if(nvect) {
            p <- list(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
              s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
              o1 <- (O - 1L) * ds_length
              if(any(lengths(s) == 0L)) report_points_unc(O, lengths(s))
              return(data.table::rbindlist(lapply(1:ds_length, function(D) crd[s[[D]],][, g := o1 + D]), use.names = FALSE))
            }), use.names = FALSE)[, c("g", "x", "y")]), data.frame(origin = rep.int(on, rep.int(ds_length, os_length)), destination = rep.int(dn,
              os_length)))
          } else {
            p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
              s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
              o1 <- (O - 1L) * ds_length
              if(any(lengths(s) == 0L)) report_points_unc(O, lengths(s))
              return(data.table::rbindlist(lapply(1:ds_length, function(D) crd[s[[D]],][, g := o1 + D]), use.names = FALSE))
            }), use.names = FALSE)[, c("g", "x", "y")]), type = "line", atts = data.frame(origin = rep.int(on, rep.int(ds_length, os_length)),
              destination = rep.int(dn, os_length)), crs = r_crs)
          }
        }
      }
    # Distances output
    } else {
      if(pairwise) {
        if(copy) {
          if(NCORESG1) {
            p <- function(O) {
              s <- lapply(O, function(o) {
                S <- igraph::shortest_paths(rst, ORIGINS[o], DESTINATIONS[o], output = "both", algorithm = "dijkstra")
                if(length(S$vpath[[1L]]) == 0L) stop("Origin ", o, " is not connected to its respective destination")
                return(list(sum(igraph::edge_attr(rst, "weight", S$epath[[1L]])), S$vpath[[1L]]))
              })
              return(list(vapply(s, `[[`, 1L, numeric(1L), USE.NAMES = FALSE), data.table::data.table(g = rep.int(O, sapply(s, function(o) length(o[[2L]]),
                USE.NAMES = FALSE)), cls = do.call(c, lapply(s, `[[`, 2L)))))
            }
            p <- parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)), p)
            p <- list(data.table::data.table(origin = on, origin_c = ORIGINS, destination = dn, destination_c = DESTINATIONS, distance = do.call(c,
              lapply(p, `[[`, 1L))), data.table::rbindlist(lapply(p, `[[`, 2L), use.names = FALSE))
          } else {
            p <- lapply(1:os_length, function(O) {
              S <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS[O], output = "both", algorithm = "dijkstra")
              if(length(S$vpath[[1L]]) == 0L) stop("Origin ", O, " is not connected to its respective destination")
              return(list(sum(igraph::edge_attr(rst, "weight", S$epath[[1L]])), S$vpath[[1L]]))
            })
            p <- list(data.table::data.table(origin = on, origin_c = ORIGINS, destination = dn, destination_c = DESTINATIONS, distance = vapply(p,
              `[[`, 1L, numeric(1L), USE.NAMES = FALSE)), data.table::data.table(g = rep.int(O, sapply(p, function(o) length(o[[2L]]), USE.NAMES = FALSE)),
              cls = do.call(c, lapply(p, `[[`, 2L))))
          }
        } else {
          if(NCORESG1) {
            p <- function(O) {
              s <- vapply(O, function(o) igraph::distances(rst, ORIGINS[o], DESTINATIONS[o], mode = "out", algorithm = "dijkstra"), numeric(1L),
                USE.NAMES = FALSE)
              if(any(is.infinite(s))) report_points_unc(O, as.integer(is.finite(s)), TRUE)
              return(s)
            }
            p <- data.table::data.table(origin = on, destination = dn, distance = do.call(c, parallel::parLapply(cl, split(1:os_length, cut(1:os_length,
              ncores, labels = FALSE)), p)))
          } else {
            p <- data.table::data.table(origin = on, destination = dn, distance = vapply(1:os_length, function(O) igraph::distances(rst, ORIGINS[O],
              DESTINATIONS[O], mode = "out", algorithm = "dijkstra"), numeric(1L), USE.NAMES = FALSE))
            if(any(is.infinite(p[["distance"]]))) report_points_unc(1:os_length, as.integer(is.finite(p[["distance"]])), TRUE)
          }
        }
      } else {
        if(copy) {
          if(NCORESG1) {
            if(os_length >= ncores || os_length > ds_length) {
              p <- function(O) {
                P <- lapply(O, function(o) {
                  s <- igraph::shortest_paths(rst, o, DESTINATIONS, output = "both", algorithm = "dijkstra")
                  if(any(lengths(s$vpath) == 0L)) report_points_unc(NULL, lengths(s$vpath), both = FALSE)
                  return(list(vapply(s$epath, function(S) sum(igraph::edge_attr(rst, "weight", S)), numeric(1L), USE.NAMES = FALSE), s$vpath))
                })
                return(list(do.call(c, lapply(P, `[[`, 1L)), do.call(c, lapply(P, `[[`, 2L))))
              }
              p <- parallel::parLapply(cl, split(ORIGINS, cut(seq_along(ORIGINS), ncores, labels = FALSE)), p)
              p <- list(data.table::data.table(origin = rep(on, each = ds_length), origin_c = rep(ORIGINS, each = ds_length), destination = rep.int(dn,
                os_length), destination_c = rep.int(DESTINATIONS, os_length), distance = do.call(c, lapply(p, `[[`, 1L))), do.call(c, lapply(P, `[[`, 2L)))
              p[[2L]] <- data.table::data.table(g = rep.int(1:(os_length * ds_length), lengths(p[[2L]])), cls = do.call(c, p[[2L]]))
            } else {
              p <- function(D) {
                P <- lapply(D, function(d) {
                  s <- igraph::shortest_paths(rst, ORIGINS, d, output = "both", algorithm = "dijkstra")
                  if(any(lengths(s$vpath) == 0L)) report_points_unc(NULL, lengths(s$vpath), O = FALSE, both = FALSE)
                  return(list(vapply(s$epath, function(S) sum(igraph::edge_attr(rst, "weight", S)), numeric(1L), USE.NAMES = FALSE), s$vpath))
                })
                return(list(do.call(c, lapply(P, `[[`, 1L)), do.call(c, lapply(P, `[[`, 2L))))
              }
              p <- parallel::parLapply(cl, split(DESTINATIONS, cut(seq_along(DESTINATIONS), ncores, labels = FALSE)), p)
              s <- seq.int(1L, by = ds_length, length.out = os_length) + rep(0:(ds_length - 1L), each = os_length)
              p <- list(data.table::data.table(origin = rep(on, each = ds_length), origin_c = rep(ORIGINS, each = ds_length), destination = rep.int(dn,
                os_length), destination_c = rep.int(DESTINATIONS, os_length), distance = do.call(c, lapply(p, `[[`, 1L))[s]), do.call(c, lapply(P, `[[`,
                2L))[[s]])
              p[[2L]] <- data.table::data.table(g = rep.int(1:(os_length * ds_length), lengths(p[[2L]])), cls = do.call(c, p[[2L]]))
            }
          } else {
            p <- lapply(1:os_length, function(O) {
              s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "both", algorithm = "dijkstra")
              if(any(lengths(s$vpath) == 0L)) report_points_unc(O, lengths(s$vpath))
              return(list(vapply(s$epath, function(S) sum(igraph::edge_attr(rst, "weight", S)), numeric(1L), USE.NAMES = FALSE), s$vpath))
            })
            p <- list(data.table::data.table(origin = rep(on, each = ds_length), origin_c = rep(ORIGINS, each = ds_length), destination = rep.int(dn,
              os_length), destination_c = rep.int(DESTINATIONS, os_length), distance = do.call(c, lapply(p, `[[`, 1L))), do.call(c, lapply(p, `[[`, 2L)))
            s <- os_length * ds_length
            if(s > 1L) {
              p[[2L]] <- data.table::data.table(g = rep.int(1:s, lengths(p[[2L]])), cls = do.call(c, p[[2L]]))
            } else {
              p[[2L]] <- data.table::data.table(g = 1L, cls = p[[2L]])
            }
          }
        } else {
          if(NCORESG1) {
            if(os_length >= ncores || os_length > ds_length) {
              p <- function(O) {
                s <- igraph::distances(rst, O, DESTINATIONS, mode = "out", algorithm = "dijkstra")
                if(any(is.infinite(s))) {
                  report_points_unc(NULL, as.integer(is.finite(s[min(which(is.infinite(s)) %% length(O)) + 1L,])), both = FALSE)
                }
                return(data.table::as.data.table(s, na.rm = FALSE))
              }
              p <- data.table::rbindlist(parallel::parLapply(cl, split(ORIGINS, cut(seq_along(ORIGINS), ncores, labels = FALSE)), p), use.names = FALSE)
            } else {
              p <- function(D) {
                s <- igraph::distances(rst, ORIGINS, D, mode = "out", algorithm = "dijkstra")
                if(any(is.infinite(s))) {
                  report_points_unc(NULL, as.integer(is.finite(s[, ceiling(min(which(is.infinite(s)) / os_length))])), O = FALSE, both = FALSE)
                }
                return(s)
              }
              p <- data.table::as.data.table(do.call(cbind, parallel::parLapply(cl, split(DESTINATIONS, cut(seq_along(DESTINATIONS), ncores, labels = FALSE)),
                p)), na.rm = FALSE)
            }
          } else {
            p <- data.table::as.data.table(igraph::distances(rst, ORIGINS, DESTINATIONS, mode = "out", algorithm = "dijkstra"), na.rm = FALSE)
            if(any(is.infinite(p))) {
              o <- min(which(is.infinite(p)) %% os_length) + 1L
              report_points_unc(o, as.integer(is.finite(unlist(p[o,], use.names = FALSE))))
            }
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
    }
  # Shortest paths when destinations are not specified
  } else {
    # Lines output
    if(output_lines) {
      o1 <- os_length - 1L
      if(copy) {
        if(NCORESG1) {
          p <- function(O) {
            return(data.table::rbindlist(lapply(split(O, by = "V1"), function(o) {
              s <- o[1L, "V1"][["V1"]]
              i <- os_length - s
              i <- sum(o1:i) - i
              s <- igraph::shortest_paths(rst, ORIGINS[s], ORIGINS[o[["V2"]]], output = "vpath", algorithm = "dijkstra")$vpath
              pl <- lengths(s)
              if(any(pl == 0L)) report_points_unc(s, pl, dest_specified = FALSE, d = o[["V2"]])
              return(data.table::data.table(g = i + rep.int(o[["V2"]], pl), cls = do.call(c, s)))
            }), use.names = FALSE))
          }
          p <- list(data.table::rbindlist(parallel::parLapply(cl, split(data.table::as.data.table(t(utils::combn(1:os_length, 2L))),
            cut(1:(os_length * o1 / 2L), ncores, labels = FALSE)), p), use.names = FALSE)[, c("x", "y") := crd[cls,]][, c("g", "x", "y", "cls")],
            stats::setNames(as.data.frame(t(utils::combn(ORIGINS, 2L))), c("origin_c", "destination_c")), stats::setNames(as.data.frame(t(utils::combn(on,
            2L))), c("origin", "destination")))
        } else {
          p <- list(data.table::rbindlist(lapply(1:o1, function(O) {
            npO <- os_length - O
            i <- sum(o1:npO) - npO
            s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
            if(any(lengths(s) == 0L)) report_points_unc(O, lengths(s), dest_specified = FALSE, d = (O + 1L):os_length)
            return(data.table::rbindlist(lapply(1:npO, function(D) crd[s[[D]],][, c("g", "cls") := list(i + D, s[[D]])]), use.names = FALSE))
          }), use.names = FALSE)[, c("g", "x", "y", "cls")], stats::setNames(as.data.frame(t(utils::combn(ORIGINS, 2L))), c("origin_c", "destination_c")),
          stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")))
        }
      } else {
        if(NCORESG1) {
          p <- function(O) {
            return(data.table::rbindlist(lapply(split(O, by = "V1"), function(o) {
              s <- o[1L, "V1"][["V1"]]
              i <- os_length - s
              i <- sum(o1:i) - i
              s <- igraph::shortest_paths(rst, ORIGINS[s], ORIGINS[o[["V2"]]], output = "vpath", algorithm = "dijkstra")$vpath
              pl <- lengths(s)
              if(any(pl == 0L)) report_points_unc(s, pl, dest_specified = FALSE, d = o[["V2"]])
              return(data.table::data.table(g = i + rep.int(o[["V2"]], pl), cls = do.call(c, s)))
            }), use.names = FALSE))
          }
          p <- terra::vect(as.matrix(data.table::rbindlist(parallel::parLapply(cl, split(data.table::as.data.table(t(utils::combn(1:os_length, 2L))),
            cut(1:(os_length * o1 / 2L), ncores, labels = FALSE)), p), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line",
            atts = stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")), crs = r_crs)
        } else if(nvect) {
          p <- list(as.matrix(data.table::rbindlist(lapply(1:o1, function(O) {
            npO <- os_length - O
            i <- sum(o1:npO) - npO
            s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
            if(any(lengths(s) == 0L)) report_points_unc(O, lengths(s), dest_specified = FALSE, d = (O + 1L):os_length)
            return(data.table::rbindlist(lapply(1:npO, function(D) crd[s[[D]],][, g := i + D]), use.names = FALSE))
          }), use.names = FALSE)[, c("g", "x", "y")]), stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")))
        } else {
          p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:o1, function(O) {
            npO <- os_length - O
            i <- sum(o1:npO) - npO
            s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
            if(any(lengths(s) == 0L)) report_points_unc(O, lengths(s), dest_specified = FALSE, d = (O + 1L):os_length)
            return(data.table::rbindlist(lapply(1:npO, function(D) crd[s[[D]],][, g := i + D]), use.names = FALSE))
          }), use.names = FALSE)[, c("g", "x", "y")]), type = "line", atts = stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin",
            "destination")), crs = r_crs)
        }
      }
    # Distances output
    } else {
      if(copy) {
        if(NCORESG1) {
          p <- function(O) {
            P <- lapply(split(O, by = "V1"), function(o) {
              g <- o[1L, "V1"][["V1"]]
              i <- os_length - g
              i <- sum(o1:i) - i
              s <- igraph::shortest_paths(rst, ORIGINS[g], ORIGINS[o[["V2"]]], output = "both", algorithm = "dijkstra")
              pl <- lengths(s$vpath)
              if(any(pl == 0L)) report_points_unc(g, pl, dest_specified = FALSE, d = o[["V2"]])
              return(list(data.table::data.table(origin = on[g], origin_c = ORIGINS[g], destination = on[o[["V2"]]],
                destination_c = ORIGINS[o[["V2"]]], distance = vapply(s$epath, function(S) sum(igraph::edge_attr(rst, "weight", S)), numeric(1L),
                USE.NAMES = FALSE)), data.table::data.table(g = i + rep.int(o[["V2"]], pl), cls = do.call(c, s$vpath))))
            })
            return(list(data.table::rbindlist(lapply(P, `[[`, 1L), use.names = FALSE), data.table::rbindlist(lapply(P, `[[`, 2L), use.names = FALSE)))
          }
          p <- data.table::rbindlist(parallel::parLapply(cl, split(data.table::as.data.table(t(utils::combn(1:os_length, 2L))),
            cut(1:(os_length * (os_length - 1L) / 2L), ncores, labels = FALSE)), p), use.names = FALSE)
          p <- list(data.table::rbindlist(lapply(p, `[[`, 1L)), data.table::rbindlist(lapply(p, `[[`, 2L)))
        } else {
          p <- lapply(1:o1, function(O) {
            npO <- os_length - O
            i <- sum(o1:npO) - npO
            D <- (O + 1L):os_length
            s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[D], output = "both", algorithm = "dijkstra")
            pl <- lengths(s$vpath)
            if(any(pl == 0L)) report_points_unc(O, pl, dest_specified = FALSE, d = D)
            return(list(data.table::data.table(origin = on[O], origin_c = ORIGINS[O], destination = on[D], destination_c = ORIGINS[D],
              distance = vapply(s$epath, function(S) sum(igraph::edge_attr(rst, "weight", S)), numeric(1L), USE.NAMES = FALSE)),
              data.table::data.table(g = i + rep.int(1:npO, pl), cls = do.call(c, s$vpath))))
          })
          p <- list(data.table::rbindlist(lapply(p, `[[`, 1L), use.names = FALSE), data.table::rbindlist(lapply(p, `[[`, 2L), use.names = FALSE))
        }
      } else {
        on_n <- is.numeric(on)
        if(on_n) on_i <- is.integer(on)
        if(NCORESG1) {
          p <- function(O) {
            return(data.table::rbindlist(lapply(split(O, by = "V1"), function(o) {
              g <- o[1L, "V1"][["V1"]]
              s <- data.table::as.data.table(igraph::distances(rst, ORIGINS[g], ORIGINS[o[["V2"]]], mode = "out", algorithm = "dijkstra"), na.rm = FALSE)
              if(any(is.infinite(s))) report_points_unc(g, as.integer(is.finite(unlist(s, use.names = FALSE))), dest_specified = FALSE, d = o[["V2"]])
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
            s <- data.table::as.data.table(igraph::distances(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], mode = "out", algorithm = "dijkstra"),
              na.rm = FALSE)
            if(any(is.infinite(s))) {
              report_points_unc(O, as.integer(is.finite(unlist(s, use.names = FALSE))), dest_specified = FALSE, d = (O + 1L):os_length)
            }
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
  }
  return(p)
}

# Report points on cells masked by update_rst
report_points_ust <- function(v, u, o, cl) {
  v <- which(v %in% u)
  v_length <- length(v)
  if(v_length > 1L) {
    if(v_length > 2L) {
      v <- paste0("cells of ", ifelse(o, "origin", "destination"), " points ", paste0(v[1:(v_length - 1L)]), ", and ", v[v_length])
    } else {
      v <- paste0("cells of ", ifelse(o, "origin", "destination"), " points ", paste0(v, collapse = " and "))
    }
  } else {
    v <- paste0("cell of ", ifelse(o, "origin", "destination"), " point ", v)
  }
  if(!is.null(cl)) parallel::stopCluster(cl)
  stop("update_rst masks the ", v)
}

# Report unconnected origins and destinations
report_points_unc <- function(o, pl, pairwise = FALSE, dest_specified = TRUE, d = NULL, O = TRUE, both = TRUE) {
  if(pairwise) {
    s <- o[pl == 0L]
    sl <- length(s)
    if(sl > 1L) {
      if(sl > 2L) {
        s <- paste0("s ", paste0(s[1:(sl - 1L)], collapse = ", "), ", and ", s[sl])
      } else {
        s <- paste0("s ", paste0(s, collapse = " and "))
      }
      s <- paste0(s, " are not connected to their respective destinations")
    } else {
      s <- paste0(" ", s, " is not connected to its respective destination")
    }
    s <- paste0("Origin", s)
  } else {
    if(dest_specified) {
      s <- which(pl == 0L)
    } else {
      s <- d[pl == 0L]
    }
    sl <- length(s)
    if(sl > 1L) {
      if(sl > 2L) {
        s <- paste0("s ", paste0(s[1:(sl - 1L)], collapse = ", "), ", and ", s[sl])
      } else {
        s <- paste0("s ", paste0(s, collapse = " and "))
      }
      if(!both) s <- paste0(s, " are")
    } else {
      s <- paste0(" ", s)
      if(!both) s <- paste0(s, " is")
    }
    if(both) {
      s <- paste0("Origin ", o, " is not connected to ", ifelse(dest_specified, "destination", "origin"), s)
    } else {
      if(O) {
        s <- paste0("Destination", s, " not connected to all origins")
      } else {
        s <- paste0("Origin", s, " not connected to all destinations")
      }
    }
  }
  stop(s)
}

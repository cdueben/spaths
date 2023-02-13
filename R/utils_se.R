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
compute_spaths1 <- function(rst, crd, origins, destinations, dest_specified, origin_nms_specified, destination_nms_specified, origin_list, dest_list,
  r_crs, output_lines, pairwise, NCORESG1, ncores = NULL, par_lvl = NULL, nfork = NULL, cl = NULL, copy = FALSE, WRITE_DISK = FALSE, file_type = NULL,
  file_type_rds = NULL, nm = NULL, wp = NULL, unconnected_error = TRUE, tr_directed = FALSE) {
  if(WRITE_DISK) {
    write_file <- function(o, P) {
      if(output_lines) {
        terra::writeVector(o, paste0(nm, formatC(P, flag = "0", width = wp), ".", file_type))
      } else if(file_type_rds) {
        saveRDS(o, paste0(nm, formatC(P, flag = "0", width = wp), ".rds"))
      } else {
        data.table::fwrite(o, paste0(nm, formatC(P, flag = "0", width = wp), ".csv"))
      }
      return(NULL)
    }
  }
  if(origin_list) {
    if(dest_specified) {
      if(dest_list) {
        if(NCORESG1 && par_lvl == 2L) {
          if(WRITE_DISK) {
            paths <- function(O, D, P) {
              return(write_file(compute_spaths2(O, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, FALSE, NULL, NULL, NULL, FALSE,
                FALSE, unconnected_error, tr_directed, D, destination_nms_specified), P))
            }
            if(nfork) {
              parallel::clusterMap(cl, paths, origins, destinations, 1:length(origins), USE.NAMES = FALSE, .scheduling = "dynamic")
            } else {
              parallel::mcmapply(paths, origins, destinations, 1:length(origins), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE,
                mc.silent = TRUE, mc.cores = ncores)
            }
          } else {
            paths <- function(O, D) {
              return(compute_spaths2(O, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, FALSE, NULL, NULL, NULL, TRUE, copy,
                unconnected_error, tr_directed, D, destination_nms_specified))
            }
            if(nfork) {
              paths <- parallel::clusterMap(cl, paths, origins, destinations, USE.NAMES = FALSE, .scheduling = "dynamic")
            } else {
              paths <- parallel::mcmapply(paths, origins, destinations, SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE, mc.silent = TRUE,
                mc.cores = ncores)
            }
            if(output_lines && !copy) paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]], crs = r_crs))
          }
        } else {
          if(WRITE_DISK) {
            mapply(function(O, D, P) write_file(compute_spaths2(O, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores,
              nfork, cl, FALSE, FALSE, unconnected_error, tr_directed, D, destination_nms_specified), P), origins, destinations, 1:length(origins),
              SIMPLIFY = FALSE, USE.NAMES = FALSE)
          } else {
            paths <- mapply(compute_spaths2, ORIGINS = origins, DESTINATIONS = destinations, MoreArgs = list(rst = rst, crd = crd, dest_specified = TRUE,
              origin_nms_specified = origin_nms_specified, r_crs = r_crs, output_lines = output_lines, pairwise = pairwise, NCORESG1 = NCORESG1,
              ncores = ncores, nfork = nfork, cl = cl, nvect = is.null(par_lvl), copy = copy, unconnected_error = unconnected_error,
              tr_directed = tr_directed, destination_nms_specified = destination_nms_specified), SIMPLIFY = FALSE, USE.NAMES = FALSE)
          }
        }
      } else {
        if(NCORESG1 && par_lvl == 2L) {
          if(WRITE_DISK) {
            paths <- function(O, P) {
              return(write_file(compute_spaths2(O, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, FALSE, NULL, NULL, NULL, FALSE,
                FALSE, unconnected_error, tr_directed, destinations, destination_nms_specified), P))
            }
            if(nfork) {
              parallel::clusterMap(cl, paths, origins, 1:length(origins), USE.NAMES = FALSE, .scheduling = "dynamic")
            } else {
              parallel::mcmapply(paths, origins, 1:length(origins), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE, mc.silent = TRUE,
                mc.cores = ncores)
            }
          } else {
            paths <- function(O) {
              return(compute_spaths2(O, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, FALSE, NULL, NULL, NULL, TRUE, copy,
                unconnected_error, tr_directed, destinations, destination_nms_specified))
            }
            if(nfork) {
              paths <- parallel::parLapplyLB(cl, origins, paths)
            } else {
              paths <- parallel::mclapply(origins, paths, mc.preschedule = FALSE, mc.silent = TRUE, mc.cores = ncores)
            }
            if(output_lines && !copy) paths <- lapply(paths, function(O) terra::vect(O[[1L]], type = "line", atts = O[[2L]], crs = r_crs))
          }
        } else {
          if(WRITE_DISK) {
            mapply(function(O, P) write_file(compute_spaths2(O, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores,
              nfork, cl, FALSE, FALSE, unconnected_error, tr_directed, destinations, destination_nms_specified), P), origins, 1:length(origins),
              SIMPLIFY = FALSE, USE.NAMES = FALSE)
          } else {
            paths <- lapply(origins, compute_spaths2, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, nfork,
              cl, is.null(par_lvl), copy, unconnected_error, tr_directed, destinations, destination_nms_specified)
          }
        }
      }
    } else {
      if(WRITE_DISK) {
        mapply(function(O, P) write_file(compute_spaths2(O, rst, crd, FALSE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, nfork,
          cl, FALSE, FALSE, unconnected_error, tr_directed), P), origins, 1:length(origins), SIMPLIFY = FALSE, USE.NAMES = FALSE)
      } else {
        paths <- lapply(origins, compute_spaths2, rst, crd, FALSE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, nfork, cl,
          is.null(par_lvl), copy, unconnected_error, tr_directed)
      }
    }
  } else {
    if(dest_specified) {
      if(dest_list) {
        if(NCORESG1 && par_lvl == 2L) {
          if(WRITE_DISK) {
            paths <- function(D, P) {
              return(write_file(compute_spaths2(origins, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, FALSE, NULL, NULL, NULL,
                FALSE, FALSE, unconnected_error, tr_directed, D, destination_nms_specified), P))
            }
            if(nfork) {
              parallel::clusterMap(cl, paths, destinations, 1:length(destinations), USE.NAMES = FALSE, .scheduling = "dynamic")
            } else {
              parallel::mcmapply(paths, destinations, 1:length(destinations), SIMPLIFY = FALSE, USE.NAMES = FALSE, mc.preschedule = FALSE,
                mc.silent = TRUE, mc.cores = ncores)
            }
          } else {
            paths <- function(D) {
              return(compute_spaths2(origins, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, FALSE, NULL, NULL, NULL, TRUE, copy,
                unconnected_error, tr_directed, D, destination_nms_specified))
            }
            if(nfork) {
              paths <- parallel::parLapplyLB(cl, destinations, paths)
            } else {
              paths <- parallel::mclapply(destinations, paths, mc.preschedule = FALSE, mc.silent = TRUE, mc.cores = ncores)
            }
            if(output_lines && !copy) paths <- lapply(paths, function(D) terra::vect(D[[1L]], type = "line", atts = D[[2L]], crs = r_crs))
          }
        } else {
          if(WRITE_DISK) {
            mapply(function(d, P) write_file(compute_spaths2(origins, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise,
              NCORESG1, ncores, nfork, cl, FALSE, FALSE, unconnected_error, tr_directed, d, destination_nms_specified), P), destinations,
              1:length(destinations), SIMPLIFY = FALSE, USE.NAMES = FALSE)
          } else {
            paths <- lapply(destinations, function(d) compute_spaths2(origins, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise,
              NCORESG1, ncores, nfork, cl, is.null(par_lvl), copy, unconnected_error, tr_directed, d, destination_nms_specified))
          }
        }
      } else {
        if(WRITE_DISK) {
          paths <- compute_spaths2(origins, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, nfork, cl, FALSE,
            FALSE, unconnected_error, tr_directed, destinations, destination_nms_specified)
          if(output_lines) {
            terra::writeVector(paths, paste0(nm, ".", file_type))
          } else if(file_type_rds) {
            saveRDS(paths, paste0(nm, ".rds"))
          } else {
            data.table::fwrite(paths, paste0(nm, ".csv"))
          }
        } else {
          paths <- compute_spaths2(origins, rst, crd, TRUE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, nfork, cl,
            is.null(par_lvl), copy, unconnected_error, tr_directed, destinations, destination_nms_specified)
        }
      }
    } else {
      if(WRITE_DISK) {
        paths <- compute_spaths2(origins, rst, crd, FALSE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, nfork, cl, FALSE, FALSE,
          unconnected_error, tr_directed)
        if(output_lines) {
          terra::writeVector(paths, paste0(nm, ".", file_type))
        } else if(file_type_rds) {
          saveRDS(paths, paste0(nm, ".rds"))
        } else {
          data.table::fwrite(paths, paste0(nm, ".csv"))
        }
      } else {
        paths <- compute_spaths2(origins, rst, crd, FALSE, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, nfork, cl,
          is.null(par_lvl), copy, unconnected_error, tr_directed)
      }
    }
  }
  if(WRITE_DISK) paths <- NULL
  return(paths)
}

# Compute shortest paths
compute_spaths2 <- function(ORIGINS, rst, crd, dest_specified, origin_nms_specified, r_crs, output_lines, pairwise, NCORESG1, ncores, nfork, cl, nvect,
  copy, unconnected_error, tr_directed, DESTINATIONS = NULL, destination_nms_specified = TRUE) {
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
              if(min(pl, na.rm = TRUE) == 0L) return(NA)
              return(data.table::data.table(g = rep.int(O, pl), cls = do.call(c, s)))
            }
            if(nfork) {
              p <- list(data.table::rbindlist(report_na(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)), p), TRUE),
                use.names = FALSE)[, c("x", "y") := crd[cls,]][, c("g", "x", "y", "cls")], data.frame(origin_c = ORIGINS, destination_c = DESTINATIONS),
                data.frame(origin = on, destination = dn))
            } else {
              p <- list(data.table::rbindlist(report_na(parallel::mclapply(split(1:os_length, cut(1:os_length, ncores, labels = FALSE)), p,
                mc.silent = TRUE, mc.cores = ncores), TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, c("g", "x", "y", "cls")],
                data.frame(origin_c = ORIGINS, destination_c = DESTINATIONS), data.frame(origin = on, destination = dn))
            }
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
              if(min(pl, na.rm = TRUE) == 0L) return(NA)
              return(data.table::data.table(g = rep.int(O, pl), cls = do.call(c, s)))
            }
            if(nfork) {
              p <- terra::vect(as.matrix(data.table::rbindlist(report_na(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores,
                labels = FALSE)), p), TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line", atts = data.frame(origin = on,
                destination = dn), crs = r_crs)
            } else {
              p <- terra::vect(as.matrix(data.table::rbindlist(report_na(parallel::mclapply(split(1:os_length, cut(1:os_length, ncores, labels = FALSE)),
                p, mc.silent = TRUE, mc.cores = ncores), TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line",
                atts = data.frame(origin = on, destination = dn), crs = r_crs)
            }
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
                if(min(pl, na.rm = TRUE) == 0L) return(NA)
                return(data.table::data.table(g = o1 + rep.int(1:ds_length, pl), cls = do.call(c, s)))
              }), use.names = FALSE))
            }
            if(nfork) {
              p <- list(data.table::rbindlist(report_na(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)), p), TRUE),
                use.names = FALSE)[, c("x", "y") := crd[cls,]][, c("g", "x", "y", "cls")], data.frame(origin_c = rep.int(ORIGINS, rep.int(ds_length,
                os_length)), destination_c = rep.int(DESTINATIONS, os_length)), data.frame(origin = rep.int(on, rep.int(ds_length, os_length)),
                destination = rep.int(dn, os_length)))
            } else {
              p <- list(data.table::rbindlist(report_na(parallel::mclapply(split(1:os_length, cut(1:os_length, ncores, labels = FALSE)), p,
                mc.silent = TRUE, mc.cores = ncores), TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, c("g", "x", "y", "cls")],
                data.frame(origin_c = rep.int(ORIGINS, rep.int(ds_length, os_length)), destination_c = rep.int(DESTINATIONS, os_length)),
                data.frame(origin = rep.int(on, rep.int(ds_length, os_length)), destination = rep.int(dn, os_length)))
            }
          } else {
            p <- list(data.table::rbindlist(lapply(1:os_length, function(O) {
              s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
              o1 <- (O - 1L) * ds_length
              if(min(lengths(s), na.rm = TRUE) == 0L) report_points_unc(O, lengths(s))
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
                if(min(pl, na.rm = TRUE) == 0L) return(NA)
                return(data.table::data.table(g = o1 + rep.int(1:ds_length, pl), cls = do.call(c, s)))
              }), use.names = FALSE))
            }
            if(nfork) {
              p <- terra::vect(as.matrix(data.table::rbindlist(report_na(parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores,
                labels = FALSE)), p), TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line",
                atts = data.frame(origin = rep.int(on, rep.int(ds_length, os_length)), destination = rep.int(dn, os_length)), crs = r_crs)
            } else {
              p <- terra::vect(as.matrix(data.table::rbindlist(report_na(parallel::mclapply(split(1:os_length, cut(1:os_length, ncores, labels = FALSE)),
                p, mc.silent = TRUE, mc.cores = ncores), TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line",
                atts = data.frame(origin = rep.int(on, rep.int(ds_length, os_length)), destination = rep.int(dn, os_length)), crs = r_crs)
            }
          } else if(nvect) {
            p <- list(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
              s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
              o1 <- (O - 1L) * ds_length
              if(min(lengths(s), na.rm = TRUE) == 0L) report_points_unc(O, lengths(s))
              return(data.table::rbindlist(lapply(1:ds_length, function(D) crd[s[[D]],][, g := o1 + D]), use.names = FALSE))
            }), use.names = FALSE)[, c("g", "x", "y")]), data.frame(origin = rep.int(on, rep.int(ds_length, os_length)), destination = rep.int(dn,
              os_length)))
          } else {
            p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:os_length, function(O) {
              s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "vpath", algorithm = "dijkstra")$vpath
              o1 <- (O - 1L) * ds_length
              if(min(lengths(s), na.rm = TRUE) == 0L) report_points_unc(O, lengths(s))
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
                if(length(S$vpath[[1L]]) == 0L) return(NA)
                return(list(sum(igraph::edge_attr(rst, "weight", S$epath[[1L]])), S$vpath[[1L]]))
              })
              if(anyNA(s)) return(NA)
              return(list(vapply(s, `[[`, 1L, numeric(1L), USE.NAMES = FALSE), data.table::data.table(g = rep.int(O, sapply(s,
                function(o) length(o[[2L]]), USE.NAMES = FALSE)), cls = do.call(c, lapply(s, `[[`, 2L)))))
            }
            if(nfork) {
              p <- parallel::parLapply(cl, split(1:os_length, cut(1:os_length, ncores, labels = FALSE)), p)
            } else {
              p <- parallel::mclapply(split(1:os_length, cut(1:os_length, ncores, labels = FALSE)), p, mc.silent = TRUE, mc.cores = ncores)
            }
            if(anyNA(p)) stop("Not all origins are connected to their respective destinations")
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
              s <- diag(igraph::distances(rst, ORIGINS[O], DESTINATIONS[O], mode = "out", algorithm = "dijkstra"), names = FALSE)
              if(unconnected_error && is.infinite(max(s, na.rm = TRUE))) return(NA)
              return(s)
            }
            if(nfork) {
              p <- data.table::data.table(origin = on, destination = dn, distance = do.call(c, report_na(parallel::parLapply(cl, split(1:os_length,
                cut(1:os_length, ncores, labels = FALSE)), p), unconnected_error)))
            } else {
              p <- data.table::data.table(origin = on, destination = dn, distance = do.call(c, report_na(parallel::mclapply(split(1:os_length, cut(1:os_length,
                ncores, labels = FALSE)), p, mc.silent = TRUE, mc.cores = ncores), unconnected_error)))
            }
          } else {
            p <- data.table::data.table(origin = on, destination = dn, distance = diag(igraph::distances(rst, ORIGINS, DESTINATIONS, mode = "out",
              algorithm = "dijkstra"), names = FALSE))
            if(unconnected_error && is.infinite(max(p[["distance"]], na.rm = TRUE))) {
              report_points_unc(1:os_length, as.integer(is.finite(p[["distance"]])), TRUE)
            }
          }
        }
      } else {
        if(copy) {
          if(NCORESG1) {
            p <- function(O) {
              P <- lapply(O, function(o) {
                s <- igraph::shortest_paths(rst, o, DESTINATIONS, output = "both", algorithm = "dijkstra")
                if(min(lengths(s$vpath), na.rm = TRUE) == 0L) return(NA)
                return(list(vapply(s$epath, function(S) sum(igraph::edge_attr(rst, "weight", S)), numeric(1L), USE.NAMES = FALSE), s$vpath))
              })
              if(anyNA(P)) return(NA)
              return(list(do.call(c, lapply(P, `[[`, 1L)), do.call(c, lapply(P, `[[`, 2L))))
            }
            if(nfork) {
              p <- parallel::parLapply(cl, split(ORIGINS, cut(seq_along(ORIGINS), ncores, labels = FALSE)), p)
            } else {
              p <- parallel::mclapply(split(ORIGINS, cut(seq_along(ORIGINS), ncores, labels = FALSE)), p, mc.silent = TRUE, mc.cores = ncores)
            }
            if(anyNA(p)) stop("Not all origins are connected to their respective destinations")
            p <- list(data.table::data.table(origin = rep(on, each = ds_length), origin_c = rep(ORIGINS, each = ds_length), destination = rep.int(dn,
              os_length), destination_c = rep.int(DESTINATIONS, os_length), distance = do.call(c, lapply(p, `[[`, 1L))), do.call(c, lapply(p, `[[`, 2L)))
            p[[2L]] <- data.table::data.table(g = rep.int(1:(os_length * ds_length), lengths(p[[2L]])), cls = do.call(c, p[[2L]]))
          } else {
            p <- lapply(1:os_length, function(O) {
              s <- igraph::shortest_paths(rst, ORIGINS[O], DESTINATIONS, output = "both", algorithm = "dijkstra")
              if(min(lengths(s$vpath), na.rm = TRUE) == 0L) report_points_unc(O, lengths(s$vpath))
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
                if(unconnected_error && is.infinite(max(s, na.rm = TRUE))) return(NA)
                return(data.table::as.data.table(s, na.rm = FALSE))
              }
              if(nfork) {
                p <- data.table::rbindlist(report_na(parallel::parLapply(cl, split(ORIGINS, cut(seq_along(ORIGINS), ncores, labels = FALSE)), p),
                  unconnected_error), use.names = FALSE)
              } else {
                p <- data.table::rbindlist(report_na(parallel::mclapply(split(ORIGINS, cut(seq_along(ORIGINS), ncores, labels = FALSE)), p,
                  mc.silent = TRUE, mc.cores = ncores), unconnected_error), use.names = FALSE)
              }
            } else {
              p <- function(D) {
                s <- igraph::distances(rst, ORIGINS, D, mode = "out", algorithm = "dijkstra")
                if(unconnected_error && is.infinite(max(s, na.rm = TRUE))) return(NA)
                return(s)
              }
              if(nfork) {
                p <- data.table::as.data.table(do.call(cbind, report_na(parallel::parLapply(cl, split(DESTINATIONS, cut(seq_along(DESTINATIONS), ncores,
                  labels = FALSE)), p), unconnected_error)), na.rm = FALSE)
              } else {
                p <- data.table::as.data.table(do.call(cbind, report_na(parallel::mclapply(split(DESTINATIONS, cut(seq_along(DESTINATIONS), ncores,
                  labels = FALSE)), p, mc.silent = TRUE, mc.cores = ncores), unconnected_error)), na.rm = FALSE)
              }
            }
          } else {
            p <- data.table::as.data.table(igraph::distances(rst, ORIGINS, DESTINATIONS, mode = "out", algorithm = "dijkstra"), na.rm = FALSE)
            if(unconnected_error && is.infinite(max(p, na.rm = TRUE))) {
              o <- min(which(is.infinite(p)) %% os_length, na.rm = TRUE) + 1L
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
    o1 <- os_length - 1L
    # Lines output
    if(output_lines) {
      if(copy) {
        if(NCORESG1) {
          if(tr_directed) {
            p <- function(O) {
              P <- lapply(O, function(o) {
                s <- igraph::shortest_paths(rst, ORIGINS[o], ORIGINS[-o], output = "vpath", algorithm = "dijkstra")$vpath
                pl <- lengths(s)
                if(min(pl, na.rm = TRUE) == 0L) return(NA)
                return(data.table::data.table(g = (o - 1L) * o1 + rep.int(1:o1, pl), cls = do.call(c, s)))
              })
              if(anyNA(P)) return(NA)
              return(data.table::rbindlist(P, use.names = FALSE))
            }
            if(nfork) {
              p <- list(data.table::rbindlist(report_na(parallel::parLapply(cl, split(1:o1, cut(1:o1, ncores, labels = FALSE)), p), TRUE),
                use.names = FALSE)[, c("x", "y") := crd[cls,]][, c("g", "x", "y", "cls")], data.frame(origin_c = rep(ORIGINS, each = o1),
                destination_c = rep.int(ORIGINS, os_length)[-seq.int(1L, by = os_length + 1L, length.out = os_length)]), data.frame(origin = rep(on,
                each = o1), destination = rep.int(on, os_length)[-seq.int(1L, by = os_length + 1L, length.out = os_length)]))
            } else {
              p <- list(data.table::rbindlist(report_na(parallel::mclapply(split(1:o1, cut(1:o1, ncores, labels = FALSE)), p, mc.silent = TRUE,
                mc.cores = ncores), TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, c("g", "x", "y", "cls")], data.frame(origin_c = rep(ORIGINS,
                each = o1), destination_c = rep.int(ORIGINS, os_length)[-seq.int(1L, by = os_length + 1L, length.out = os_length)]),
                data.frame(origin = rep(on, each = o1), destination = rep.int(on, os_length)[-seq.int(1L, by = os_length + 1L, length.out = os_length)]))
            }
          } else {
            p <- function(O) {
              P <- lapply(O, function(o) {
                s <- igraph::shortest_paths(rst, ORIGINS[o], ORIGINS[(o + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
                pl <- lengths(s)
                if(min(pl, na.rm = TRUE) == 0L) return(NA)
                i <- os_length - o
                return(data.table::data.table(g = sum(o1:i) - i + rep.int(1:i, pl), cls = do.call(c, s)))
              })
              if(anyNA(P)) return(NA)
              return(data.table::rbindlist(P, use.names = FALSE))
            }
            if(nfork) {
              p <- list(data.table::rbindlist(report_na(parallel::parLapply(cl, split(1:o1, cut(1:o1, ncores, labels = FALSE)), p), TRUE),
                use.names = FALSE)[, c("x", "y") := crd[cls,]][, c("g", "x", "y", "cls")], stats::setNames(as.data.frame(t(utils::combn(ORIGINS, 2L))),
                c("origin_c", "destination_c")), stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")))
            } else {
              p <- list(data.table::rbindlist(report_na(parallel::mclapply(split(1:o1, cut(1:o1, ncores, labels = FALSE)), p, mc.silent = TRUE,
                mc.cores = ncores), TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, c("g", "x", "y", "cls")],
                stats::setNames(as.data.frame(t(utils::combn(ORIGINS, 2L))), c("origin_c", "destination_c")),
                stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")))
            }
          }
        } else {
          if(tr_directed) {
            p <- list(data.table::rbindlist(lapply(1:o1, function(O) {
              i <- (O - 1L) * o1
              s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[-O], output = "vpath", algorithm = "dijkstra")$vpath
              if(min(lengths(s), na.rm = TRUE) == 0L) report_points_unc(O, lengths(s), dest_specified = FALSE, d = (1:os_length)[-O])
              return(data.table::rbindlist(lapply(1:o1, function(D) crd[s[[D]],][, c("g", "cls") := list(i + D, s[[D]])]), use.names = FALSE))
            }), use.names = FALSE)[, c("g", "x", "y", "cls")], data.frame(origin_c = rep(ORIGINS, each = o1), destination_c = rep.int(ORIGINS,
              os_length)[-seq.int(1L, by = os_length + 1L, length.out = os_length)]), data.frame(origin = rep(on, each = o1), destination = rep.int(on,
              os_length)[-seq.int(1L, by = os_length + 1L, length.out = os_length)]))
          } else {
            p <- list(data.table::rbindlist(lapply(1:o1, function(O) {
              npO <- os_length - O
              i <- sum(o1:npO) - npO
              s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
              if(min(lengths(s), na.rm = TRUE) == 0L) report_points_unc(O, lengths(s), dest_specified = FALSE, d = (O + 1L):os_length)
              return(data.table::rbindlist(lapply(1:npO, function(D) crd[s[[D]],][, c("g", "cls") := list(i + D, s[[D]])]), use.names = FALSE))
            }), use.names = FALSE)[, c("g", "x", "y", "cls")], stats::setNames(as.data.frame(t(utils::combn(ORIGINS, 2L))), c("origin_c",
              "destination_c")), stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")))
          }
        }
      } else {
        if(NCORESG1) {
          if(tr_directed) {
            p <- function(O) {
              P <- lapply(O, function(o) {
                s <- igraph::shortest_paths(rst, ORIGINS[o], ORIGINS[-o], output = "vpath", algorithm = "dijkstra")$vpath
                pl <- lengths(s)
                if(min(pl, na.rm = TRUE) == 0L) return(NA)
                return(data.table::data.table(g = (o - 1L) * o1 + rep.int(1:o1, pl), cls = do.call(c, s)))
              })
              if(anyNA(P)) return(NA)
              return(data.table::rbindlist(P, use.names = FALSE))
            }
            if(nfork) {
              p <- terra::vect(as.matrix(data.table::rbindlist(report_na(parallel::parLapply(cl, split(1:o1, cut(1:o1, ncores, labels = FALSE)), p),
                TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line", atts = data.frame(origin = rep(on, each = o1),
                destination = rep.int(on, os_length)[-seq.int(1L, by = os_length + 1L, length.out = os_length)]), crs = r_crs)
            } else {
              p <- terra::vect(as.matrix(data.table::rbindlist(report_na(parallel::mclapply(split(1:o1, cut(1:o1, ncores, labels = FALSE)), p,
                mc.silent = TRUE, mc.cores = ncores), TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line",
                atts = data.frame(origin = rep(on, each = o1), destination = rep.int(on, os_length)[-seq.int(1L, by = os_length + 1L,
                length.out = os_length)]), crs = r_crs)
            }
          } else {
            p <- function(O) {
              P <- lapply(O, function(o) {
                s <- igraph::shortest_paths(rst, ORIGINS[o], ORIGINS[(o + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
                pl <- lengths(s)
                if(min(pl, na.rm = TRUE) == 0L) return(NA)
                i <- os_length - o
                return(data.table::data.table(g = sum(o1:i) - i + rep.int(1:i, pl), cls = do.call(c, s)))
              })
              if(anyNA(P)) return(NA)
              return(data.table::rbindlist(P, use.names = FALSE))
            }
            if(nfork) {
              p <- terra::vect(as.matrix(data.table::rbindlist(report_na(parallel::parLapply(cl, split(1:o1, cut(1:o1, ncores, labels = FALSE)), p),
                TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line",
                atts = stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")), crs = r_crs)
            } else {
              p <- terra::vect(as.matrix(data.table::rbindlist(report_na(parallel::mclapply(split(1:o1, cut(1:o1, ncores, labels = FALSE)), p,
                mc.silent = TRUE, mc.cores = ncores), TRUE), use.names = FALSE)[, c("x", "y") := crd[cls,]][, cls := NULL]), type = "line",
                atts = stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")), crs = r_crs)
            }
          }
        } else if(nvect) {
          if(tr_directed) {
            p <- list(as.matrix(data.table::rbindlist(lapply(1:o1, function(O) {
              i <- (O - 1L) * o1
              s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[-O], output = "vpath", algorithm = "dijkstra")$vpath
              if(min(lengths(s), na.rm = TRUE) == 0L) report_points_unc(O, lengths(s), dest_specified = FALSE, d = (1:os_length)[-O])
              return(data.table::rbindlist(lapply(1:o1, function(D) crd[s[[D]],][, g := i + D]), use.names = FALSE))
            }), use.names = FALSE)[, c("g", "x", "y")]), data.frame(origin = rep(on, each = o1), destination = rep.int(on, os_length)[-seq.int(1L,
              by = os_length + 1L, length.out = os_length)]))
          } else {
            p <- list(as.matrix(data.table::rbindlist(lapply(1:o1, function(O) {
              npO <- os_length - O
              i <- sum(o1:npO) - npO
              s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
              if(min(lengths(s), na.rm = TRUE) == 0L) report_points_unc(O, lengths(s), dest_specified = FALSE, d = (O + 1L):os_length)
              return(data.table::rbindlist(lapply(1:npO, function(D) crd[s[[D]],][, g := i + D]), use.names = FALSE))
            }), use.names = FALSE)[, c("g", "x", "y")]), stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin", "destination")))
          }
        } else {
          if(tr_directed) {
            p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:o1, function(O) {
              i <- (O - 1L) * o1
              s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[-O], output = "vpath", algorithm = "dijkstra")$vpath
              if(min(lengths(s), na.rm = TRUE) == 0L) report_points_unc(O, lengths(s), dest_specified = FALSE, d = (1:os_length)[-O])
              return(data.table::rbindlist(lapply(1:o1, function(D) crd[s[[D]],][, g := i + D]), use.names = FALSE))
            }), use.names = FALSE)[, c("g", "x", "y")]), type = "line", atts = data.frame(origin = rep(on, each = o1), destination = rep.int(on,
              os_length)[-seq.int(1L, by = os_length + 1L, length.out = os_length)]), crs = r_crs)
          } else {
            p <- terra::vect(as.matrix(data.table::rbindlist(lapply(1:o1, function(O) {
              npO <- os_length - O
              i <- sum(o1:npO) - npO
              s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[(O + 1L):os_length], output = "vpath", algorithm = "dijkstra")$vpath
              if(min(lengths(s), na.rm = TRUE) == 0L) report_points_unc(O, lengths(s), dest_specified = FALSE, d = (O + 1L):os_length)
              return(data.table::rbindlist(lapply(1:npO, function(D) crd[s[[D]],][, g := i + D]), use.names = FALSE))
            }), use.names = FALSE)[, c("g", "x", "y")]), type = "line", atts = stats::setNames(as.data.frame(t(utils::combn(on, 2L))), c("origin",
              "destination")), crs = r_crs)
          }
        }
      }
    # Distances output
    } else {
      if(copy) {
        if(NCORESG1) {
          if(tr_directed) {
            p <- function(O) {
              P <- lapply(O, function(o) {
                s <- igraph::shortest_paths(rst, ORIGINS[o], ORIGINS[-o], output = "both", algorithm = "dijkstra")
                pl <- lengths(s$vpath)
                if(min(pl, na.rm = TRUE) == 0L) return(NA)
                return(list(data.table::data.table(origin = on[o], origin_c = ORIGINS[o], destination = on[-o], destination_c = ORIGINS[-o],
                  distance = vapply(s$epath, function(S) sum(igraph::edge_attr(rst, "weight", S)), numeric(1L), USE.NAMES = FALSE)),
                  data.table::data.table(g = (o - 1L) * o1 + rep.int(1:o1, pl), cls = do.call(c, s$vpath))))
              })
              if(anyNA(P)) return(NA)
              return(list(data.table::rbindlist(lapply(P, `[[`, 1L), use.names = FALSE), data.table::rbindlist(lapply(P, `[[`, 2L), use.names = FALSE)))
            }
          } else {
            p <- function(O) {
              P <- lapply(O, function(o) {
                s <- igraph::shortest_paths(rst, ORIGINS[o], ORIGINS[(o + 1L):os_length], output = "both", algorithm = "dijkstra")
                pl <- lengths(s$vpath)
                if(min(pl, na.rm = TRUE) == 0L) return(NA)
                i <- os_length - o
                return(list(data.table::data.table(origin = on[o], origin_c = ORIGINS[o], destination = on[(o + 1L):os_length],
                  destination_c = ORIGINS[(o + 1L):os_length], distance = vapply(s$epath, function(S) sum(igraph::edge_attr(rst, "weight", S)),
                  numeric(1L), USE.NAMES = FALSE)), data.table::data.table(g = sum(o1:i) - i + rep.int(1:i, pl), cls = do.call(c, s$vpath))))
              })
              if(anyNA(P)) return(NA)
              return(list(data.table::rbindlist(lapply(P, `[[`, 1L), use.names = FALSE), data.table::rbindlist(lapply(P, `[[`, 2L), use.names = FALSE)))
            }
          }
          if(nfork) {
            p <- parallel::parLapply(cl, split(1:o1, cut(1:o1, ncores, labels = FALSE)), p)
          } else {
            p <- parallel::mclapply(split(1:o1, cut(1:o1, ncores, labels = FALSE)), p, mc.silent = TRUE, mc.cores = ncores)
          }
          if(anyNA(p)) stop("Not all origins are connected to their respective destinations")
        } else {
          if(tr_directed) {
            p <- lapply(1:o1, function(O) {
              s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[-O], output = "both", algorithm = "dijkstra")
              pl <- lengths(s$vpath)
              if(min(pl, na.rm = TRUE) == 0L) report_points_unc(O, pl, dest_specified = FALSE, d = (1:os_length)[-O])
              return(list(data.table::data.table(origin = on[O], origin_c = ORIGINS[O], destination = on[-O], destination_c = ORIGINS[-O],
                distance = vapply(s$epath, function(S) sum(igraph::edge_attr(rst, "weight", S)), numeric(1L), USE.NAMES = FALSE)),
                data.table::data.table(g = (O - 1L) * o1 + rep.int(1:o1, pl), cls = do.call(c, s$vpath))))
            })
          } else {
            p <- lapply(1:o1, function(O) {
              npO <- os_length - O
              i <- sum(o1:npO) - npO
              D <- (O + 1L):os_length
              s <- igraph::shortest_paths(rst, ORIGINS[O], ORIGINS[D], output = "both", algorithm = "dijkstra")
              pl <- lengths(s$vpath)
              if(min(pl, na.rm = TRUE) == 0L) report_points_unc(O, pl, dest_specified = FALSE, d = D)
              return(list(data.table::data.table(origin = on[O], origin_c = ORIGINS[O], destination = on[D], destination_c = ORIGINS[D],
                distance = vapply(s$epath, function(S) sum(igraph::edge_attr(rst, "weight", S)), numeric(1L), USE.NAMES = FALSE)),
                data.table::data.table(g = i + rep.int(1:npO, pl), cls = do.call(c, s$vpath))))
            })
          }
        }
        p <- list(data.table::rbindlist(lapply(p, `[[`, 1L), use.names = FALSE), data.table::rbindlist(lapply(p, `[[`, 2L), use.names = FALSE))
      } else {
        if(NCORESG1) {
          if(tr_directed) {
            p <- function(o) {
              s <- as.vector(igraph::distances(rst, ORIGINS[o], ORIGINS, mode = "out", algorithm = "dijkstra"))
              l <- length(o)
              s <- data.table::data.table(origin = rep.int(o, os_length), destination = rep(on, each = l), distance = s)[-seq.int((min(o,
                na.rm = TRUE) - 1L) * l + 1L, by = l + 1L, length.out = l),]
              data.table::setorder(s, origin)
              if(origin_nms_specified) s[, origin := on[origin]]
              if(unconnected_error && is.infinite(max(s[["distance"]], na.rm = TRUE))) return(NA)
              return(s)
            }
          } else {
            p <- function(o) {
              m1 <- min(o, na.rm = TRUE) + 1L
              s <- igraph::distances(rst, ORIGINS[o], ORIGINS[m1:os_length], mode = "out", algorithm = "dijkstra")
              ut <- upper.tri(s, TRUE)
              s <- data.table::data.table(origin = o[row(s)[ut]], destination = on[(m1:os_length)[col(s)[ut]]], distance = s[ut])
              rm(m1, ut)
              data.table::setorder(s, origin)
              if(origin_nms_specified) s[, origin := on[origin]]
              if(unconnected_error && is.infinite(max(s[["distance"]], na.rm = TRUE))) return(NA)
              return(s)
            }
          }
          if(nfork) {
            p <- data.table::rbindlist(report_na(parallel::parLapply(cl, split(1:o1, cut(1:o1, ncores, labels = FALSE)), p), unconnected_error),
              use.names = FALSE)
          } else {
            p <- data.table::rbindlist(report_na(parallel::mclapply(split(1:o1, cut(1:o1, ncores, labels = FALSE)), p, mc.silent = TRUE,
              mc.cores = ncores), unconnected_error), use.names = FALSE)
          }
        } else {
          if(tr_directed) {
            p <- as.vector(igraph::distances(rst, ORIGINS, ORIGINS, mode = "out", algorithm = "dijkstra"))
            p <- data.table::data.table(origin = rep.int(1:os_length, os_length), destination = rep(on, each = os_length), distance = p)[-seq.int(1L,
              by = os_length, length.out = os_length)]
          } else {
            p <- igraph::distances(rst, ORIGINS, ORIGINS[2:os_length], mode = "out", algorithm = "dijkstra")
            p <- p[upper.tri(p, TRUE)]
            p <- data.table::data.table(origin = list_origins(os_length), destination = on[rep.int(2:os_length, 1:o1)], distance = p)
          }
          data.table::setorder(p, origin)
          if(origin_nms_specified) p[, origin := on[origin]]
          if(unconnected_error && is.infinite(max(p[["distance"]], na.rm = TRUE))) {
            p <- p[p[which.max(p[["distance"]]), "origin"], nomatch = NULL, on = "origin"]
            report_points_unc(p[1L, "origin"][["origin"]], as.integer(is.finite(p[["distance"]])), dest_specified = FALSE, d = p[["destination"]])
          }
        }
      }
    }
  }
  return(p)
}

# Report points on cells masked by update_rst
report_points_ust <- function(v, u, o) {
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
  stop("update_rst masks the ", v)
}

# Report elements of update_rst masking points
report_update_m <- function(u) {
  p <- length(u)
  if(p > 1L) {
    if(p > 2L) {
      p <- paste0("s ", paste0(u[1:(p - 1L)], collapse = ", "), ", and ", u[p])
    } else {
      p <- paste0("s ", u[1L], " and ", u[2L])
    }
    p <- c(p, "")
  } else {
    p <- c(paste0(" ", u), "s")
  }
  stop("Element", p[1L], " of update_rst mask", p[2L], " cells of origin or destination points")
}

# Report update_rst elements not allowing to connect all points
report_update_rst <- function(u) {
  if(is.null(u)) {
    p <- "Not all points are connected when updating rst with update_rst"
  } else {
    p <- length(u)
    if(p > 1L) {
      if(p > 2L) {
        p <- paste0("s ", paste0(u[1:(p - 1L)], collapse = ", "), ", and ", u[p])
      } else {
        p <- paste0("s ", u[1L], " and ", u[2L])
      }
    } else {
      p <- paste0(" ", u)
    }
    p <- paste0("Not all points are connected when updating rst with element", p, " of update_rst")
  }
  stop(p)
}

# Report NA values
report_na <- function(p, unconnected_error) {
  if(unconnected_error && anyNA(p)) stop("Not all origins are connected to their respective destinations")
  return(p)
}

# Rbind or return NA
rbind_na <- function(p) {
  if(anyNA(p)) return(NA)
  return(data.table::rbindlist(p, use.names = FALSE))
}


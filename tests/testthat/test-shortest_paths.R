if(as.logical(Sys.getenv("NOT_CRAN", "false"))) {
  set.seed(2L)
  tr <- terra::rast(crs = "epsg:4326", resolution = 2, vals = sample(c(1L, NA_integer_), 16200L, TRUE, c(0.8, 0.2)))
  pts_o <- rnd_locations(5, output_type = "SpatVector")
  pts_o$name <- sample(letters, 5)
  pts_d <- rnd_locations(5, output_type = "SpatVector")
  pts_d$name <- sample(LETTERS, 5)
  urst <- terra::vect("POLYGON ((-179 1, 30 1, 30 0, -179 0, -179 1))", crs = "epsg:4326")
  urstl <- list(urst, terra::vect("POLYGON ((0 0, 0 89, 1 89, 1 0, 0 0))", crs = "epsg:4326"))
  tr_m <- matrix(terra::values(tr)[, 1L], nrow = terra::nrow(tr), byrow = TRUE)
  pts_o_m <- terra::crds(pts_o)
  pts_d_m <- terra::crds(pts_d)
  urst_v <- terra::extract(tr, urst, cells = TRUE, ID = FALSE, touches = TRUE)$cell
  urst_m <- 1:16200
  urst_m[urst_v] <- NA_integer_
  urst_m <- matrix(urst_m, nrow = 90L, byrow = TRUE)
  urstl_v <- lapply(urstl, function(x) terra::extract(tr, x, cells = TRUE, ID = FALSE, touches = TRUE)$cell)
  urstl_m <- lapply(urstl_v, function(x) {
    u_m <- 1:16200
    u_m[x] <- NA_integer_
    u_m <- matrix(u_m, nrow = 90L, byrow = TRUE)
    return(u_m)
  })
  
  para_comb_inner <- expand.grid(pre = c(TRUE, FALSE), early_stopping = c(TRUE, FALSE), ncores = c(1L, 2L), par_lvl = c("update_rst", "points"),
    path_type = c("int", "unsigned short int"), stringsAsFactors = FALSE)
  para_comb_outer <- expand.grid(dt = c("double", "float", "int", "unsigned short int"), ct = c("queen", "rook"), dc = c("spaths", "terra"), pw = c(TRUE,
    FALSE), stringsAsFactors = FALSE)
  
  # distances, SpatRaster
  mapply(function(dt, ct, dc, pw) {
    sapply(list(NULL, urst, urstl), function(u) {
      sapply(list(NULL, pts_d), function(y) {
        if(!(is.null(y) && pw)) {
          test_baseline <- new.env()
          test_that(paste0("spat,distances,", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",", (!is.null(u)) + is.list(u)), {
            assign("r", expect_no_error(shortest_paths(tr, pts_o, y, dist_comp = dc, pairwise = pw, contiguity = ct, update_rst = u, distance_type = dt)),
              envir = test_baseline)
          })
  
          mapply(function(pr, es, nc, pl, pt) {
            if(!(is.null(u) && pt == "points")) {
              test_that(paste0("spat,distances,", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",", pr, ",", es, ",", (!is.null(u)) + is.list(u), ",",
                nc, ",", pt), {
                expect_identical(shortest_paths(tr, pts_o, y, pairwise = pw, contiguity = ct, dist_comp = dc, pre = pr, early_stopping = es, update_rst = u,
                  ncores = nc, par_lvl = pl, path_type = pt, distance_type = dt), test_baseline$r)
              })
            }
          }, para_comb_inner$pre, para_comb_inner$early_stopping, para_comb_inner$ncores, para_comb_inner$par_lvl,
            para_comb_inner$path_type, USE.NAMES = FALSE)
          test_that(paste0("spat,distances,", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",", (!is.null(u)) + is.list(u), ",data.table"), {
            expect_s3_class(test_baseline$r, "data.table")
          })
          if(is.null(u)) ncls <- 3L else ncls <- 4L
          test_that(paste0("spat,distances,", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",", (!is.null(u)) + is.list(u), ",ncol"), {
            expect_identical(ncol(test_baseline$r), ncls)
          })
          if(is.null(y)) {
            nrws <- 10L
          } else if(pw) {
            nrws <- 5L
          } else {
            nrws <- 25L
          }
          if(!is.null(u)) {
            if(is.list(u)) {
              nrws <- nrws * 3L
            } else {
              nrws <- nrws * 2L
            }
          }
          test_that(paste0("spat,distances,", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",", (!is.null(u)) + is.list(u), ",nrow"), {
            expect_identical(nrow(test_baseline$r), nrws)
          })
          nms <- c("origins", "destinations", "distances")
          if(!is.null(u)) nms <- c(nms, "layer")
          test_that(paste0("spat,distances,", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",", (!is.null(u)) + is.list(u), ",names"), {
            expect_named(test_baseline$r, nms)
          })
          test_that(paste0("spat,distances,", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",", (!is.null(u)) + is.list(u), ",>0"), {
            expect_true(all(test_baseline$r[["distances"]] > 0))
          })
          test_that(paste0("spat,distances,", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",", (!is.null(u)) + is.list(u), ",!Inf"), {
            expect_true(all(is.finite(test_baseline$r[["distances"]])))
          })
        }
      }, USE.NAMES = FALSE)
    }, USE.NAMES = FALSE)
  }, para_comb_outer$dt, para_comb_outer$ct, para_comb_outer$dc, para_comb_outer$pw, USE.NAMES = FALSE)
  
  para_comb_outer <- expand.grid(dt = c("double", "float", "int", "unsigned short int"), ct = c("queen", "rook"), pw = c(TRUE, FALSE),
    stringsAsFactors = FALSE)
  
  # distances, matrix
  mapply(function(dt, ct, pw) {
    mapply(function(u, u_m) {
      mapply(function(y, y_m) {
        if(!(is.null(y) && pw)) {
          test_baseline <- new.env()
          test_that(paste0("matrix,distances,distance_type:", dt, ",contiguity:", ct, ",pairwise:", pw, ",without destinations:", is.null(y), ",u:",
            typeof(u_m)), {
            assign("r", expect_no_error(shortest_paths(tr_m, pts_o_m, y_m, extent = c(-180, 180, -90, 90), pairwise = pw, contiguity = ct, update_rst = u_m,
              distance_type = dt)), envir = test_baseline)
          })

          test_that(paste0("matrix,distances,distance_type:", dt, ",contiguity:", ct, ",pairwise:", pw, ",without destinations:", is.null(y), ",u:",
            typeof(u_m), ",rst_type"), {
            expect_identical(shortest_paths(tr, pts_o, y, pairwise = pw, contiguity = ct, update_rst = u, distance_type = dt), test_baseline$r)
          })

          mapply(function(pr, es, nc, pl, pt) {
            test_that(paste0("matrix,distances,distance_type:", dt, ",contiguity:", ct, ",pairwise:", pw, ",without destinations:", is.null(y), ",pre:", pr,
              ",early_stopping:", es, ",ncores:", nc, ",path_type:", pt, ",u:", typeof(u)), {
              expect_identical(shortest_paths(tr_m, pts_o_m, y_m, extent = c(-180, 180, -90, 90), pairwise = pw, contiguity = ct, pre = pr,
                early_stopping = es, update_rst = u_m, ncores = nc, par_lvl = pl, path_type = pt, distance_type = dt), test_baseline$r)
            })
          }, para_comb_inner$pre, para_comb_inner$early_stopping, para_comb_inner$ncores, para_comb_inner$par_lvl, para_comb_inner$path_type,
            USE.NAMES = FALSE)
          test_that(paste0("matrix,distances,distance_type:", dt, ",contiguity:", ct, ",pairwise:", pw, ",without destinations:", is.null(y),
            ",data.table"), {
            expect_s3_class(test_baseline$r, "data.table")
          })
          if(is.null(u)) ncls <- 3L else ncls <- 4L
          test_that(paste0("matrix,distances,distance_type:", dt, ",contiguity:", ct, ",pairwise:", pw, ",without destinations:", is.null(y), ",u:",
            typeof(u_m), ",ncol"), {
            expect_identical(ncol(test_baseline$r), ncls)
          })
          nms <- c("origins", "destinations", "distances")
          if(!is.null(u)) nms <- c(nms, "layer")
          test_that(paste0("matrix,distances,distance_type:", dt, ",contiguity:", ct, ",pairwise:", pw, ",without destinations:", is.null(y), ",u:",
            typeof(u_m), ",names"), {
            expect_named(test_baseline$r, nms)
          })
          test_that(paste0("matrix,distances,distance_type:", dt, ",contiguity:", ct, ",pairwise:", pw, ",without destinations:", is.null(y), ",u:",
            typeof(u_m), ",>0"), {
            expect_true(all(test_baseline$r[["distances"]] > 0))
          })
          test_that(paste0("matrix,distances,distance_type:", dt, ",contiguity:", ct, ",pairwise:", pw, ",without destinations:", is.null(y), ",u:",
            typeof(u_m), ",!Inf"), {
            expect_true(all(is.finite(test_baseline$r[["distances"]])))
          })
        }
      }, list(NULL, pts_d), list(NULL, pts_d_m), USE.NAMES = FALSE)
    }, list(NULL, urst, urst, urstl, urstl), list(NULL, urst_m, urst_v, urstl_m, urstl_v), USE.NAMES = FALSE)
  }, para_comb_outer$dt, para_comb_outer$ct, para_comb_outer$pw, USE.NAMES = FALSE)

  # lines, SpatRaster
  para_comb_outer <- expand.grid(o = c("lines", "both"), dt = c("double", "float", "int", "unsigned short int"), ct = c("queen", "rook"), dc = c("spaths",
    "terra"), pw = c(TRUE, FALSE), stringsAsFactors = FALSE)

  mapply(function(o, dt, ct, dc, pw) {
    sapply(list(NULL, pts_d), function(y) {
      if(!(is.null(y) && pw)) {
        test_baseline <- new.env()
        test_that(paste0("spat,output:", o, ",distance_type:", dt, ",contiguity:", ct, ",dist_comp:", dc, ",pairwise:", pw, ",without destinations:",
          is.null(y)), {
          r <- expect_no_error(shortest_paths(tr, pts_o, y, output = o, dist_comp = dc, pairwise = pw, contiguity = ct, distance_type = dt))
          assign("r", r, envir = test_baseline)
          assign("r_geom", expect_no_error(terra::geom(r)), envir = test_baseline)
          assign("r_values", expect_no_error(terra::values(r)), envir = test_baseline)
        })

        mapply(function(pr, es, nc, pt) {
          test_that(paste0("spat,output:", o, ",distance_type:", dt, ",contiguity:", ct, ",dist_comp:", dc, ",pairwise:", pw, ",without destinations:",
            is.null(y), ",", pr, ",", es, ",", nc, ",", pt), {
            r_x <- expect_no_error(shortest_paths(tr, pts_o, y, output = o, pairwise = pw, contiguity = ct, dist_comp = dc, tr_directed = td, pre = pr,
              early_stopping = es, ncores = nc, path_type = pt, distance_type = dt))
            expect_identical(terra::geom(r_x), test_baseline$r_geom)
            expect_identical(terra::values(r_x), test_baseline$r_values)
          })
        }, para_comb_inner$pre, para_comb_inner$early_stopping, para_comb_inner$ncores, para_comb_inner$path_type,
          USE.NAMES = FALSE)
        test_that(paste0("spat,", o, ",", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",SpatVector"), {
          expect_s4_class(test_baseline$r, "SpatVector")
        })
        test_that(paste0("spat,", o, ",", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",>0"), {
          expect_true(all(test_baseline$r_values[["distances"]] > 0))
        })
        test_that(paste0("spat,", o, ",", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",!Inf"), {
          expect_true(all(is.finite(test_baseline$r_values[["distances"]])))
        })
        if(o == "both") {
          test_that(paste0("spat,", o, ",", dt, ",", ct, ",", dc, ",", pw, ",", is.null(y), ",distances"), {
            expect_identical(data.table::as.data.table(test_baseline$r_values), shortest_paths(tr, pts_o, y, dist_comp = dc, pairwise = pw, contiguity = ct,
              distance_type = dt))
          })
        }
      }
    }, USE.NAMES = FALSE)
  }, para_comb_outer$o, para_comb_outer$dt, para_comb_outer$ct, para_comb_outer$dc, para_comb_outer$pw, USE.NAMES = FALSE)

  # lines, matrix
  para_comb_outer <- expand.grid(o = c("lines", "both"), dt = c("double", "float", "int", "unsigned short int"), ct = c("queen", "rook"), pw = c(TRUE,
    FALSE), stringsAsFactors = FALSE)

  mapply(function(o, dt, ct, pw) {
    mapply(function(y, y_m) {
      if(!(is.null(y) && pw)) {
        test_baseline <- new.env()
        test_that(paste0("matrix,output:", o, ",contiguity:", ct, ",pairwise:", pw, ",distance_type:", dt, ",without destinations:",
          is.null(y)), {
          assign("r", expect_no_error(shortest_paths(tr_m, pts_o_m, y_m, output = o, extent = c(-180, 180, -90, 90), pairwise = pw, contiguity = ct,
            distance_type = dt)), envir = test_baseline)
        })

        test_that(paste0("matrix,output:", o, ",contiguity:", ct, ",pairwise:", pw, ",distance_type:", dt, ",without destinations:",
          is.null(y), ",rst_type"), {
          expect_identical(shortest_paths(tr, pts_o, y, output = o, output_class = "list", pairwise = pw, contiguity = ct,
            distance_type = dt)[c("attributes", "lines")], test_baseline$r[c("attributes", "lines")])
        })

        mapply(function(pr, es, nc, pt) {
          test_that(paste0("matrix,output:", o, ",contiguity:", ct, ",pairwise:", pw, ",distance_type:", dt, ",without destinations:",
          is.null(y), ",pre:", pr, ",early_stopping:", es, ",ncores:", nc, ",path_type:", pt), {
            expect_identical(shortest_paths(tr_m, pts_o_m, y_m, output = o, extent = c(-180, 180, -90, 90), pairwise = pw, contiguity = ct, pre = pr,
              early_stopping = es, ncores = nc, path_type = pt, distance_type = dt), test_baseline$r)
          })
        }, para_comb_inner$pre, para_comb_inner$early_stopping, para_comb_inner$ncores, para_comb_inner$path_type, USE.NAMES = FALSE)
        test_that(paste0("matrix,output:", o, ",contiguity:", ct, ",pairwise:", pw, ",distance_type:", dt, ",without destinations:",
          is.null(y), ",list"), {
          expect_type(test_baseline$r, "list")
        })
        test_that(paste0("matrix,output:", o, ",contiguity:", ct, ",pairwise:", pw, ",distance_type:", dt, ",without destinations:",
          is.null(y), ",>0"), {
          expect_true(all(test_baseline$r[["attributes"]][["distances"]] > 0))
        })
        test_that(paste0("matrix,output:", o, ",contiguity:", ct, ",pairwise:", pw, ",distance_type:", dt, ",without destinations:",
          is.null(y), ",!Inf"), {
          expect_true(all(is.finite(test_baseline$r[["attributes"]][["distances"]])))
        })
        if(o == "both") {
          test_that(paste0("matrix,output:", o, ",contiguity:", ct, ",pairwise:", pw, ",distance_type:", dt, ",without destinations:",
          is.null(y), ",distances"), {
            expect_identical(data.table::as.data.table(test_baseline$r[["attributes"]]), shortest_paths(tr, pts_o, y, pairwise = pw, contiguity = ct,
              distance_type = dt))
          })
        }
      }
    }, list(NULL, pts_d), list(NULL, pts_d_m), USE.NAMES = FALSE)
  }, para_comb_outer$o, para_comb_outer$dt, para_comb_outer$ct, para_comb_outer$pw, USE.NAMES = FALSE)

  lapply(c("distances", "lines", "both"), function(o) {
    test_that(paste0("bidirectional,", o), {
      expect_no_error(shortest_paths(tr, pts_o, output = o, bidirectional = TRUE))
    })
  })

  test_that("names", {
    r <- expect_no_error(shortest_paths(tr, pts_o, origin_names = "name"))
    expect_in(r[["origins"]], letters)
    expect_in(r[["destinations"]], letters)
    r <- expect_no_error(shortest_paths(tr, pts_o, pts_d, origin_names = "name", destination_names = "name"))
    expect_in(r[["origins"]], letters)
    expect_in(r[["destinations"]], LETTERS)
    r <- expect_no_error(shortest_paths(tr, pts_o, pts_d, origin_names = "name", destination_names = "name", pairwise = TRUE))
    expect_in(r[["origins"]], letters)
    expect_in(r[["destinations"]], LETTERS)
    pts_oi <- pts_o
    pts_oi$name <- sample.int(10L, 5L)
    r <- expect_no_error(shortest_paths(tr, pts_oi, origin_names = "name"))
    expect_type(r[["origins"]], "integer")
    expect_type(r[["destinations"]], "integer")
    r <- expect_no_error(shortest_paths(tr, pts_oi, pts_d, origin_names = "name", destination_names = "name"))
    expect_type(r[["origins"]], "integer")
    expect_in(r[["destinations"]], LETTERS)
    r <- expect_no_error(shortest_paths(tr, pts_oi, pts_d, origin_names = "name", destination_names = "name", pairwise = TRUE))
    expect_type(r[["origins"]], "integer")
    expect_in(r[["destinations"]], LETTERS)
    pts_di <- pts_d
    pts_di$name <- sample.int(10L, 5L)
    r <- expect_no_error(shortest_paths(tr, pts_oi, pts_di, origin_names = "name", destination_names = "name"))
    expect_type(r[["origins"]], "integer")
    expect_type(r[["destinations"]], "integer")
    r <- expect_no_error(shortest_paths(tr, pts_oi, pts_di, origin_names = "name", destination_names = "name", pairwise = TRUE))
    expect_type(r[["origins"]], "integer")
    expect_type(r[["destinations"]], "integer")
    pts_oi$name <- as.numeric(pts_oi$name)
    pts_di$name <- as.numeric(pts_di$name)
    r <- expect_no_error(shortest_paths(tr, pts_oi, origin_names = "name"))
    expect_type(r[["origins"]], "double")
    expect_type(r[["destinations"]], "double")
    r <- expect_no_error(shortest_paths(tr, pts_oi, pts_di, origin_names = "name", destination_names = "name"))
    expect_type(r[["origins"]], "double")
    expect_type(r[["destinations"]], "double")
    r <- expect_no_error(shortest_paths(tr, pts_oi, pts_di, origin_names = "name", destination_names = "name", pairwise = TRUE))
    expect_type(r[["origins"]], "double")
    expect_type(r[["destinations"]], "double")
    expect_no_error(shortest_paths(tr, pts_o, output = "lines", origin_names = "name"))
    expect_no_error(shortest_paths(tr, pts_o, output = "both", origin_names = "name"))
    expect_no_error(shortest_paths(tr, pts_o, pts_d, output = "lines", origin_names = "name"))
    expect_no_error(shortest_paths(tr, pts_o, pts_d, output = "both", origin_names = "name"))
    expect_no_error(shortest_paths(tr, pts_o, pts_d, output = "lines", origin_names = "name", destination_names = "name"))
    expect_no_error(shortest_paths(tr, pts_o, pts_d, output = "both", origin_names = "name", destination_names = "name"))
  })

  test_that("errors", {
    expect_error(shortest_paths(tr, pts_o, pairwise = TRUE), "destinations must not be NULL, if pairwise is TRUE")
  })

  para_comb <- expand.grid(op = c("distances", "lines", "both"), es = c(TRUE, FALSE), stringsAsFactors = FALSE)

  mapply(function(op, es) {
    mapply(function(y, y_m) {
      test_baseline <- new.env()
      test_that(paste0("spat,transition_function,output:", op, ",early_stopping:", es, ",tr_directed:TRUE"), {
        assign("rts", expect_no_error(shortest_paths(tr, pts_o, y, output = op,
          tr_fun = function(d, x1, x2, y1, y2, v1, v2, nc) abs(d + x1 + x2 + y1 %% y2 + v1 + v2 + nc), early_stopping = es)), envir = test_baseline)
      })
      test_that(paste0("matrix,transition_function,output:", op, ",without destinations:", is.null(y), ",early_stopping:", es, ",tr_directed:TRUE"), {
        assign("rtm", expect_no_error(shortest_paths(tr_m, pts_o_m, y_m, output = op, extent = c(-180, 180, -90, 90),
          tr_fun = function(d, x1, x2, y1, y2, v1, v2, nc) abs(d + x1 + x2 + y1 %% y2 + v1 + v2 + nc), early_stopping = es)), envir = test_baseline)
      })
      test_that(paste0("spat,transition_function,output:", op, ",without destinations:", is.null(y), ",early_stopping:", es, ",tr_directed:FALSE"), {
        assign("rfs", expect_no_error(shortest_paths(tr, pts_o, y, output = op,
          tr_fun = function(d, x1, x2, y1, y2, v1, v2, nc) abs(d + x1 + x2 + y1 + y2 + v1 + v2 + nc), tr_directed = FALSE, early_stopping = es)),
          envir = test_baseline)
      })
      test_that(paste0("matrix,transition_function,output:", op, ",without destinations:", is.null(y), ",early_stopping:", es, ",tr_directed:FALSE"), {
        assign("rfm", expect_no_error(shortest_paths(tr_m, pts_o_m, y_m, output = op, extent = c(-180, 180, -90, 90),
          tr_fun = function(d, x1, x2, y1, y2, v1, v2, nc) abs(d + x1 + x2 + y1 + y2 + v1 + v2 + nc), tr_directed = FALSE, early_stopping = es)),
          envir = test_baseline)
      })
      if(op == "distances") {
        test_that(paste0("transition_function,output:", op, ",early_stopping:", es, ",tr_directed:TRUE"), {
          expect_identical(test_baseline$rts, test_baseline$rtm)
        })
        test_that(paste0("transition_function,output:", op, ",early_stopping:", es, ",tr_directed:FALSE"), {
          expect_identical(test_baseline$rfs, test_baseline$rfm)
        })
      } else {
        test_that(paste0("transition_function,output:", op, ",early_stopping:", es, ",tr_directed:TRUE"), {
          expect_identical(terra::values(test_baseline$rts), as.data.frame(test_baseline$rtm[["attributes"]]))
        })
        test_that(paste0("transition_function,output:", op, ",early_stopping:", es, ",tr_directed:FALSE"), {
          expect_identical(terra::values(test_baseline$rfs), as.data.frame(test_baseline$rfm[["attributes"]]))
        })
      }
    }, list(NULL, pts_d), list(NULL, pts_d_m), USE.NAMES = FALSE)
  }, para_comb$op, para_comb$es, USE.NAMES = FALSE)
}
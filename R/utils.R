# Functions called by both spaths_earth and spaths_general

# Update points' cell numbers
update_points <- function(v, v_list, crd, nms_specified) {
  if(v_list) {
    if(nms_specified) {
      v <- lapply(v, function(V) {
        return(data.table::setnames(crd[V, c("c_n_c", "nms"), on = "c_n==cls"], "c_n_c", "cls"))
      })
    } else {
      v <- lapply(v, function(V) {
        return(crd[.(V), "c_n_c", on = "c_n"][["c_n_c"]])
      })
    }
  } else {
    if(nms_specified) {
      v <- data.table::setnames(crd[v, c("c_n_c", "nms"), on = "c_n==cls"], "c_n_c", "cls")
    } else {
      v <- crd[.(v), "c_n_c", on = "c_n"][["c_n_c"]]
    }
  }
  return(v)
}

# Report unconnected origins and destinations
report_points_unc <- function(o, pl, pairwise = FALSE, dest_specified = TRUE, d = NULL, O = TRUE, both = TRUE, u = NULL) {
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
  if(!is.null(u)) {
    s <- paste0(s, " when updating rst with element ", u, " of update_rst")
  }
  stop(s)
}

# Avoid R CMD check note
utils::globalVariables(c(".", "O", "c_n", "c_n_c", "cls", "g", "origin", "destination", "distance", "from", "to"))

# Function called when loading the package (circumvents current igraph RAM bug)
.onLoad <- function(libname, pkgname) {
  igraph::igraph_options(return.vs.es = FALSE)
}

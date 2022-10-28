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
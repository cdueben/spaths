The package computes shortest paths between places taking barriers - such as landmasses for ships - or cost surfaces - e.g. based on terrain topography - into account. Apart from relating locations on Earth, spaths can also compute shortest paths more generally on spheres and planes.

`spaths` uses the `igraph` package's implementation of Dijkstra's (1959) algorithm in identifying the paths. It does not rely on or recycle code from any other spatial least cost path packages.

`spaths` emphasizes computational efficiency, compatibility with different generations of spatial classes in R, customizability, and user friendliness.

Install the package via `devtools::install_github("ChristianDueben/spaths", build_vignettes = T)` and consult the vignette with `vignette("spaths_introduction", "spaths")` for an introduction.

Let me know, if you come across any errors.

Extensions to this package are certainly welcome. You can either suggest modifications to `spaths` by sending me an email or submitting a pull request, or you can build a package that calls `spaths`.

**A major update is coming in late 2023. It drops the igraph dependency and moves a lot of R code into C++.**

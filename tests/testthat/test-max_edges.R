set.seed(2L)
input_grid_unprojected <- terra::rast(crs = "epsg:4326", resolution = 2, vals = sample(c(1L, NA_integer_), 16200L, TRUE, c(0.8, 0.2)))

test_that("max_edges_spatraster_unprojected", {
  expect_equal(max_edges(input_grid_unprojected), 103004)
})

test_that("max_edges_spatraster_projected", {
  expect_equal(max_edges(terra::project(input_grid_unprojected, "ESRI:54009")), 407084)
})

test_that("max_edges_matrix_unprojected", {
  expect_equal(max_edges(matrix(terra::values(input_grid_unprojected, FALSE), nrow = 90L, byrow = TRUE), extent = c(-180, 180, -90, 90)), 103004)
})

test_that("max_edges_list_unprojected", {
  grid <- matrix(terra::values(input_grid_unprojected, FALSE), nrow = 90L, byrow = TRUE)
  expect_equal(max_edges(list(grid, grid), extent = c(-180, 180, -90, 90)), 103004)
})

devtools::install(upgrade = "never", quick = TRUE)
Sys.setenv(NOT_CRAN = "true")
testthat::test_dir("./tests/testthat/", load_package = "installed", package = "scorecompdir")


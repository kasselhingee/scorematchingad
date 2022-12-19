devtools::install(upgrade = "never")
Sys.setenv(NOT_CRAN = "true")
testthat::test_dir("./tests/testthat/", load_package = "installed", package = "scorecompdir")


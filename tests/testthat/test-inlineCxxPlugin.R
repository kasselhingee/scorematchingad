test_that("getting inline plugin passes", {
  inline::registerPlugin("scorematchingad", inlineCxxPlugin)
  expect_no_error(inline::getPlugin("scorematchingad"))
})

test_that("simple cxx compiles and runs", {
  inline::registerPlugin("scorematchingad", inlineCxxPlugin)
  expect_no_error({fx <- inline::cxxfunction(signature(x = "integer", y = "numeric"), 
                             "return wrap( as<int>(x) * as<double>(y));",
                             plugin = "scorematchingad" )})
  expect_no_error(fx(2L, 5))
  #write(strsplit(fx@code, "\n")[[1]], file = "tmp.cpp")
})


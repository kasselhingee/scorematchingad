#' @noRd 
inlineCxxPlugin <- function(...) {
    plugin <-
        Rcpp::Rcpp.plugin.maker(
                  include.before = c("#include <cppad/cppad.hpp>", "#include <RcppEigen.h>", '#include <scorematchingad.h>'),
                  package        = "scorematchingad",
                  Makevars = "PKG_CXXFLAGS = -I../inst/include",
                  LinkingTo = c("RcppEigen", "scorematchingad", "Rcpp")
              )
    settings <- plugin()
    settings
}


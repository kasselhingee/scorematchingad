#' @noRd 
inlineCxxPlugin <- function(...) {
    plugin <-
        Rcpp::Rcpp.plugin.maker(
                  include.before = c("#include <cppad/cppad.hpp>", "#include <RcppEigen.h>"),
                  package        = "scorematchingad",
                  LinkingTo = c("RcppEigen", "scorematchingad", "Rcpp")
              )
    settings <- plugin()
    settings
}


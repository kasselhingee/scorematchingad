#' @noRd 
inlineCxxPlugin <- function(...) {
    plugin <-
        Rcpp::Rcpp.plugin.maker(
                  include.before = "#include <RcppEigen.h>",
                  package        = "scorematchingad",
                  LinkingTo = c("RcppEigen", "scorematchingad", "Rcpp")
              )
    settings <- plugin()
    settings
}


#' @noRd 
inlineCxxPlugin <- function(...) {
    plugin <-
        Rcpp::Rcpp.plugin.maker(
                  include.before = c('#include <scorematchingad.h>'),
                  include.after = paste("#include", c("<utils/PrintFor.hpp>", "<tapell.cpp>", "<utils/wrapas.hpp>", "<manifoldtransforms/transforms.hpp>", "<likelihoods/likelihoods.hpp>")),
                  package        = "scorematchingad",
                  LinkingTo = c("RcppEigen", "scorematchingad", "Rcpp")
              )
    settings <- plugin()
    settings
}


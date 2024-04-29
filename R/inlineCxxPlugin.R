#' @noRd 
inlineCxxPlugin <- function(...) {
    plugin <-
        Rcpp::Rcpp.plugin.maker(
                  include.before = c('#include <scorematchingad.h>'),
                  include.after = paste("#include", c("<utils/PrintFor.hpp>", "<tapellcpp.h>", "<utils/wrapas.hpp>", "<manifoldtransforms/transforms.hpp>", "<likelihoods/likelihoods.hpp>"), collapse = "\n"),
                  package        = "scorematchingad",
                  LinkingTo = c("RcppEigen", "scorematchingad", "Rcpp")
              )
    settings <- plugin()
    settings
}


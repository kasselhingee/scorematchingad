# R6 ADFun class. Making it here because templating of the Cpp class means it is hard to expose via Rcpp
# Using R6 because it has modify in place semantics which match the ADFun objects
# following https://adv-r.hadley.nz/r6.html
# @examples

ADFun <- R6::R6Class("ADFun",
  private = list( #private == not modifiable
    .ptr = NULL,
    .name = NULL,
    .xtape = vector("numeric"),
    .dyntape = vector("numeric")
  ),
  public = list(
    initialize = function(ptr, name = NULL, xtape = vector("numeric"), dyntape = vector("numeric")){
      stopifnot(isa(ptr, "externalptr"))
      stopifnot(is.null(name) | isa(name, "character"))
      stopifnot(isa(xtape, "numeric"))
      stopifnot(isa(dyntape, "numeric"))
      private$.ptr <- ptr
      private$.name <- name
      private$.xtape <- xtape
      private$.dyntape <- dyntape
    }
  ),
  active = list(
    ptr = function(value){if (missing(value)){private$.ptr} else {stop("`$ptr' is read only", call. = FALSE)}},
    name = function(value){if (missing(value)){private$.name} else {stop("`$name' is read only", call. = FALSE)}},
    xtape = function(value){if (missing(value)){private$.xtape} else {stop("`$xtape' is read only", call. = FALSE)}},
    dyntape = function(value){if (missing(value)){private$.dyntape} else {stop("`$dyntape' is read only", call. = FALSE)}}
  ),
  cloneable = FALSE
)



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{tape_gradoffset}
\alias{tape_gradoffset}
\title{Tape the Gradient Offset of a Quadratic CppAD Tape}
\usage{
tape_gradoffset(pfun)
}
\arguments{
\item{pfun}{An \code{Rcpp_ADFun} object.}
}
\value{
An \code{Rcpp_ADFun} object. The independent argument to the function are the dynamic parameters of \code{pfun}.
}
\description{
Tape the Gradient Offset of a Quadratic CppAD Tape
}
\details{
A quadratic function can be written as
\deqn{f(x;\theta) = \frac{1}{2} x^T W(\theta) x + b(\theta)^Tx + c.}
The gradient of \eqn{f(x; \theta)} with respect to \eqn{x} is
\deqn{\Delta f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T)x + b(\theta).}
The Hessian is
\deqn{H f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T),}
which does not depend on \eqn{x},
so the gradient of the function can be rewritten as
\deqn{\Delta f(x;\theta) = H f(x; \theta) x + b(\theta)^T.}
The tape calculates \eqn{b(\theta)} as
\deqn{b(\theta) = \Delta f(x;\theta) - H f(x; \theta) x,}
which does not depend on \eqn{x}.

For creating this tape, the values of \code{pfun$xtape} and \code{pfun$dyntape} are used.
}
\seealso{
Other tape builders: 
\code{\link{tape_Hessian}()},
\code{\link{tape_Jacobian}()},
\code{\link{tape_logJacdet}()},
\code{\link{tape_smd}()},
\code{\link{tape_swap}()},
\code{\link{tape_uld}()}
}
\concept{tape builders}

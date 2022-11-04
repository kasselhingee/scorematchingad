#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#' @useDynLib scorecompdir
#' @details
#' This package includes score matching estimators for particular distributions and a general capacity to implement additional score matching estimators.
#' Score matching is a popular estimation technique when normalising constants for the proposed model are difficult to calculate or compute.
#' Score matching was first developed by \insertCite{hyvarinen2005es}{scorecompdir} and was further developed for subsets of Euclidean space \insertCite{@hyvarinen2007ex; @yu2019ge; @yu2020ge, @liu2021es}{scorecompdir}, Riemannian manifolds \insertCite{@mardia2016sc; @mardia2018ne}{scorecompdir},
#' and Riemannian manifolds with boundary \insertCite{scealy2022sc}{scorecompdir}.
#' In the most general form (Riemannian manifolds with boundary) score matching minimises the weighted Hyvarinen divergence \insertCite{@Equation 7 @scealy2022sc}{scorecompdir} 
#' \deqn{\phi(f, f_0) =  \frac{1}{2} \int_M f_0(z)h(z)^2 \left\lVert P(z)\Big(\nabla_z \log(f) - \nabla_z \log(f_0)\Big)\right\rVert^2 dM(z),}
#' where
#' + \eqn{M} is the manifold, isometrically embedded in Euclidean space, and \eqn{dM(z)} is the unnormalised uniform measure on \eqn{M}.
#' + \eqn{P(z)} is the matrix that projects points onto the tangent space of the manifold at \eqn{z}.
#' + \eqn{f_0} is the *true* density of the data-generating process, defined with respect to \eqn{dM(z)}.
#' + \eqn{f} is the density of a posited model, again defined with respect to \eqn{dM(z)}.
#' + \eqn{h(z)} is a function, termed the *divergence weight function*, that is zero on the boundary of \eqn{M} \insertCite{@Section 3.2, @scealy2022sc}{scorecompdir}.
#' + \eqn{\nabla_z} is the Euclidean gradient operator.
#' + \eqn{\lVert \cdot \rVert} is the Euclidean norm.
#' 
#' Note that, because \eqn{P(z)} is the projection matrix, \eqn{\left\lVert P(z)\Big(\nabla_z \log(f) - \nabla_z \log(f_0)\Big)\right\rVert^2} is the natural inner product of the gradient of the log ratio of \eqn{f} and \eqn{f_0}.
#'
#' When the density functions \eqn{f} and \eqn{f_0} are smooth and positive inside \eqn{M},
#' and the divergence weight function is smooth or of particular forms for specific manifolds \insertCite{@Section 3.2, @scealy2022sc}{scorecompdir},
#' then minimising the weighted Hyvarinen divergence \eqn{\phi(f, f_0)} is equivalent to minimising \insertCite{@Theorem 1 @scealy2022sc}{scorecompdir}
#'  \deqn{\psi(f, f_0) = \frac{1}{2} \int f_0(z)\big(A(z_i) + B(z_i) + C(z_i)\big)dM(z),}
#' where 
#' \deqn{A(z_i) = \left(\nabla_z \log(f)\right)^T P(z) \left(\nabla_z \log(f)\right),}
#' \deqn{B(z_i) = h(z)^2\Delta_z\log(f),}
#' \deqn{C(z_i) = \left(\nabla_z h(z)^2)\right)^T P(z) \left(\nabla_z \log(f)\right).}
#' When \eqn{n} iid observations from \eqn{f_0} are available, the integration in \eqn{\psi(f, f_0)} can be approximated by an average over the observations, 
#'  \deqn{\psi(f, f_0) \approx \hat\psi(f, f_0) = \frac{1}{2n} \sum_{i = 1}^n A(z_i) + B(z_i) + C(z_i).}
#' The *score matching estimator* minimises \eqn{\hat\psi(f, f_0)}.
#' 
#' This package's major contributions are
#'  * Implementations of score matching estimators that use algorithmic differentiation to avoid tedious by-hand algebraic calculations.
#' The package uses `CppAD` and `Eigen` to differentiate model densities and an empiricial score matching objective function.
#' The score matching objective is minimised by a modern implementation of conjugate gradient descent available through [`Rcgmin::Rcgmin()`].
#'  * Function [`ppi()`], containing score matching estimators for the Polynomially-Tilted Pairwise Interaction (PPI) Model \insertCite{scealy2022sc}{scorecompdir}.
#'  * Score matching and hybrid score matching estimators for some directional distributions \insertCite{mardia2016sc}{scorecompdir}.
#'  * Function [`Windham()`], that implements a modification of Windham's robustifying method \insertCite{windham1995ro}{scorecompdir} for many exponential family distributions. 
#' For some models the density approaches infinity at some locations, creating difficulties for the weights in Windham's original method \insertCite{@scealy2022ro}{scorecompdir}.
#' \insertNoCite{*}{scorecompdir}
#' @references
#' \insertAllCited{}
"_PACKAGE"


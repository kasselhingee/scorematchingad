// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "cdabyppi_types.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// corels
bool corels();
RcppExport SEXP _cdabyppi_corels() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(corels());
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _cdabyppi_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// pmanifold
XPtr< manifold<a1type> > pmanifold(std::string manifoldname);
RcppExport SEXP _cdabyppi_pmanifold(SEXP manifoldnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type manifoldname(manifoldnameSEXP);
    rcpp_result_gen = Rcpp::wrap(pmanifold(manifoldname));
    return rcpp_result_gen;
END_RCPP
}
// ptapesmo
XPtr< CppAD::ADFun<double> > ptapesmo(svecd u, svecd theta, size_t n, XPtr< CppAD::ADFun<double> > pll, XPtr< manifold<a1type> > pman, std::string weightname, const double acut);
RcppExport SEXP _cdabyppi_ptapesmo(SEXP uSEXP, SEXP thetaSEXP, SEXP nSEXP, SEXP pllSEXP, SEXP pmanSEXP, SEXP weightnameSEXP, SEXP acutSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< svecd >::type u(uSEXP);
    Rcpp::traits::input_parameter< svecd >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< size_t >::type n(nSEXP);
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pll(pllSEXP);
    Rcpp::traits::input_parameter< XPtr< manifold<a1type> > >::type pman(pmanSEXP);
    Rcpp::traits::input_parameter< std::string >::type weightname(weightnameSEXP);
    Rcpp::traits::input_parameter< const double >::type acut(acutSEXP);
    rcpp_result_gen = Rcpp::wrap(ptapesmo(u, theta, n, pll, pman, weightname, acut));
    return rcpp_result_gen;
END_RCPP
}
// psmo
double psmo(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain);
RcppExport SEXP _cdabyppi_psmo(SEXP pfunSEXP, SEXP uSEXP, SEXP betainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type u(uSEXP);
    Rcpp::traits::input_parameter< svecd >::type betain(betainSEXP);
    rcpp_result_gen = Rcpp::wrap(psmo(pfun, u, betain));
    return rcpp_result_gen;
END_RCPP
}
// psmograd
svecd psmograd(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain);
RcppExport SEXP _cdabyppi_psmograd(SEXP pfunSEXP, SEXP uSEXP, SEXP betainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type u(uSEXP);
    Rcpp::traits::input_parameter< svecd >::type betain(betainSEXP);
    rcpp_result_gen = Rcpp::wrap(psmograd(pfun, u, betain));
    return rcpp_result_gen;
END_RCPP
}
// ptapell
XPtr< CppAD::ADFun<double> > ptapell(svecd z, svecd theta, std::string llname, XPtr< manifold<a1type> > pman);
RcppExport SEXP _cdabyppi_ptapell(SEXP zSEXP, SEXP thetaSEXP, SEXP llnameSEXP, SEXP pmanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< svecd >::type z(zSEXP);
    Rcpp::traits::input_parameter< svecd >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< std::string >::type llname(llnameSEXP);
    Rcpp::traits::input_parameter< XPtr< manifold<a1type> > >::type pman(pmanSEXP);
    rcpp_result_gen = Rcpp::wrap(ptapell(z, theta, llname, pman));
    return rcpp_result_gen;
END_RCPP
}
// swapDynamic
XPtr< CppAD::ADFun<double> > swapDynamic(XPtr< CppAD::ADFun<double> > pfun, svecd newvalue, svecd newdynparam);
RcppExport SEXP _cdabyppi_swapDynamic(SEXP pfunSEXP, SEXP newvalueSEXP, SEXP newdynparamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type newvalue(newvalueSEXP);
    Rcpp::traits::input_parameter< svecd >::type newdynparam(newdynparamSEXP);
    rcpp_result_gen = Rcpp::wrap(swapDynamic(pfun, newvalue, newdynparam));
    return rcpp_result_gen;
END_RCPP
}
// pJacobian
svecd pJacobian(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain);
RcppExport SEXP _cdabyppi_pJacobian(SEXP pfunSEXP, SEXP uSEXP, SEXP betainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type u(uSEXP);
    Rcpp::traits::input_parameter< svecd >::type betain(betainSEXP);
    rcpp_result_gen = Rcpp::wrap(pJacobian(pfun, u, betain));
    return rcpp_result_gen;
END_RCPP
}
// pForward0
double pForward0(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain);
RcppExport SEXP _cdabyppi_pForward0(SEXP pfunSEXP, SEXP uSEXP, SEXP betainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type u(uSEXP);
    Rcpp::traits::input_parameter< svecd >::type betain(betainSEXP);
    rcpp_result_gen = Rcpp::wrap(pForward0(pfun, u, betain));
    return rcpp_result_gen;
END_RCPP
}
// pHessian
svecd pHessian(XPtr< CppAD::ADFun<double> > pfun, svecd u, svecd betain);
RcppExport SEXP _cdabyppi_pHessian(SEXP pfunSEXP, SEXP uSEXP, SEXP betainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type u(uSEXP);
    Rcpp::traits::input_parameter< svecd >::type betain(betainSEXP);
    rcpp_result_gen = Rcpp::wrap(pHessian(pfun, u, betain));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cdabyppi_corels", (DL_FUNC) &_cdabyppi_corels, 0},
    {"_cdabyppi_rcpp_hello_world", (DL_FUNC) &_cdabyppi_rcpp_hello_world, 0},
    {"_cdabyppi_pmanifold", (DL_FUNC) &_cdabyppi_pmanifold, 1},
    {"_cdabyppi_ptapesmo", (DL_FUNC) &_cdabyppi_ptapesmo, 7},
    {"_cdabyppi_psmo", (DL_FUNC) &_cdabyppi_psmo, 3},
    {"_cdabyppi_psmograd", (DL_FUNC) &_cdabyppi_psmograd, 3},
    {"_cdabyppi_ptapell", (DL_FUNC) &_cdabyppi_ptapell, 4},
    {"_cdabyppi_swapDynamic", (DL_FUNC) &_cdabyppi_swapDynamic, 3},
    {"_cdabyppi_pJacobian", (DL_FUNC) &_cdabyppi_pJacobian, 3},
    {"_cdabyppi_pForward0", (DL_FUNC) &_cdabyppi_pForward0, 3},
    {"_cdabyppi_pHessian", (DL_FUNC) &_cdabyppi_pHessian, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_cdabyppi(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

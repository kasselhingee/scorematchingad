// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "scorecompdir_types.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pmanifold
XPtr< manifold<a1type> > pmanifold(std::string manifoldname);
RcppExport SEXP _scorecompdir_pmanifold(SEXP manifoldnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type manifoldname(manifoldnameSEXP);
    rcpp_result_gen = Rcpp::wrap(pmanifold(manifoldname));
    return rcpp_result_gen;
END_RCPP
}
// testmanifold
int testmanifold(XPtr< manifold<a1type> > pman, svecd u);
RcppExport SEXP _scorecompdir_testmanifold(SEXP pmanSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< manifold<a1type> > >::type pman(pmanSEXP);
    Rcpp::traits::input_parameter< svecd >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(testmanifold(pman, u));
    return rcpp_result_gen;
END_RCPP
}
// ptoM
svecd ptoM(XPtr< manifold<a1type> > pman, svecd u);
RcppExport SEXP _scorecompdir_ptoM(SEXP pmanSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< manifold<a1type> > >::type pman(pmanSEXP);
    Rcpp::traits::input_parameter< svecd >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(ptoM(pman, u));
    return rcpp_result_gen;
END_RCPP
}
// ptapesmo
XPtr< CppAD::ADFun<double> > ptapesmo(svecd u, svecd theta, XPtr< CppAD::ADFun<double> > pll, XPtr< manifold<a1type> > pman, std::string weightname, const double acut, bool verbose);
RcppExport SEXP _scorecompdir_ptapesmo(SEXP uSEXP, SEXP thetaSEXP, SEXP pllSEXP, SEXP pmanSEXP, SEXP weightnameSEXP, SEXP acutSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< svecd >::type u(uSEXP);
    Rcpp::traits::input_parameter< svecd >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pll(pllSEXP);
    Rcpp::traits::input_parameter< XPtr< manifold<a1type> > >::type pman(pmanSEXP);
    Rcpp::traits::input_parameter< std::string >::type weightname(weightnameSEXP);
    Rcpp::traits::input_parameter< const double >::type acut(acutSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(ptapesmo(u, theta, pll, pman, weightname, acut, verbose));
    return rcpp_result_gen;
END_RCPP
}
// ptapell
XPtr< CppAD::ADFun<double> > ptapell(svecd z, svecd theta, std::string llname, XPtr< manifold<a1type> > pman, std::vector<int> fixedtheta, bool verbose);
RcppExport SEXP _scorecompdir_ptapell(SEXP zSEXP, SEXP thetaSEXP, SEXP llnameSEXP, SEXP pmanSEXP, SEXP fixedthetaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< svecd >::type z(zSEXP);
    Rcpp::traits::input_parameter< svecd >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< std::string >::type llname(llnameSEXP);
    Rcpp::traits::input_parameter< XPtr< manifold<a1type> > >::type pman(pmanSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type fixedtheta(fixedthetaSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(ptapell(z, theta, llname, pman, fixedtheta, verbose));
    return rcpp_result_gen;
END_RCPP
}
// swapDynamic
XPtr< CppAD::ADFun<double> > swapDynamic(XPtr< CppAD::ADFun<double> > pfun, svecd newvalue, svecd newdynparam);
RcppExport SEXP _scorecompdir_swapDynamic(SEXP pfunSEXP, SEXP newvalueSEXP, SEXP newdynparamSEXP) {
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
vecd pJacobian(XPtr< CppAD::ADFun<double> > pfun, svecd value, svecd theta);
RcppExport SEXP _scorecompdir_pJacobian(SEXP pfunSEXP, SEXP valueSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type value(valueSEXP);
    Rcpp::traits::input_parameter< svecd >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(pJacobian(pfun, value, theta));
    return rcpp_result_gen;
END_RCPP
}
// pForward0
svecd pForward0(XPtr< CppAD::ADFun<double> > pfun, svecd value, svecd theta);
RcppExport SEXP _scorecompdir_pForward0(SEXP pfunSEXP, SEXP valueSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type value(valueSEXP);
    Rcpp::traits::input_parameter< svecd >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(pForward0(pfun, value, theta));
    return rcpp_result_gen;
END_RCPP
}
// pHessian
svecd pHessian(XPtr< CppAD::ADFun<double> > pfun, svecd value, svecd theta);
RcppExport SEXP _scorecompdir_pHessian(SEXP pfunSEXP, SEXP valueSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type value(valueSEXP);
    Rcpp::traits::input_parameter< svecd >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(pHessian(pfun, value, theta));
    return rcpp_result_gen;
END_RCPP
}
// pTaylorApprox
svecd pTaylorApprox(XPtr< CppAD::ADFun<double> > pfun, svecd value, svecd centre, svecd theta, size_t order);
RcppExport SEXP _scorecompdir_pTaylorApprox(SEXP pfunSEXP, SEXP valueSEXP, SEXP centreSEXP, SEXP thetaSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type value(valueSEXP);
    Rcpp::traits::input_parameter< svecd >::type centre(centreSEXP);
    Rcpp::traits::input_parameter< svecd >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< size_t >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(pTaylorApprox(pfun, value, centre, theta, order));
    return rcpp_result_gen;
END_RCPP
}
// pTapeJacobianSwap
XPtr< CppAD::ADFun<double> > pTapeJacobianSwap(XPtr< CppAD::ADFun<double> > pfun, svecd value, svecd theta);
RcppExport SEXP _scorecompdir_pTapeJacobianSwap(SEXP pfunSEXP, SEXP valueSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type value(valueSEXP);
    Rcpp::traits::input_parameter< svecd >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(pTapeJacobianSwap(pfun, value, theta));
    return rcpp_result_gen;
END_RCPP
}
// pTapeHessianSwap
XPtr< CppAD::ADFun<double> > pTapeHessianSwap(XPtr< CppAD::ADFun<double> > pfun, svecd value, svecd theta);
RcppExport SEXP _scorecompdir_pTapeHessianSwap(SEXP pfunSEXP, SEXP valueSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr< CppAD::ADFun<double> > >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< svecd >::type value(valueSEXP);
    Rcpp::traits::input_parameter< svecd >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(pTapeHessianSwap(pfun, value, theta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scorecompdir_pmanifold", (DL_FUNC) &_scorecompdir_pmanifold, 1},
    {"_scorecompdir_testmanifold", (DL_FUNC) &_scorecompdir_testmanifold, 2},
    {"_scorecompdir_ptoM", (DL_FUNC) &_scorecompdir_ptoM, 2},
    {"_scorecompdir_ptapesmo", (DL_FUNC) &_scorecompdir_ptapesmo, 7},
    {"_scorecompdir_ptapell", (DL_FUNC) &_scorecompdir_ptapell, 6},
    {"_scorecompdir_swapDynamic", (DL_FUNC) &_scorecompdir_swapDynamic, 3},
    {"_scorecompdir_pJacobian", (DL_FUNC) &_scorecompdir_pJacobian, 3},
    {"_scorecompdir_pForward0", (DL_FUNC) &_scorecompdir_pForward0, 3},
    {"_scorecompdir_pHessian", (DL_FUNC) &_scorecompdir_pHessian, 3},
    {"_scorecompdir_pTaylorApprox", (DL_FUNC) &_scorecompdir_pTaylorApprox, 5},
    {"_scorecompdir_pTapeJacobianSwap", (DL_FUNC) &_scorecompdir_pTapeJacobianSwap, 3},
    {"_scorecompdir_pTapeHessianSwap", (DL_FUNC) &_scorecompdir_pTapeHessianSwap, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_scorecompdir(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

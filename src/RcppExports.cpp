// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/scorematchingad.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// abort_recording
void abort_recording();
RcppExport SEXP _scorematchingad_abort_recording() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    abort_recording();
    return R_NilValue;
END_RCPP
}
// taylorApprox_currentdynparam
vecd taylorApprox_currentdynparam(pADFun& pfun, vecd x, vecd centre, const size_t order);
RcppExport SEXP _scorematchingad_taylorApprox_currentdynparam(SEXP pfunSEXP, SEXP xSEXP, SEXP centreSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< pADFun& >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< vecd >::type x(xSEXP);
    Rcpp::traits::input_parameter< vecd >::type centre(centreSEXP);
    Rcpp::traits::input_parameter< const size_t >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(taylorApprox_currentdynparam(pfun, x, centre, order));
    return rcpp_result_gen;
END_RCPP
}
// taylorApprox
vecd taylorApprox(pADFun& pfun, vecd x, vecd centre, vecd dynparam, const size_t order);
RcppExport SEXP _scorematchingad_taylorApprox(SEXP pfunSEXP, SEXP xSEXP, SEXP centreSEXP, SEXP dynparamSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< pADFun& >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< vecd >::type x(xSEXP);
    Rcpp::traits::input_parameter< vecd >::type centre(centreSEXP);
    Rcpp::traits::input_parameter< vecd >::type dynparam(dynparamSEXP);
    Rcpp::traits::input_parameter< const size_t >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(taylorApprox(pfun, x, centre, dynparam, order));
    return rcpp_result_gen;
END_RCPP
}
// set_cppad_error_handler
void set_cppad_error_handler();
RcppExport SEXP _scorematchingad_set_cppad_error_handler() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    set_cppad_error_handler();
    return R_NilValue;
END_RCPP
}
// test_Rcpphandler
void test_Rcpphandler();
RcppExport SEXP _scorematchingad_test_Rcpphandler() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test_Rcpphandler();
    return R_NilValue;
END_RCPP
}
// fixdynamic
pADFun fixdynamic(pADFun& pfun, veca1 theta, Eigen::Matrix<int, Eigen::Dynamic, 1> isfixed);
RcppExport SEXP _scorematchingad_fixdynamic(SEXP pfunSEXP, SEXP thetaSEXP, SEXP isfixedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< pADFun& >::type pfun(pfunSEXP);
    Rcpp::traits::input_parameter< veca1 >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Eigen::Matrix<int, Eigen::Dynamic, 1> >::type isfixed(isfixedSEXP);
    rcpp_result_gen = Rcpp::wrap(fixdynamic(pfun, theta, isfixed));
    return rcpp_result_gen;
END_RCPP
}
// reembed
pADFun reembed(pADFun& uld, transform<a1type>& tran);
RcppExport SEXP _scorematchingad_reembed(SEXP uldSEXP, SEXP tranSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< pADFun& >::type uld(uldSEXP);
    Rcpp::traits::input_parameter< transform<a1type>& >::type tran(tranSEXP);
    rcpp_result_gen = Rcpp::wrap(reembed(uld, tran));
    return rcpp_result_gen;
END_RCPP
}
// tape_Jacobian
pADFun tape_Jacobian(pADFun& pfun);
RcppExport SEXP _scorematchingad_tape_Jacobian(SEXP pfunSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< pADFun& >::type pfun(pfunSEXP);
    rcpp_result_gen = Rcpp::wrap(tape_Jacobian(pfun));
    return rcpp_result_gen;
END_RCPP
}
// tape_Hessian
pADFun tape_Hessian(pADFun& pfun);
RcppExport SEXP _scorematchingad_tape_Hessian(SEXP pfunSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< pADFun& >::type pfun(pfunSEXP);
    rcpp_result_gen = Rcpp::wrap(tape_Hessian(pfun));
    return rcpp_result_gen;
END_RCPP
}
// tape_gradoffset
pADFun tape_gradoffset(pADFun& pfun);
RcppExport SEXP _scorematchingad_tape_gradoffset(SEXP pfunSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< pADFun& >::type pfun(pfunSEXP);
    rcpp_result_gen = Rcpp::wrap(tape_gradoffset(pfun));
    return rcpp_result_gen;
END_RCPP
}
// tape_logJacdet
pADFun tape_logJacdet(pADFun& pfun);
RcppExport SEXP _scorematchingad_tape_logJacdet(SEXP pfunSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< pADFun& >::type pfun(pfunSEXP);
    rcpp_result_gen = Rcpp::wrap(tape_logJacdet(pfun));
    return rcpp_result_gen;
END_RCPP
}
// tape_swap
pADFun tape_swap(pADFun& pfun);
RcppExport SEXP _scorematchingad_tape_swap(SEXP pfunSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< pADFun& >::type pfun(pfunSEXP);
    rcpp_result_gen = Rcpp::wrap(tape_swap(pfun));
    return rcpp_result_gen;
END_RCPP
}
// ptapell2
pADFun ptapell2(veca1 z_ad, veca1 theta_ad, Rcpp::XPtr<llPtr> llfXPtr, transform_a1type& tran, Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta, bool verbose);
RcppExport SEXP _scorematchingad_ptapell2(SEXP z_adSEXP, SEXP theta_adSEXP, SEXP llfXPtrSEXP, SEXP tranSEXP, SEXP fixedthetaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< veca1 >::type z_ad(z_adSEXP);
    Rcpp::traits::input_parameter< veca1 >::type theta_ad(theta_adSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<llPtr> >::type llfXPtr(llfXPtrSEXP);
    Rcpp::traits::input_parameter< transform_a1type& >::type tran(tranSEXP);
    Rcpp::traits::input_parameter< Eigen::Matrix<int, Eigen::Dynamic, 1> >::type fixedtheta(fixedthetaSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(ptapell2(z_ad, theta_ad, llfXPtr, tran, fixedtheta, verbose));
    return rcpp_result_gen;
END_RCPP
}
// getllptr
Rcpp::XPtr<llPtr> getllptr(std::string llname);
RcppExport SEXP _scorematchingad_getllptr(SEXP llnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type llname(llnameSEXP);
    rcpp_result_gen = Rcpp::wrap(getllptr(llname));
    return rcpp_result_gen;
END_RCPP
}
// tape_uld_inbuilt
pADFun tape_uld_inbuilt(std::string name, veca1 x, veca1 theta);
RcppExport SEXP _scorematchingad_tape_uld_inbuilt(SEXP nameSEXP, SEXP xSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type name(nameSEXP);
    Rcpp::traits::input_parameter< veca1 >::type x(xSEXP);
    Rcpp::traits::input_parameter< veca1 >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(tape_uld_inbuilt(name, x, theta));
    return rcpp_result_gen;
END_RCPP
}
// tapesmd
pADFun tapesmd(pADFun& uldtape, transform<a1type>& tran, manifold<a1type>& M, std::string weightname, const double& acut, bool verbose);
RcppExport SEXP _scorematchingad_tapesmd(SEXP uldtapeSEXP, SEXP tranSEXP, SEXP MSEXP, SEXP weightnameSEXP, SEXP acutSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< pADFun& >::type uldtape(uldtapeSEXP);
    Rcpp::traits::input_parameter< transform<a1type>& >::type tran(tranSEXP);
    Rcpp::traits::input_parameter< manifold<a1type>& >::type M(MSEXP);
    Rcpp::traits::input_parameter< std::string >::type weightname(weightnameSEXP);
    Rcpp::traits::input_parameter< const double& >::type acut(acutSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(tapesmd(uldtape, tran, M, weightname, acut, verbose));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_manifolds();
RcppExport SEXP _rcpp_module_boot_ADFun();

static const R_CallMethodDef CallEntries[] = {
    {"_scorematchingad_abort_recording", (DL_FUNC) &_scorematchingad_abort_recording, 0},
    {"_scorematchingad_taylorApprox_currentdynparam", (DL_FUNC) &_scorematchingad_taylorApprox_currentdynparam, 4},
    {"_scorematchingad_taylorApprox", (DL_FUNC) &_scorematchingad_taylorApprox, 5},
    {"_scorematchingad_set_cppad_error_handler", (DL_FUNC) &_scorematchingad_set_cppad_error_handler, 0},
    {"_scorematchingad_test_Rcpphandler", (DL_FUNC) &_scorematchingad_test_Rcpphandler, 0},
    {"_scorematchingad_fixdynamic", (DL_FUNC) &_scorematchingad_fixdynamic, 3},
    {"_scorematchingad_reembed", (DL_FUNC) &_scorematchingad_reembed, 2},
    {"_scorematchingad_tape_Jacobian", (DL_FUNC) &_scorematchingad_tape_Jacobian, 1},
    {"_scorematchingad_tape_Hessian", (DL_FUNC) &_scorematchingad_tape_Hessian, 1},
    {"_scorematchingad_tape_gradoffset", (DL_FUNC) &_scorematchingad_tape_gradoffset, 1},
    {"_scorematchingad_tape_logJacdet", (DL_FUNC) &_scorematchingad_tape_logJacdet, 1},
    {"_scorematchingad_tape_swap", (DL_FUNC) &_scorematchingad_tape_swap, 1},
    {"_scorematchingad_ptapell2", (DL_FUNC) &_scorematchingad_ptapell2, 6},
    {"_scorematchingad_getllptr", (DL_FUNC) &_scorematchingad_getllptr, 1},
    {"_scorematchingad_tape_uld_inbuilt", (DL_FUNC) &_scorematchingad_tape_uld_inbuilt, 3},
    {"_scorematchingad_tapesmd", (DL_FUNC) &_scorematchingad_tapesmd, 6},
    {"_rcpp_module_boot_manifolds", (DL_FUNC) &_rcpp_module_boot_manifolds, 0},
    {"_rcpp_module_boot_ADFun", (DL_FUNC) &_rcpp_module_boot_ADFun, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_scorematchingad(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

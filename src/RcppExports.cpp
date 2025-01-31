// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// c_bound
int c_bound(arma::vec pik);
RcppExport SEXP _StratifiedSampling_c_bound(SEXP pikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pik(pikSEXP);
    rcpp_result_gen = Rcpp::wrap(c_bound(pik));
    return rcpp_result_gen;
END_RCPP
}
// c_bound2
bool c_bound2(arma::vec pik);
RcppExport SEXP _StratifiedSampling_c_bound2(SEXP pikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type pik(pikSEXP);
    rcpp_result_gen = Rcpp::wrap(c_bound2(pik));
    return rcpp_result_gen;
END_RCPP
}
// calibRaking
Rcpp::List calibRaking(arma::mat Xs, arma::vec d, arma::vec total, arma::vec q, int max_iter, double tol);
RcppExport SEXP _StratifiedSampling_calibRaking(SEXP XsSEXP, SEXP dSEXP, SEXP totalSEXP, SEXP qSEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type total(totalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(calibRaking(Xs, d, total, q, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}
// gencalibRaking
Rcpp::List gencalibRaking(arma::mat Xs, arma::mat Zs, arma::vec d, arma::vec total, arma::vec q, int max_iter, double tol);
RcppExport SEXP _StratifiedSampling_gencalibRaking(SEXP XsSEXP, SEXP ZsSEXP, SEXP dSEXP, SEXP totalSEXP, SEXP qSEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xs(XsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Zs(ZsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type total(totalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(gencalibRaking(Xs, Zs, d, total, q, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}
// disj
Rcpp::IntegerMatrix disj(Rcpp::IntegerVector strata_input);
RcppExport SEXP _StratifiedSampling_disj(SEXP strata_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type strata_input(strata_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(disj(strata_input));
    return rcpp_result_gen;
END_RCPP
}
// ncat
Rcpp::NumericVector ncat(Rcpp::IntegerMatrix Xcat_input);
RcppExport SEXP _StratifiedSampling_ncat(SEXP Xcat_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type Xcat_input(Xcat_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(ncat(Xcat_input));
    return rcpp_result_gen;
END_RCPP
}
// disjMatrix
Rcpp::IntegerMatrix disjMatrix(Rcpp::IntegerMatrix strata_input);
RcppExport SEXP _StratifiedSampling_disjMatrix(SEXP strata_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type strata_input(strata_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(disjMatrix(strata_input));
    return rcpp_result_gen;
END_RCPP
}
// distUnitk
arma::vec distUnitk(arma::mat X, int k, bool tore, double toreBound);
RcppExport SEXP _StratifiedSampling_distUnitk(SEXP XSEXP, SEXP kSEXP, SEXP toreSEXP, SEXP toreBoundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type tore(toreSEXP);
    Rcpp::traits::input_parameter< double >::type toreBound(toreBoundSEXP);
    rcpp_result_gen = Rcpp::wrap(distUnitk(X, k, tore, toreBound));
    return rcpp_result_gen;
END_RCPP
}
// ffphase
Rcpp::NumericVector ffphase(arma::mat Xbal, arma::vec prob, bool order);
RcppExport SEXP _StratifiedSampling_ffphase(SEXP XbalSEXP, SEXP probSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xbal(XbalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prob(probSEXP);
    Rcpp::traits::input_parameter< bool >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(ffphase(Xbal, prob, order));
    return rcpp_result_gen;
END_RCPP
}
// inclprob
arma::vec inclprob(arma::vec& x, const double& n);
RcppExport SEXP _StratifiedSampling_inclprob(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(inclprob(x, n));
    return rcpp_result_gen;
END_RCPP
}
// qfromw
NumericMatrix qfromw(NumericVector& w, const int& n);
RcppExport SEXP _StratifiedSampling_qfromw(SEXP wSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(qfromw(w, n));
    return rcpp_result_gen;
END_RCPP
}
// sfromq
IntegerVector sfromq(const NumericMatrix& q);
RcppExport SEXP _StratifiedSampling_sfromq(SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(sfromq(q));
    return rcpp_result_gen;
END_RCPP
}
// pikfromq
NumericVector pikfromq(NumericMatrix& q);
RcppExport SEXP _StratifiedSampling_pikfromq(SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(pikfromq(q));
    return rcpp_result_gen;
END_RCPP
}
// piktfrompik
NumericVector piktfrompik(NumericVector& pik, int max_iter, double tol);
RcppExport SEXP _StratifiedSampling_piktfrompik(SEXP pikSEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(piktfrompik(pik, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}
// cps
IntegerVector cps(NumericVector& pik, double eps);
RcppExport SEXP _StratifiedSampling_cps(SEXP pikSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(cps(pik, eps));
    return rcpp_result_gen;
END_RCPP
}
// maxentpi2
NumericMatrix maxentpi2(NumericVector pikr);
RcppExport SEXP _StratifiedSampling_maxentpi2(SEXP pikrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pikr(pikrSEXP);
    rcpp_result_gen = Rcpp::wrap(maxentpi2(pikr));
    return rcpp_result_gen;
END_RCPP
}
// osod
IntegerVector osod(NumericVector pikr, bool full);
RcppExport SEXP _StratifiedSampling_osod(SEXP pikrSEXP, SEXP fullSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pikr(pikrSEXP);
    Rcpp::traits::input_parameter< bool >::type full(fullSEXP);
    rcpp_result_gen = Rcpp::wrap(osod(pikr, full));
    return rcpp_result_gen;
END_RCPP
}
// vEst
arma::mat vEst(arma::mat Xauxs, arma::vec piks, arma::vec ys);
RcppExport SEXP _StratifiedSampling_vEst(SEXP XauxsSEXP, SEXP piksSEXP, SEXP ysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xauxs(XauxsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type piks(piksSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ys(ysSEXP);
    rcpp_result_gen = Rcpp::wrap(vEst(Xauxs, piks, ys));
    return rcpp_result_gen;
END_RCPP
}
// vDBS
double vDBS(arma::mat Xauxs, arma::mat Xspreads, arma::vec piks, arma::vec ys);
RcppExport SEXP _StratifiedSampling_vDBS(SEXP XauxsSEXP, SEXP XspreadsSEXP, SEXP piksSEXP, SEXP ysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xauxs(XauxsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xspreads(XspreadsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type piks(piksSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ys(ysSEXP);
    rcpp_result_gen = Rcpp::wrap(vDBS(Xauxs, Xspreads, piks, ys));
    return rcpp_result_gen;
END_RCPP
}
// vApp
arma::mat vApp(arma::mat Xaux, arma::vec pik, arma::vec y);
RcppExport SEXP _StratifiedSampling_vApp(SEXP XauxSEXP, SEXP pikSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xaux(XauxSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(vApp(Xaux, pik, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StratifiedSampling_c_bound", (DL_FUNC) &_StratifiedSampling_c_bound, 1},
    {"_StratifiedSampling_c_bound2", (DL_FUNC) &_StratifiedSampling_c_bound2, 1},
    {"_StratifiedSampling_calibRaking", (DL_FUNC) &_StratifiedSampling_calibRaking, 6},
    {"_StratifiedSampling_gencalibRaking", (DL_FUNC) &_StratifiedSampling_gencalibRaking, 7},
    {"_StratifiedSampling_disj", (DL_FUNC) &_StratifiedSampling_disj, 1},
    {"_StratifiedSampling_ncat", (DL_FUNC) &_StratifiedSampling_ncat, 1},
    {"_StratifiedSampling_disjMatrix", (DL_FUNC) &_StratifiedSampling_disjMatrix, 1},
    {"_StratifiedSampling_distUnitk", (DL_FUNC) &_StratifiedSampling_distUnitk, 4},
    {"_StratifiedSampling_ffphase", (DL_FUNC) &_StratifiedSampling_ffphase, 3},
    {"_StratifiedSampling_inclprob", (DL_FUNC) &_StratifiedSampling_inclprob, 2},
    {"_StratifiedSampling_qfromw", (DL_FUNC) &_StratifiedSampling_qfromw, 2},
    {"_StratifiedSampling_sfromq", (DL_FUNC) &_StratifiedSampling_sfromq, 1},
    {"_StratifiedSampling_pikfromq", (DL_FUNC) &_StratifiedSampling_pikfromq, 1},
    {"_StratifiedSampling_piktfrompik", (DL_FUNC) &_StratifiedSampling_piktfrompik, 3},
    {"_StratifiedSampling_cps", (DL_FUNC) &_StratifiedSampling_cps, 2},
    {"_StratifiedSampling_maxentpi2", (DL_FUNC) &_StratifiedSampling_maxentpi2, 1},
    {"_StratifiedSampling_osod", (DL_FUNC) &_StratifiedSampling_osod, 2},
    {"_StratifiedSampling_vEst", (DL_FUNC) &_StratifiedSampling_vEst, 3},
    {"_StratifiedSampling_vDBS", (DL_FUNC) &_StratifiedSampling_vDBS, 4},
    {"_StratifiedSampling_vApp", (DL_FUNC) &_StratifiedSampling_vApp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_StratifiedSampling(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

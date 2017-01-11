// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// HiddenVarRidgei
double HiddenVarRidgei(int ii, Rcpp::NumericMatrix tX, double aRand, double bRand, double cSigma, double dSigma, int maxiter, arma::colvec varSel);
RcppExport SEXP ShrinkNet_HiddenVarRidgei(SEXP iiSEXP, SEXP tXSEXP, SEXP aRandSEXP, SEXP bRandSEXP, SEXP cSigmaSEXP, SEXP dSigmaSEXP, SEXP maxiterSEXP, SEXP varSelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ii(iiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX(tXSEXP);
    Rcpp::traits::input_parameter< double >::type aRand(aRandSEXP);
    Rcpp::traits::input_parameter< double >::type bRand(bRandSEXP);
    Rcpp::traits::input_parameter< double >::type cSigma(cSigmaSEXP);
    Rcpp::traits::input_parameter< double >::type dSigma(dSigmaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type varSel(varSelSEXP);
    rcpp_result_gen = Rcpp::wrap(HiddenVarRidgei(ii, tX, aRand, bRand, cSigma, dSigma, maxiter, varSel));
    return rcpp_result_gen;
END_RCPP
}
// HiddenEdgeBFprime
arma::colvec HiddenEdgeBFprime(Rcpp::NumericVector idx, Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX);
RcppExport SEXP ShrinkNet_HiddenEdgeBFprime(SEXP idxSEXP, SEXP thematSEXP, SEXP tXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type themat(thematSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX(tXSEXP);
    rcpp_result_gen = Rcpp::wrap(HiddenEdgeBFprime(idx, themat, tX));
    return rcpp_result_gen;
END_RCPP
}
// HiddenEstimatep0
double HiddenEstimatep0(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX);
RcppExport SEXP ShrinkNet_HiddenEstimatep0(SEXP thematSEXP, SEXP tXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type themat(thematSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX(tXSEXP);
    rcpp_result_gen = Rcpp::wrap(HiddenEstimatep0(themat, tX));
    return rcpp_result_gen;
END_RCPP
}
// HiddenEdgeSelection
Rcpp::List HiddenEdgeSelection(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX, double p0, double lfdrcut, int maxNbEdges);
RcppExport SEXP ShrinkNet_HiddenEdgeSelection(SEXP thematSEXP, SEXP tXSEXP, SEXP p0SEXP, SEXP lfdrcutSEXP, SEXP maxNbEdgesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type themat(thematSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX(tXSEXP);
    Rcpp::traits::input_parameter< double >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< double >::type lfdrcut(lfdrcutSEXP);
    Rcpp::traits::input_parameter< int >::type maxNbEdges(maxNbEdgesSEXP);
    rcpp_result_gen = Rcpp::wrap(HiddenEdgeSelection(themat, tX, p0, lfdrcut, maxNbEdges));
    return rcpp_result_gen;
END_RCPP
}
// HiddenVarRidgeiGetKappa
arma::mat HiddenVarRidgeiGetKappa(int ii, Rcpp::List SVDs, arma::mat tX, double aRand, double bRand, arma::colvec bRandStarInit, arma::colvec dSigmaStarInit);
RcppExport SEXP ShrinkNet_HiddenVarRidgeiGetKappa(SEXP iiSEXP, SEXP SVDsSEXP, SEXP tXSEXP, SEXP aRandSEXP, SEXP bRandSEXP, SEXP bRandStarInitSEXP, SEXP dSigmaStarInitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ii(iiSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type SVDs(SVDsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type tX(tXSEXP);
    Rcpp::traits::input_parameter< double >::type aRand(aRandSEXP);
    Rcpp::traits::input_parameter< double >::type bRand(bRandSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type bRandStarInit(bRandStarInitSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type dSigmaStarInit(dSigmaStarInitSEXP);
    rcpp_result_gen = Rcpp::wrap(HiddenVarRidgeiGetKappa(ii, SVDs, tX, aRand, bRand, bRandStarInit, dSigmaStarInit));
    return rcpp_result_gen;
END_RCPP
}
// HiddenVarAlgo
Rcpp::List HiddenVarAlgo(Rcpp::List SVDs, arma::mat tX, double aRand, double bRand, int maxiter, int globalShrink, double tol, bool verbose);
RcppExport SEXP ShrinkNet_HiddenVarAlgo(SEXP SVDsSEXP, SEXP tXSEXP, SEXP aRandSEXP, SEXP bRandSEXP, SEXP maxiterSEXP, SEXP globalShrinkSEXP, SEXP tolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type SVDs(SVDsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type tX(tXSEXP);
    Rcpp::traits::input_parameter< double >::type aRand(aRandSEXP);
    Rcpp::traits::input_parameter< double >::type bRand(bRandSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< int >::type globalShrink(globalShrinkSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(HiddenVarAlgo(SVDs, tX, aRand, bRand, maxiter, globalShrink, tol, verbose));
    return rcpp_result_gen;
END_RCPP
}
// getSVD
Rcpp::List getSVD(int ii, Rcpp::NumericMatrix tX);
RcppExport SEXP ShrinkNet_getSVD(SEXP iiSEXP, SEXP tXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ii(iiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX(tXSEXP);
    rcpp_result_gen = Rcpp::wrap(getSVD(ii, tX));
    return rcpp_result_gen;
END_RCPP
}

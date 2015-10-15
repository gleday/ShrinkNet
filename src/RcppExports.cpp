// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// varRidgei
double varRidgei(int ii, Rcpp::NumericMatrix tX, double aRand, double bRand, double cSigma, double dSigma, int maxiter, arma::colvec varSel);
RcppExport SEXP ShrinkNet_varRidgei(SEXP iiSEXP, SEXP tXSEXP, SEXP aRandSEXP, SEXP bRandSEXP, SEXP cSigmaSEXP, SEXP dSigmaSEXP, SEXP maxiterSEXP, SEXP varSelSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type ii(iiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX(tXSEXP);
    Rcpp::traits::input_parameter< double >::type aRand(aRandSEXP);
    Rcpp::traits::input_parameter< double >::type bRand(bRandSEXP);
    Rcpp::traits::input_parameter< double >::type cSigma(cSigmaSEXP);
    Rcpp::traits::input_parameter< double >::type dSigma(dSigmaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type varSel(varSelSEXP);
    __result = Rcpp::wrap(varRidgei(ii, tX, aRand, bRand, cSigma, dSigma, maxiter, varSel));
    return __result;
END_RCPP
}
// estimatep0
double estimatep0(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX, int maxedges);
RcppExport SEXP ShrinkNet_estimatep0(SEXP thematSEXP, SEXP tXSEXP, SEXP maxedgesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type themat(thematSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX(tXSEXP);
    Rcpp::traits::input_parameter< int >::type maxedges(maxedgesSEXP);
    __result = Rcpp::wrap(estimatep0(themat, tX, maxedges));
    return __result;
END_RCPP
}
// edgeSelection
arma::mat edgeSelection(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX, double p0, int maxedges, double lfdrcut);
RcppExport SEXP ShrinkNet_edgeSelection(SEXP thematSEXP, SEXP tXSEXP, SEXP p0SEXP, SEXP maxedgesSEXP, SEXP lfdrcutSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type themat(thematSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX(tXSEXP);
    Rcpp::traits::input_parameter< double >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< int >::type maxedges(maxedgesSEXP);
    Rcpp::traits::input_parameter< double >::type lfdrcut(lfdrcutSEXP);
    __result = Rcpp::wrap(edgeSelection(themat, tX, p0, maxedges, lfdrcut));
    return __result;
END_RCPP
}
// varRidgeiOneIter
Rcpp::List varRidgeiOneIter(int ii, Rcpp::List SVDs, double aRand, double bRand, arma::colvec bRandStarInit, arma::colvec dSigmaStarInit);
RcppExport SEXP ShrinkNet_varRidgeiOneIter(SEXP iiSEXP, SEXP SVDsSEXP, SEXP aRandSEXP, SEXP bRandSEXP, SEXP bRandStarInitSEXP, SEXP dSigmaStarInitSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type ii(iiSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type SVDs(SVDsSEXP);
    Rcpp::traits::input_parameter< double >::type aRand(aRandSEXP);
    Rcpp::traits::input_parameter< double >::type bRand(bRandSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type bRandStarInit(bRandStarInitSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type dSigmaStarInit(dSigmaStarInitSEXP);
    __result = Rcpp::wrap(varRidgeiOneIter(ii, SVDs, aRand, bRand, bRandStarInit, dSigmaStarInit));
    return __result;
END_RCPP
}
// mydigamma
arma::colvec mydigamma(arma::colvec vec);
RcppExport SEXP ShrinkNet_mydigamma(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::colvec >::type vec(vecSEXP);
    __result = Rcpp::wrap(mydigamma(vec));
    return __result;
END_RCPP
}
// varAlgo
Rcpp::List varAlgo(Rcpp::List SVDs, double aRand, double bRand, int maxiter, int globalShrink, double tol);
RcppExport SEXP ShrinkNet_varAlgo(SEXP SVDsSEXP, SEXP aRandSEXP, SEXP bRandSEXP, SEXP maxiterSEXP, SEXP globalShrinkSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::List >::type SVDs(SVDsSEXP);
    Rcpp::traits::input_parameter< double >::type aRand(aRandSEXP);
    Rcpp::traits::input_parameter< double >::type bRand(bRandSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< int >::type globalShrink(globalShrinkSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    __result = Rcpp::wrap(varAlgo(SVDs, aRand, bRand, maxiter, globalShrink, tol));
    return __result;
END_RCPP
}
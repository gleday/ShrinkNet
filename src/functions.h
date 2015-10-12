#ifndef FUNCTIONS_H
#define FUNCTIONS_H

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rcpp.h>

// Variational algorithm for local shrinkage - update only once variational parameters
Rcpp::List varRidgeiOneIter(int ii,  Rcpp::List SVDs, double aRand, double bRand, arma::colvec bRandStarInit, arma::colvec dSigmaStarInit);

// Variational algorithm for global-local shrinkage
Rcpp::List varAlgo(Rcpp::List SVDs, double aRand, double bRand, int maxiter, int globalShrink);

// Variational algorithm for local shrinkage
double varRidgei(int ii,  Rcpp::NumericMatrix tX, double aRand, double bRand, double cSigma, double dSigma, int maxiter, arma::colvec varSel);

// Estimate proportion of null hypothesis
double estimatep0(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX, int maxedges);

// Edge selection using Bayesian analogue of local false discovery rate
arma::mat edgeSelection(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX, double p0, int maxedges, double lfdrcut);

// digamma function with input vector
arma::colvec mydigamma(arma::colvec vec);

#endif // FUNCTIONS_H

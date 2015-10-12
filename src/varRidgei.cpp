#include "functions.h"

// [[Rcpp::export]]
double varRidgei(int ii,  Rcpp::NumericMatrix tX, double aRand, double bRand, double cSigma, double dSigma, int maxiter, arma::colvec varSel){

  //Rcpp::Rcout << "ii= " << ii << std::endl;

  using arma::trans;

  // Data
  arma::mat tXbis(tX.begin(),tX.nrow(),tX.ncol(),false);
  arma::colvec myy = trans(tXbis.row(ii-1));
  arma::mat tXbis2 = tXbis.rows(arma::find(arma::linspace(1,tXbis.n_rows, tXbis.n_rows)!=ii));

  // Variable subset
  arma::mat myX = arma::ones(myy.n_elem);
  if(sum(varSel)>0){
    arma::colvec tp=varSel;
    arma::uvec idxs = arma::find(tp);
    myX.insert_cols(myX.n_cols, trans(tXbis2.rows(idxs)));
  }
  arma::mat XTX = myX.t() * myX;

  // Initialization
  int then = myX.n_rows;
  int thep = myX.n_cols;
  double aRandStar = aRand+0.5*thep;
  double cSigmaStar = cSigma+0.5*then+0.5*thep;
  arma::colvec ratesRand(maxiter+1);
  arma::colvec ratesSig(maxiter+1);
  arma::colvec L(maxiter+1);
  ratesRand(0) = bRand;
  ratesSig(0) = dSigma;
  int ct = 1;
  bool mybool = true;
  double expTau, expSig, Lrand, Lsig, theplus;
  float Ldiff;
  arma::colvec postMean(thep);
  arma::mat postSigma(thep,thep);

  // Algo
  while(mybool){

    // Expectations
    expTau = aRandStar/ratesRand(ct-1);
    expSig = cSigmaStar/ratesSig(ct-1);

    // Update sigma
    postSigma = arma::inv(expSig*(XTX + expTau*arma::eye<arma::mat>(thep, thep)));

    // Update beta
    postMean = expSig*postSigma*trans(myX)*myy;

    // Update Gamma rate for \sigma_{-2}
    ratesSig(ct) = dSigma + 0.5*(arma::as_scalar(trans(myy-myX*postMean)*(myy-myX*postMean))+arma::trace(XTX*postSigma)) + 0.5*expTau*(arma::as_scalar(trans(postMean)*(postMean)) + arma::trace(postSigma));

    // Update Gamma rate for \tau_{-2}
    ratesRand(ct) = bRand + 0.5*expSig*(arma::as_scalar(trans(postMean)*(postMean)) + arma::trace(postSigma));

    // Lower bound marginal likelihood
    Lrand = aRand*log(bRand)-aRandStar*log(ratesRand(ct))+lgamma(aRandStar)-lgamma(aRand);
    Lsig = cSigma*log(dSigma)-cSigmaStar*log(ratesSig(ct))+lgamma(cSigmaStar)-lgamma(cSigma);
    theplus = 0.5*(expSig)*(expTau)*(arma::as_scalar(trans(postMean)*postMean) + arma::trace(postSigma));
    try{
      L(ct) = 0.5*thep - 0.5*then*log(2*arma::datum::pi) + 0.5*sum(arma::log(arma::eig_sym(postSigma))) + Lrand + Lsig + theplus;
    }
    catch(...){
      L(ct) = arma::datum::nan;
    }

    // Monitor convergence
    if(ct==maxiter){
      mybool = FALSE;
    }
    if(ct>2){
      Ldiff = std::abs(L(ct)-L(ct-1));
      if(Ldiff<0.001){
        mybool = FALSE;
      }
    }
    ct++;
  }

  return L(ct-1);
}

#include "functions.h"

// [[Rcpp::export]]
Rcpp::List varRidgeiOneIter(int ii,  Rcpp::List SVDs, double aRand, double bRand, arma::colvec bRandStarInit, arma::colvec dSigmaStarInit){

  //Rcpp::Rcout << "ii= " << ii << std::endl;

  using arma::trans;
  double cSigma = 0.001;
  double dSigma = 0.001;

  // Data
  Rcpp::List tplist = Rcpp::as<Rcpp::List>(SVDs[ii-1]);
  arma::colvec myy = Rcpp::as<arma::colvec>(tplist["myy"]);
  arma::mat myX = Rcpp::as<arma::mat>(tplist["myX"]);
  arma::mat XTX = myX.t() * myX;
  int then = myX.n_rows;
  int thep = myX.n_cols;

  // Update posterior shape parameters
  double aRandStar = aRand + 0.5*thep;
  double cSigmaStar = cSigma + 0.5*then + 0.5*thep;

  // Expectations
  double expTau, expSig;
  if(all(bRandStarInit==0)){
    expTau = aRandStar/bRand;
    expSig = cSigmaStar/dSigma;
  }else{
    expTau = aRandStar/bRandStarInit(ii-1);
    expSig = cSigmaStar/dSigmaStarInit(ii-1);
  }

  // Update sigma
  //postSigma <- (1/expSig)*SVDs[[ii]]$v %*% solve(SVDs[[ii]]$FTF + expTau*diag(thep)) %*% t(SVDs[[ii]]$v)
  arma::mat postSigma = arma::inv(expSig*(XTX + expTau*arma::eye<arma::mat>(thep, thep)));

  // Update beta
  arma::colvec postMean = expSig*postSigma*trans(myX)*myy;

  // Update posterior rate for sigma^2
  double dSigmaStar = dSigma + 0.5*(arma::as_scalar(trans(myy-myX*postMean)*(myy-myX*postMean))+arma::trace(XTX*postSigma)) + 0.5*expTau*(arma::as_scalar(trans(postMean)*(postMean)) + arma::trace(postSigma));

  // Update posterior rates for tau^2
  double bRandStar = bRand + 0.5*expSig*(arma::as_scalar(trans(postMean)*(postMean)) + arma::trace(postSigma));

  //    if(ii==1){
  //      Rcpp::Rcout << "expTau = " << expTau << std::endl;
  //      Rcpp::Rcout << "bRand = " << bRand << std::endl;
  //      Rcpp::Rcout << "expSig = " << expSig << std::endl;
  //      Rcpp::Rcout << "arma::as_scalar(trans(postMean)*(postMean)) = " << arma::as_scalar(trans(postMean)*(postMean)) << std::endl;
  //      Rcpp::Rcout << "arma::trace(postSigma) = " << arma::trace(postSigma) << std::endl;
  //    }

  // Lower bound marginal likelihood
  double Lrand = aRand*log(bRand)-aRandStar*log(bRandStar)+lgamma(aRandStar)-lgamma(aRand);
  double Lsig = cSigma*log(dSigma)-cSigmaStar*log(dSigmaStar)+lgamma(cSigmaStar)-lgamma(cSigma);
  double theplus = 0.5*(expSig)*(expTau)*(arma::as_scalar(trans(postMean)*postMean) + arma::trace(postSigma));
  double L;
  try{
    L = 0.5*thep - 0.5*then*log(2*arma::datum::pi) + 0.5*sum(arma::log(arma::eig_sym(postSigma))) + Lrand + Lsig + theplus;
  }
  catch(...){
    L = arma::datum::nan;
  }
  //Rcpp::Rcout << "Debug X5 " << std::endl;
  // Output
  arma::colvec postSd = sqrt(arma::diagvec(postSigma));
  arma::colvec ratio = arma::abs(postMean)/postSd;
  arma::mat postBeta;
  postBeta.insert_cols(postBeta.n_cols, postMean);
  postBeta.insert_cols(postBeta.n_cols, postSd);
  postBeta.insert_cols(postBeta.n_cols, ratio);
  Rcpp::NumericVector priorRand(2);
  priorRand(0) = aRand;
  priorRand(1) = bRand;
  Rcpp::NumericVector priorSig(2);
  priorSig(0) = cSigma;
  priorSig(1) = dSigma;
  Rcpp::NumericVector postRand(2);
  postRand(0) = aRandStar;
  postRand(1) = bRandStar;
  Rcpp::NumericVector postSig(2);
  postSig(0) = cSigmaStar;
  postSig(1) = dSigmaStar;
  //Rcpp::Rcout << "Debug X6 " << std::endl;
  return Rcpp::List::create(Rcpp::Named("postBeta") = postBeta, Rcpp::Named("L") = L, Rcpp::Named("priorRand") = priorRand, Rcpp::Named("priorSig") = priorSig, Rcpp::Named("postRand") = postRand, Rcpp::Named("postSig") = postSig);
}

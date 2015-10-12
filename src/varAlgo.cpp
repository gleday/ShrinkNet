#include "functions.h"

// [[Rcpp::export]]
Rcpp::List varAlgo(Rcpp::List SVDs, double aRand, double bRand, int maxiter, int globalShrink){

  // Initialization
  arma::mat allmargs(maxiter+1, SVDs.size());
  arma::mat parTau(maxiter+1, 2);
  parTau(0,0) = aRand;
  parTau(0,1) = bRand;
  int ct = 0;
  bool mybool = true;
  arma::colvec allaRandStar(SVDs.size());
  arma::colvec allbRandStar(SVDs.size());
  arma::colvec alldSigmaStar(SVDs.size());
  allaRandStar.zeros();
  allbRandStar.zeros();
  alldSigmaStar.zeros();
  arma::colvec allbRandStarnew(SVDs.size());
  allbRandStarnew.zeros();
  arma::colvec alldSigmaStarnew(SVDs.size());
  alldSigmaStarnew.zeros();
  Rcpp::List tplist;
  arma::colvec tpvals, tpvals2;
  arma::mat tpmat;
  arma::mat matThres(SVDs.size(), SVDs.size());
  matThres.zeros();
  arma::uvec idxs;
  arma::colvec tpmat2(SVDs.size());
  tpmat2.zeros();
  arma::colvec a1(30);
  arma::colvec b1(30);

  // Algo
  while(mybool){
    Rcpp::Rcout << "iteration " << ct+1 << std::endl;

    // Fit all models
    for(int j=0; j<SVDs.size(); j++){
      tplist = Rcpp::as<Rcpp::List>(varRidgeiOneIter(j+1, SVDs, parTau(ct,0), parTau(ct,1), allbRandStar, alldSigmaStar));
      tpvals = Rcpp::as<arma::colvec>(tplist["postRand"]);
      tpvals2 = Rcpp::as<arma::colvec>(tplist["postSig"]);
      allaRandStar(j) = tpvals(0);
      allbRandStarnew(j) = tpvals(1);
      alldSigmaStarnew(j) = tpvals2(1);
      allmargs(ct,j) = Rcpp::as<double>(tplist["L"]);
      tpmat = Rcpp::as<arma::mat>(tplist["postBeta"]);
      idxs = arma::find(arma::linspace(1,SVDs.size(), SVDs.size())!=(j+1));
      tpmat2.zeros();
      tpmat2.elem(idxs) = tpmat.col(2);
      matThres.col(j) = tpmat2;
    }
    allbRandStar = allbRandStarnew;
    allbRandStarnew.zeros();
    alldSigmaStar = alldSigmaStarnew;
    alldSigmaStarnew.zeros();

    // Variational Empirical Bayes using fixed-point iteration as in Valpola and Honkela (2006)
    if(globalShrink==1){
      a1.zeros();
      b1.zeros();
      a1(0) = parTau(ct,0);
      b1(0) = parTau(ct,1);
      int cpt = 1;
      bool mybool2 = true;
      double tp = (allaRandStar.n_elem/sum(allaRandStar/allbRandStar));
      double tp2 = mean(log(allbRandStar) - mydigamma(allaRandStar))-log(tp);
      while(mybool2){
        a1(cpt) = a1(cpt-1) + 0.5*(1/(R::digamma(a1(cpt-1))-log(a1(cpt-1)))) + 0.5*(1/tp2);
        b1(cpt) = a1(cpt)*tp;
        if(cpt==29){
          mybool2 = false;
        }else{
          cpt++;
        }//end if
      }//end while
      parTau(ct+1,0) = a1(cpt);
      parTau(ct+1,1) = b1(cpt);
    }
    // Variational Empirical Bayes using approximate analytical solution as in Leday et al (2015)
    if(globalShrink==2){
      double tp = sum(allaRandStar/allbRandStar);
      double tp2 = mean(log(allbRandStar)-mydigamma(allaRandStar));
      parTau(ct+1,0) = 0.5*(1/(log(tp)+tp2-log(SVDs.size())));
      parTau(ct+1,1) = parTau(ct+1,0)*SVDs.size()*(1/tp);
    }//end if

    // Monitor convergence
    double maxDiffML;
    if(ct>2){
      // Check relative increase in total log-ML
      maxDiffML = max(abs((allmargs.row(ct)-allmargs.row(ct-1))/allmargs.row(ct-1)));
      if(maxDiffML<0.0001){
        mybool = false;
      }//end if
    }//end if
    if(ct==(maxiter)){
      mybool = false;
    }else{
      ct++;
    }
  }//end while
  Rcpp::Rcout << "DONE" << std::endl;

  // Symmetrize matThres
  matThres = (matThres + matThres.t())/2;

  return Rcpp::List::create(Rcpp::Named("matThres") = matThres, Rcpp::Named("parTau") = parTau.rows(0,ct-1), Rcpp::Named("allmargs") = allmargs.rows(0,ct-1));
}//end varAlgo

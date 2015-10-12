#include "functions.h"

// [[Rcpp::export]]
double estimatep0(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX, int maxedges){

  using arma::trans;

  arma::mat mymat(themat.begin(),themat.nrow(),themat.ncol(),false);
  arma::mat mymatL = arma::trimatl(mymat);
  arma::colvec uniquevals = arma::unique(mymatL);
  arma::colvec allvals = arma::sort(uniquevals.elem(arma::find(uniquevals!=0)), "descend");

  // Algo
  arma::mat tempGraph = arma::zeros(themat.nrow(),themat.ncol());
  arma::colvec allML(themat.nrow());
  for(int j=1; j<=tX.nrow(); j++){
    allML(j-1) = varRidgei(j, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, arma::zeros(tX.nrow()-1));
  }
  double logML01, logML02, logML11, logML12;
  double L1, L2;
  arma::mat logBFs(maxedges,2);
  logBFs.zeros();
  int cpt = 0;
  bool mybool = true;
  while(mybool){

    // Consider new edge
    arma::uvec newedge = arma::find(mymat==allvals(cpt))/themat.nrow();
    tempGraph(newedge(1), newedge(0)) = 1;
    tempGraph(newedge(0), newedge(1)) = 1;
    logML01 = allML(newedge(1));
    logML02 = allML(newedge(0));

    // Obtain ML for the two regression equations that include to this edge
    arma::colvec varSel1 = tempGraph.col(newedge(1));
    arma::colvec varSel2 = tempGraph.col(newedge(0));
    logML11 = varRidgei(newedge(1)+1, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel1.elem(arma::find(arma::linspace(1,tX.nrow(), tX.nrow())!=(newedge(1)+1))));
    logML12 = varRidgei(newedge(0)+1, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel2.elem(arma::find(arma::linspace(1,tX.nrow(), tX.nrow())!=(newedge(0)+1))));
    allML(newedge(1)) = logML11;
    allML(newedge(0)) = logML12;

    // logarithm of Bayes factors
    logBFs(cpt, 0) = logML11-logML01;
    logBFs(cpt, 1) = logML12-logML02;

    arma::uvec idx1 = arma::find(logBFs.col(0)>0);
    arma::uvec idx2 = arma::find(logBFs.col(1)>0);
    double p0 = 1-(((double)idx1.n_elem+(double)idx2.n_elem)/((double)tX.nrow()*((double)tX.nrow()-1)));

    Rcpp::Rcout << "cpt = " << cpt << " - " << p0 << std::endl;

    // Convergence
    if(cpt==(maxedges-1)){
      mybool = false;
    }else{
      cpt++;
    }
  }

  // OUTPUT
  arma::uvec idx1 = arma::find(logBFs.col(0)>0);
  arma::uvec idx2 = arma::find(logBFs.col(1)>0);
  double p0 = 1-(((double)idx1.n_elem+(double)idx2.n_elem)/((double)tX.nrow()*((double)tX.nrow()-1)));

  //return List::create(Named("p0") = p0, Named("logBFs") = logBFs.rows(0,cpt));
  return p0;
}

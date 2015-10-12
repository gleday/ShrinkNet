#include "functions.h"

// [[Rcpp::export]]
arma::mat edgeSelection(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX, double p0, int maxedges, double lfdrcut){

  using arma::trans;

  arma::mat mymat(themat.begin(),themat.nrow(),themat.ncol(),false);
  arma::mat mymatL = arma::trimatl(mymat);
  arma::colvec uniquevals = arma::unique(mymatL);
  arma::colvec allvals = arma::sort(uniquevals.elem(arma::find(uniquevals!=0)), "descend");

  // Algo
  arma::mat tempGraph = arma::zeros(themat.nrow(),themat.ncol());
  arma::mat myGraph = arma::zeros(themat.nrow(),themat.ncol());
  arma::colvec allML(themat.nrow());
  for(int j=1; j<=tX.nrow(); j++){
    allML(j-1) = varRidgei(j, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, arma::zeros(tX.nrow()-1));
  }
  double logML01, logML02, logML11, logML12;
  double L1, L2;
  int cpt = 0;
  int ctstop = 0;
  bool mybool = true;
  double maxBF, minlfdr;
  arma::colvec logBFs(2);
  logBFs.zeros();
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

    // logarithm of Bayes factors
    logBFs(0) = logML11-logML01;
    logBFs(1) = logML12-logML02;
    maxBF = exp(max(logBFs));
    minlfdr = p0/(maxBF*(1-p0)+p0);

    // Thresholding - edge selection
    if(minlfdr<=lfdrcut){
      myGraph(newedge(1), newedge(0)) = 1;
      myGraph(newedge(0), newedge(1)) = 1;
      allML(newedge(1)) = logML11;
      allML(newedge(0)) = logML12;
      ctstop = 0;
    }else{
      tempGraph(newedge(1), newedge(0)) = 0;
      tempGraph(newedge(0), newedge(1)) = 0;
      ctstop++;
    }

    // Convergence
    if(cpt==(maxedges-1) || ctstop==100){
      mybool = false;
    }else{
      cpt++;
    }
  }

  return myGraph;
}

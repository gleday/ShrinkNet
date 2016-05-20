//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rcpp.h>

using namespace arma;

// [[Rcpp::export]]
double HiddenVarRidgei(int ii,  Rcpp::NumericMatrix tX, double aRand, double bRand, double cSigma, double dSigma, int maxiter, arma::colvec varSel){
  
  //Rcpp::Rcout << "ii= " << ii << std::endl;
  
  // Data
  mat tXbis(tX.begin(),tX.nrow(),tX.ncol(),false);
  colvec myy = trans(tXbis.row(ii-1));
  mat tXbis2 = tXbis.rows(find(linspace(1,tXbis.n_rows, tXbis.n_rows)!=ii));
  
  // Variable subset
  mat myX = ones(myy.n_elem);
  if(sum(varSel)>0){
    colvec tp=varSel;
    uvec idxs = find(tp);
    myX.insert_cols(myX.n_cols, trans(tXbis2.rows(idxs)));
  }
  int then = myX.n_rows;
  int thep = myX.n_cols;
  
  double aRandStar = aRand+0.5*thep;
  double cSigmaStar = cSigma+0.5*then+0.5*thep;
  colvec ratesRand(maxiter+1);
  colvec ratesSig(maxiter+1);
  colvec L(maxiter+1);
  ratesRand(0) = bRand;
  ratesSig(0) = dSigma;
  int ct = 1;
  bool mybool = true;
  double expTau, expSig, Lrand, Lsig, theplus;
  float Ldiff;
  
  if(thep<4){
    // Initialization
    mat XTX = myX.t() * myX;
    colvec postMean(thep);
    mat postSigma(thep,thep);
    
    // Algo
    while(mybool){
      
      // Expectations
      expTau = aRandStar/ratesRand(ct-1);
      expSig = cSigmaStar/ratesSig(ct-1);
      
      // Update sigma
      postSigma = inv(expSig*(XTX + expTau*eye<mat>(thep, thep)));
      
      // Update beta
      postMean = expSig*postSigma*trans(myX)*myy;
      
      // Update Gamma rate for \sigma_{-2}
      ratesSig(ct) = dSigma + 0.5*(as_scalar(trans(myy-myX*postMean)*(myy-myX*postMean))+trace(XTX*postSigma)) + 0.5*expTau*(as_scalar(trans(postMean)*(postMean)) + trace(postSigma));
      
      // Update Gamma rate for \tau_{-2}
      ratesRand(ct) = bRand + 0.5*expSig*(as_scalar(trans(postMean)*(postMean)) + trace(postSigma));
      
      // Lower bound marginal likelihood
      Lrand = aRand*log(bRand)-aRandStar*log(ratesRand(ct))+lgamma(aRandStar)-lgamma(aRand);
      Lsig = cSigma*log(dSigma)-cSigmaStar*log(ratesSig(ct))+lgamma(cSigmaStar)-lgamma(cSigma);
      theplus = 0.5*(expSig)*(expTau)*(as_scalar(trans(postMean)*postMean) + trace(postSigma));
      try{
        L(ct) = 0.5*thep - 0.5*then*log(2*datum::pi) + 0.5*sum(log(eig_sym(postSigma))) + Lrand + Lsig + theplus;
      }
      catch(...){
        L(ct) = datum::nan;
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
  }else{
    // Fast SVD along the lines of the R function fast.svd() in package corpcor
    mat u;
    colvec d;
    mat v;
    if(then >= thep){
      mat XTX = myX.t()*myX;
      svd_econ(u,d,v,XTX,"right");
      d = sqrt(d);
      u = myX*v*diagmat(1/d);
    }
    else{
      mat XXT = myX*myX.t();
      svd_econ(u,d,v,XXT,"left");
      d = sqrt(d);
      v = myX.t()*u*diagmat(1/d);
    }
    mat myF = u*diagmat(d);
    mat FTF = myF.t() * myF;
    
    // Initialization
    colvec vec0(thep);
    vec0.zeros();
    vec0.subvec(0,d.n_elem-1) = d % d;
    double tracePostSigma, tempval, valLogDet;
    mat tempPostOmega, postOmega;
    colvec postTheta;
    
    // Algo
    while(mybool){
      
      // Expectations
      expTau = aRandStar/ratesRand(ct-1);
      expSig = cSigmaStar/ratesSig(ct-1);
      
      // Calculate postOmega
      tempPostOmega = diagmat(1/(d%d + expTau*ones(d.n_elem)));
      postOmega = (1/expSig)*tempPostOmega;
      
      // Calculate the trace of postSigma
      tracePostSigma = (1/expSig)*sum(1/(vec0+expTau*ones(vec0.n_elem)));
      
      // Calculate posterior mean of theta
      postTheta = expSig*postOmega*trans(myF)*myy;
      
      // Update posterior rate for sigma^2
      tempval = as_scalar(trans(postTheta)*(postTheta)) + tracePostSigma;
      ratesSig(ct) = dSigma + 0.5*(as_scalar(trans(myy-myF*postTheta)*(myy-myF*postTheta))+trace(FTF*postOmega)) + 0.5*expTau*tempval;
      
      // Update posterior rates for tau^2
      ratesRand(ct) = bRand + 0.5*expSig*tempval;
      
      // Lower bound marginal likelihood
      Lrand = aRand*log(bRand)-aRandStar*log(ratesRand(ct))+lgamma(aRandStar)-lgamma(aRand);
      Lsig = cSigma*log(dSigma)-cSigmaStar*log(ratesSig(ct))+lgamma(cSigmaStar)-lgamma(cSigma);
      theplus = 0.5*(expSig)*(expTau)*(as_scalar(trans(postTheta)*postTheta) + tracePostSigma);
      valLogDet = -thep*log(expSig) - sum(log(vec0+expTau*ones(vec0.n_elem)));
      try{
        L(ct) = 0.5*thep - 0.5*then*log(2*datum::pi) + 0.5*valLogDet + Lrand + Lsig + theplus;
      }
      catch(...){
        L(ct) = datum::nan;
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
  }
  //Rcpp::Rcout << "nb iter = " << ct << std::endl;
  return L(ct-1);
}

// [[Rcpp::export]]
arma::colvec HiddenEdgeBFprime(Rcpp::NumericVector idx, Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX){
  
  //Rcpp::Rcout << "edge = " << idx << std::endl;
  
  mat mymat(themat.begin(),themat.nrow(),themat.ncol(),false);

  // Kappa value
  double kappa = themat(idx(0)-1,idx(1)-1);

  // Indicator for variables to be included in regression 1 (under H0)
  colvec col1 = mymat.col(idx(0)-1);
  colvec vars1 = col1.elem(find(linspace(1,themat.nrow(), themat.nrow())!=(idx(0))));
  colvec varSel1(themat.nrow()-1);
  varSel1.zeros();
  if(any(vars1>kappa)){
    varSel1.elem(find(vars1>kappa)) += 1;
  }

  // Indicator for variables to be included in regression 2 (under H0)
  colvec col2 = mymat.col(idx(1)-1);
  colvec vars2 = col2.elem(find(linspace(1,themat.nrow(), themat.nrow())!=(idx(1))));
  colvec varSel2(themat.nrow()-1);
  varSel2.zeros();
  if(any(vars2>kappa)){
    varSel2.elem(find(vars2>kappa)) += 1;
  }
  
  // Compute variational lower bounds for the models under H0
  double logML01 = HiddenVarRidgei(idx(0), tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel1);
  double logML02 = HiddenVarRidgei(idx(1), tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel2);
  
  // Indicator for variables to be included in regression 1 (under H1)
  varSel1.zeros();
  varSel1.elem(find(vars1>=kappa)) += 1;
  
  // Indicator for variables to be included in regression 2 (under H1)
  varSel2.zeros();
  varSel2.elem(find(vars2>=kappa)) += 1;

  // Compute variational lower bounds for the models under H1
  double logML11 = HiddenVarRidgei(idx(0), tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel1);
  double logML12 = HiddenVarRidgei(idx(1), tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel2);
  
  // logarithm of Bayes factors
  colvec logBFs(2);
  logBFs(0) = logML11-logML01;
  logBFs(1) = logML12-logML02;

  return logBFs;
}

// [[Rcpp::export]]
double HiddenEstimatep0(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX){
  
  mat mymat(themat.begin(),themat.nrow(),themat.ncol(),false);
  mat mymatL = trimatl(mymat);
  colvec uniquevals = unique(mymatL);
  colvec allvals = sort(uniquevals.elem(find(uniquevals!=0)), "descend");
  
  // Algo
  mat tempGraph = zeros(themat.nrow(),themat.ncol());
  colvec allML(themat.nrow());
  for(int j=1; j<=tX.nrow(); j++){
    allML(j-1) = HiddenVarRidgei(j, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, zeros(tX.nrow()-1));
  }
  double logML01, logML02, logML11, logML12;
  int maxedges = 0.5*themat.nrow()*(themat.nrow()-1);
  mat logBFs(maxedges,2);
  logBFs.zeros();
  colvec p0(maxedges);
  p0.zeros();
  uvec idx1, idx2;
  int cpt = 0;
  bool mybool = true;
  while(mybool){
    
    // Consider new edge
    uvec newedge = find(mymat==allvals(cpt))/themat.nrow();
    tempGraph(newedge(1), newedge(0)) = 1;
    tempGraph(newedge(0), newedge(1)) = 1;
    logML01 = allML(newedge(1));
    logML02 = allML(newedge(0));
    
    // Obtain ML for the two regression equations that include to this edge
    colvec varSel1 = tempGraph.col(newedge(1));
    colvec varSel2 = tempGraph.col(newedge(0));
    logML11 = HiddenVarRidgei(newedge(1)+1, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel1.elem(find(linspace(1,tX.nrow(), tX.nrow())!=(newedge(1)+1))));
    logML12 = HiddenVarRidgei(newedge(0)+1, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel2.elem(find(linspace(1,tX.nrow(), tX.nrow())!=(newedge(0)+1))));
    
    allML(newedge(1)) = logML11;
    allML(newedge(0)) = logML12;
    
    // logarithm of Bayes factors
    logBFs(cpt, 0) = logML11-logML01;
    logBFs(cpt, 1) = logML12-logML02;
    
    // proportion of null hypothesis
    idx1 = find(logBFs.col(0)>0);
    idx2 = find(logBFs.col(1)>0);
    p0(cpt) = 1-(((double)idx1.n_elem+(double)idx2.n_elem)/((double)tX.nrow()*((double)tX.nrow()-1)));

    //Rcpp::Rcout << "cpt = " << cpt << " - " << p0(cpt) << " - " << std::endl;
    
    // Convergence
    if(cpt==(maxedges-1)){
      mybool = false;
    }else{
        cpt++;
    }
  }
  
  // OUTPUT
  //uvec idx1 = find(logBFs.col(0)>0);
  //uvec idx2 = find(logBFs.col(1)>0);
  //double p0 = 1-(((double)idx1.n_elem+(double)idx2.n_elem)/((double)tX.nrow()*((double)tX.nrow()-1)));
  
  //return Rcpp::List::create(Rcpp::Named("p0") = p0(cpt), Rcpp::Named("logBFs") = logBFs.rows(0,cpt));
  return p0(cpt);
}

// [[Rcpp::export]]
arma::mat HiddenEdgeSelection(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX, double p0, double lfdrcut, int maxNbEdges){
  
  mat mymat(themat.begin(),themat.nrow(),themat.ncol(),false);
  mat mymatL = trimatl(mymat);
  colvec uniquevals = unique(mymatL);
  colvec allvals = sort(uniquevals.elem(find(uniquevals!=0)), "descend");
  int maxedges = allvals.n_elem;

  // Algo
  mat tempGraph = zeros(themat.nrow(),themat.ncol());
  mat myGraph = zeros(themat.nrow(),themat.ncol());
  colvec allML(themat.nrow());
  for(int j=1; j<=tX.nrow(); j++){
    allML(j-1) = HiddenVarRidgei(j, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, zeros(tX.nrow()-1));
  }
  double logML01, logML02, logML11, logML12;
  int cpt = 0;
  int ctstop = 0;
  int nbSel = 0;
  bool mybool = true;
  double maxBF, minlfdr;
  colvec logBFs(2);
  logBFs.zeros();
  while(mybool){

    // Consider new edge
    uvec newedge = find(mymat==allvals(cpt))/themat.nrow();
    tempGraph(newedge(1), newedge(0)) = 1;
    tempGraph(newedge(0), newedge(1)) = 1;
    logML01 = allML(newedge(1));
    logML02 = allML(newedge(0));
    
    // Obtain ML for the two regression equations that include to this edge
    colvec varSel1 = tempGraph.col(newedge(1));
    colvec varSel2 = tempGraph.col(newedge(0));
    logML11 = HiddenVarRidgei(newedge(1)+1, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel1.elem(find(linspace(1,tX.nrow(), tX.nrow())!=(newedge(1)+1))));
    logML12 = HiddenVarRidgei(newedge(0)+1, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel2.elem(find(linspace(1,tX.nrow(), tX.nrow())!=(newedge(0)+1))));

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
      nbSel++;
      ctstop = 0;
    }else{
      tempGraph(newedge(1), newedge(0)) = 0;
      tempGraph(newedge(0), newedge(1)) = 0;
      ctstop++;
    }
    //Rcpp::Rcout << "cpt = " << cpt << "    nbSel = " << nbSel << std::endl;
    
    // Convergence
    if(cpt==(maxedges-1)){
      mybool = false;
    }else{
      // Stop if the maximum number of edges to select is reached
      if(maxNbEdges==nbSel){
        mybool = false;
      }else{
        cpt++;
      }
    }
  }

  return myGraph;
}

double myf(int myk, mat mm1, mat mm2){
  rowvec myrowvec = mm1.row(myk);
  colvec mycolvec = mm2.col(myk);
  double out = as_scalar(myrowvec*mycolvec);
  return out;
}

Rcpp::List HiddenVarRidgeiOneIter(int ii,  Rcpp::List SVDs, arma::mat tX, double aRand, double bRand, arma::colvec bRandStarInit, arma::colvec dSigmaStarInit){
  
  double cSigma = 0.001;
  double dSigma = 0.001;
  
  // Data
  Rcpp::List tplist = Rcpp::as<Rcpp::List>(SVDs[ii-1]);
  mat myu = Rcpp::as<mat>(tplist["u"]);
  colvec myd = Rcpp::as<colvec>(tplist["d"]);
  mat myv = Rcpp::as<mat>(tplist["v"]);
  mat FTF = diagmat(myd%myd);
  mat myF = myu*diagmat(myd);
  int then = myu.n_rows;
  int thep = myv.n_rows;
  colvec myy = trans(tX.row(ii-1));
  colvec vec0(thep);
  vec0.zeros();
  vec0.subvec(0,myd.n_elem-1) = myd%myd;
  
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

  // Calculate postOmega
  mat tempPostOmega = diagmat(1/(myd%myd + expTau*ones(myd.n_elem)));
  mat postOmega = (1/expSig)*tempPostOmega;

  // Calculate the trace of postSigma
  double tracePostSigma = (1/expSig)*sum(1/(vec0+expTau*ones(vec0.n_elem)));

  // Calculate posterior mean of theta
  colvec postTheta = expSig*postOmega*trans(myF)*myy;

  // Update posterior rate for sigma^2
  double tempval = as_scalar(trans(postTheta)*(postTheta)) + tracePostSigma;
  double dSigmaStar = dSigma + 0.5*(as_scalar(trans(myy-myF*postTheta)*(myy-myF*postTheta))+trace(FTF*postOmega)) + 0.5*expTau*tempval;

  // Update posterior rates for tau^2
  double bRandStar = bRand + 0.5*expSig*tempval;
  
  // Lower bound marginal likelihood
  double Lrand = aRand*log(bRand)-aRandStar*log(bRandStar)+lgamma(aRandStar)-lgamma(aRand);
  double Lsig = cSigma*log(dSigma)-cSigmaStar*log(dSigmaStar)+lgamma(cSigmaStar)-lgamma(cSigma);
  double theplus = 0.5*(expSig)*(expTau)*(as_scalar(trans(postTheta)*postTheta) + tracePostSigma);
  double valLogDet = -thep*log(expSig) - sum(log(vec0+expTau*ones(vec0.n_elem)));
  double L;
  try{
    L = 0.5*thep - 0.5*then*log(2*datum::pi) + 0.5*valLogDet + Lrand + Lsig + theplus;
  }
  catch(...){
    L = datum::nan;
  }

  // Output
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

  return Rcpp::List::create(Rcpp::Named("L") = L, Rcpp::Named("priorRand") = priorRand, Rcpp::Named("priorSig") = priorSig, Rcpp::Named("postRand") = postRand, Rcpp::Named("postSig") = postSig);
}

// [[Rcpp::export]]
arma::colvec HiddenVarRidgeiGetKappa(int ii,  Rcpp::List SVDs, arma::mat tX, double aRand, double bRand, arma::colvec bRandStarInit, arma::colvec dSigmaStarInit){
  
  double cSigma = 0.001;
  double dSigma = 0.001;
  
  // Data
  Rcpp::List tplist = Rcpp::as<Rcpp::List>(SVDs[ii-1]);
  mat myu = Rcpp::as<mat>(tplist["u"]);
  colvec myd = Rcpp::as<colvec>(tplist["d"]);
  mat myv = Rcpp::as<mat>(tplist["v"]);
  mat FTF = diagmat(myd%myd);
  int then = myu.n_rows;
  int thep = myv.n_rows;
  colvec myy = trans(tX.row(ii-1));
  mat myX = trans(tX.rows(find(linspace(1,tX.n_rows, tX.n_rows)!=ii)));
  mat myF = myu*diagmat(myd);

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
  
  // Calculate postOmega
  mat tempPostOmega = diagmat(1/(myd%myd + expTau*ones(myd.n_elem)));
  mat postOmega = (1/expSig)*tempPostOmega;
  
  // Calculate posterior mean of theta
  colvec postTheta = expSig*postOmega*trans(myF)*myy;
  
  // Initialize posterior mean, sd of beta and kappa
  colvec postSd(thep);
  postSd.zeros();
  colvec postMean(thep);
  postMean.zeros();
  colvec kappa(thep+1);
  kappa.zeros();
  
  // Calculation of posterior mean, sd and ratio of beta if light=false only
  postMean = myv*postTheta;
  if(then>=thep){
    // via svd
    mat mymm1 = myv*postOmega;
    for(int k=0; k<thep; k++){
      postSd(k) = sqrt(myf(k, mymm1, trans(myv)));
    }
  }else{
    // via the Woodbury identity (more stable)
    mat myprod = myu * tempPostOmega * trans(myu);
    mat mymm1 = trans(myX)*myprod;
    for(int k=0; k<thep; k++){
      postSd(k) = sqrt(myf(k, mymm1, myX));
    }
  }
  kappa.elem(find(linspace(1,kappa.n_elem, kappa.n_elem)!=ii)) = abs(postMean)/postSd;
  
  return kappa;
}

arma::colvec mydigamma(colvec vec){
  colvec out(vec.n_elem);
  for(int k = 0; k<vec.n_elem; k++){
    out(k) = R::digamma(vec(k));
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List HiddenVarAlgo(Rcpp::List SVDs, arma::mat tX, double aRand, double bRand, int maxiter, int globalShrink, double tol, bool verbose){
  
  // Initialization
  mat allmargs(maxiter+1, SVDs.size());
  mat parTau(maxiter+1, 2);
  parTau(0,0) = aRand;
  parTau(0,1) = bRand;
  int ct = 0;
  bool mybool = true;
  colvec allaRandStar(SVDs.size());
  colvec allbRandStar(SVDs.size());
  colvec alldSigmaStar(SVDs.size());
  allaRandStar.zeros();
  allbRandStar.zeros();
  alldSigmaStar.zeros();
  colvec allbRandStarnew(SVDs.size());
  allbRandStarnew.zeros();
  colvec alldSigmaStarnew(SVDs.size());
  alldSigmaStarnew.zeros();
  Rcpp::List tplist;
  colvec tpvals, tpvals2;
  mat tpmat;
  mat matThres(SVDs.size(), SVDs.size());
  matThres.zeros();
  uvec idxs;
  colvec tpmat2(SVDs.size());
  tpmat2.zeros();
  colvec a1(30);
  colvec b1(30);
  
  // Algo
  while(mybool){
    if(verbose){
      Rcpp::Rcout << "iteration " << ct+1 << std::endl; 
    }
    
    allbRandStar = allbRandStarnew;
    allbRandStarnew.zeros();
    alldSigmaStar = alldSigmaStarnew;
    alldSigmaStarnew.zeros();
    
    // Fit all models
    for(int j=0; j<SVDs.size(); j++){
      //Rcpp::Rcout << "j " << j+1 << std::endl;
      tplist = Rcpp::as<Rcpp::List>(HiddenVarRidgeiOneIter(j+1, SVDs, tX, parTau(ct,0), parTau(ct,1), allbRandStar, alldSigmaStar));
      tpvals = Rcpp::as<colvec>(tplist["postRand"]);
      tpvals2 = Rcpp::as<colvec>(tplist["postSig"]);
      allaRandStar(j) = tpvals(0);
      allbRandStarnew(j) = tpvals(1);
      alldSigmaStarnew(j) = tpvals2(1);
      allmargs(ct,j) = Rcpp::as<double>(tplist["L"]);
    }
    
    // Variational Empirical Bayes using fixed-point iteration as in Valpola and Honkela (2006)
    if(globalShrink==1){
      a1.zeros();
      b1.zeros();
      a1(0) = parTau(ct,0);
      b1(0) = parTau(ct,1);
      int cpt = 1;
      bool mybool2 = true;
      double tp = (allaRandStar.n_elem/sum(allaRandStar/allbRandStarnew));
      double tp2 = mean(log(allbRandStarnew) - mydigamma(allaRandStar))-log(tp);
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
      double tp = sum(allaRandStar/allbRandStarnew);
      double tp2 = mean(log(allbRandStarnew)-mydigamma(allaRandStar));
      parTau(ct+1,0) = 0.5*(1/(log(tp)+tp2-log(SVDs.size())));
      parTau(ct+1,1) = parTau(ct+1,0)*SVDs.size()*(1/tp);
    }//end if
    
    // Monitor convergence
    double maxRelDiffML, diffTotalML;
    if(ct==(maxiter-1)){
      mybool = false;
    }else{
      if(ct>2){
        // Check relative increase for each variational lower bound
        diffTotalML = sum(allmargs.row(ct))-sum(allmargs.row(ct-1));
        maxRelDiffML = max(abs((allmargs.row(ct)-allmargs.row(ct-1))/allmargs.row(ct-1)));
        //Rcpp::Rcout << ",  maxRelDiffML = " << maxRelDiffML << std::endl;
        //Rcpp::Rcout << ",  diffTotalML = " << diffTotalML << std::endl;
        if(diffTotalML>=0){
          if(maxRelDiffML<tol || diffTotalML<0){
            mybool = false;
          }else{
            ct++;
          }
        }else{
          mybool = false;
          ct--;
        }
      }else{
        ct++;
      }
    }
  }//end while
  if(verbose){
    Rcpp::Rcout << "DONE" << std::endl;
  }
  
  return Rcpp::List::create(Rcpp::Named("parTau") = parTau.rows(0,ct), Rcpp::Named("allmargs") = allmargs.rows(0,ct), Rcpp::Named("allbRandStar") = allbRandStar, Rcpp::Named("alldSigmaStar") = alldSigmaStar);
}//end varAlgo

// [[Rcpp::export]]
Rcpp::List getSVD(int ii,  Rcpp::NumericMatrix tX){
  
  // Data
  mat tXbis(tX.begin(),tX.nrow(),tX.ncol(),false);
  colvec myy = trans(tXbis.row(ii-1));
  mat myX = trans(tXbis.rows(find(linspace(1,tXbis.n_rows, tXbis.n_rows)!=ii)));
  int then = myX.n_rows;
  int thep = myX.n_cols;
  
  // Fast SVD (along the lines of the R function fast.svd() in package corpcor)
  mat u;
  colvec d;
  mat v;
  if(then >= thep){
    mat XTX = myX.t()*myX;
    svd_econ(u,d,v,XTX,"right");
    d = sqrt(d);
    u = myX*v*diagmat(1/d);
  }
  else{
    mat XXT = myX*myX.t();
    svd_econ(u,d,v,XXT,"left");
    d = sqrt(d);
    v = myX.t()*u*diagmat(1/d);
  }
  
  //mat myF = u*diagmat(d);
  //mat FTF = myF.t() * myF;

  return Rcpp::List::create(Rcpp::Named("u") = u, Rcpp::Named("d") = d, Rcpp::Named("v") = v);
}

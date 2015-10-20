//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double HiddenVarRidgei(int ii,  Rcpp::NumericMatrix tX, double aRand, double bRand, double cSigma, double dSigma, int maxiter, arma::colvec varSel){
  
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
      if(Ldiff<0.0001){
        mybool = FALSE;
      }
    }
    ct++;
  }
  
  return L(ct-1);
}

// [[Rcpp::export]]
double HiddenEstimatep0(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX, int maxedges){
  
  using arma::trans;
  
  arma::mat mymat(themat.begin(),themat.nrow(),themat.ncol(),false);
  arma::mat mymatL = arma::trimatl(mymat);
  arma::colvec uniquevals = arma::unique(mymatL);
  arma::colvec allvals = arma::sort(uniquevals.elem(arma::find(uniquevals!=0)), "descend");
  
  // Algo
  arma::mat tempGraph = arma::zeros(themat.nrow(),themat.ncol());
  arma::colvec allML(themat.nrow());
  for(int j=1; j<=tX.nrow(); j++){
    allML(j-1) = HiddenVarRidgei(j, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, arma::zeros(tX.nrow()-1));
  }
  double logML01, logML02, logML11, logML12;
  arma::mat logBFs(maxedges,2);
  logBFs.zeros();
  arma::colvec p0(maxedges);
  p0.zeros();
  arma::uvec idx1, idx2;
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
    logML11 = HiddenVarRidgei(newedge(1)+1, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel1.elem(arma::find(arma::linspace(1,tX.nrow(), tX.nrow())!=(newedge(1)+1))));
    logML12 = HiddenVarRidgei(newedge(0)+1, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel2.elem(arma::find(arma::linspace(1,tX.nrow(), tX.nrow())!=(newedge(0)+1))));
    allML(newedge(1)) = logML11;
    allML(newedge(0)) = logML12;
    
    // logarithm of Bayes factors
    logBFs(cpt, 0) = logML11-logML01;
    logBFs(cpt, 1) = logML12-logML02;
    
    // proportion of null hypothesis
    idx1 = arma::find(logBFs.col(0)>0);
    idx2 = arma::find(logBFs.col(1)>0);
    p0(cpt) = 1-(((double)idx1.n_elem+(double)idx2.n_elem)/((double)tX.nrow()*((double)tX.nrow()-1)));

    //Rcpp::Rcout << "cpt = " << cpt << " - " << p0(cpt) << std::endl;
    
    // Convergence
    if(cpt==(maxedges-1)){
      mybool = false;
    }else{
      if(cpt>100){
        if((p0(cpt-100)-p0(cpt))<=0.001){
          mybool = false;
        }else{
          cpt++;
        }
      }else{
        cpt++;
      }
    }
  }
  
  // OUTPUT
  //arma::uvec idx1 = arma::find(logBFs.col(0)>0);
  //arma::uvec idx2 = arma::find(logBFs.col(1)>0);
  //double p0 = 1-(((double)idx1.n_elem+(double)idx2.n_elem)/((double)tX.nrow()*((double)tX.nrow()-1)));
  
  //return List::create(Named("p0") = p0, Named("logBFs") = logBFs.rows(0,cpt));
  return p0(cpt);
}

// [[Rcpp::export]]
arma::mat HiddenEdgeSelection(Rcpp::NumericMatrix themat, Rcpp::NumericMatrix tX, double p0, int maxedges, double lfdrcut){
  
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
    allML(j-1) = HiddenVarRidgei(j, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, arma::zeros(tX.nrow()-1));
  }
  double logML01, logML02, logML11, logML12;
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
    logML11 = HiddenVarRidgei(newedge(1)+1, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel1.elem(arma::find(arma::linspace(1,tX.nrow(), tX.nrow())!=(newedge(1)+1))));
    logML12 = HiddenVarRidgei(newedge(0)+1, tX, 0.5, (double) tX.ncol()/2, 0.001, 0.001, 5000, varSel2.elem(arma::find(arma::linspace(1,tX.nrow(), tX.nrow())!=(newedge(0)+1))));
    
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
    if(cpt==(maxedges-1) || ctstop==200){
      mybool = false;
    }else{
      cpt++;
    }
  }
  
  return myGraph;
}

double myf(int myk, arma::mat mm1, arma::mat mm2){
  arma::rowvec myrowvec = mm1.row(myk);
  arma::colvec mycolvec = mm2.col(myk);
  double out = arma::as_scalar(myrowvec*mycolvec);
  return out;
}

// [[Rcpp::export]]
Rcpp::List HiddenVarRidgeiOneIter(int ii,  Rcpp::List SVDs, double aRand, double bRand, arma::colvec bRandStarInit, arma::colvec dSigmaStarInit, bool light){
  
  //Rcpp::Rcout << "ii= " << ii << std::endl;
  
  using arma::trans;
  double cSigma = 0.001;
  double dSigma = 0.001;
  
  // Data
  Rcpp::List tplist = Rcpp::as<Rcpp::List>(SVDs[ii-1]);
  arma::colvec myy = Rcpp::as<arma::colvec>(tplist["myy"]);
  arma::mat myu = Rcpp::as<arma::mat>(tplist["u"]);
  arma::colvec myd = Rcpp::as<arma::colvec>(tplist["d"]);
  arma::mat myv = Rcpp::as<arma::mat>(tplist["v"]);
  arma::mat FTF = Rcpp::as<arma::mat>(tplist["FTF"]);
  arma::mat myX = Rcpp::as<arma::mat>(tplist["myX"]);
  arma::mat myF = Rcpp::as<arma::mat>(tplist["myF"]);
  int then = myX.n_rows;
  int thep = myX.n_cols;
  arma::colvec vec0(thep);
  vec0.zeros();
  vec0.subvec(0,myd.n_elem-1) = myd % myd;
  
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
  arma::mat tempPostOmega = arma::diagmat(1/(myd%myd + expTau*arma::ones(myd.n_elem)));
  arma::mat postOmega = (1/expSig)*tempPostOmega;

  // Calculate the trace of postSigma
  double tracePostSigma = (1/expSig)*sum(1/(vec0+expTau*arma::ones(vec0.n_elem)));

  // Calculate posterior mean of theta
  arma::colvec postTheta = expSig*postOmega*trans(myF)*myy;

  // Update posterior rate for sigma^2
  double tempval = arma::as_scalar(trans(postTheta)*(postTheta)) + tracePostSigma;
  double dSigmaStar = dSigma + 0.5*(arma::as_scalar(trans(myy-myF*postTheta)*(myy-myF*postTheta))+arma::trace(FTF*postOmega)) + 0.5*expTau*tempval;

  // Update posterior rates for tau^2
  double bRandStar = bRand + 0.5*expSig*tempval;
  
  // Lower bound marginal likelihood
  double Lrand = aRand*log(bRand)-aRandStar*log(bRandStar)+lgamma(aRandStar)-lgamma(aRand);
  double Lsig = cSigma*log(dSigma)-cSigmaStar*log(dSigmaStar)+lgamma(cSigmaStar)-lgamma(cSigma);
  double theplus = 0.5*(expSig)*(expTau)*(arma::as_scalar(trans(postTheta)*postTheta) + tracePostSigma);
  double valLogDet = -thep*log(expSig) - sum(log(vec0+expTau*arma::ones(vec0.n_elem)));
  double L;
  try{
    L = 0.5*thep - 0.5*then*log(2*arma::datum::pi) + 0.5*valLogDet + Lrand + Lsig + theplus;
  }
  catch(...){
    L = arma::datum::nan;
  }

  // Initialize posterior mean, sd and ratio of beta
  arma::colvec postSd(thep);
  postSd.zeros();
  arma::colvec postMean(thep);
  postMean.zeros();
  arma::colvec ratio(thep);
  ratio.zeros();

  // Calculation of posterior mean, sd and ratio of beta if light=false only
  if(light ==  false){
    postMean = myv*postTheta;
    if(then>=thep){
      // via svd
      arma::mat mymm1 = myv*postOmega;
      for(int k=0; k<thep; k++){
        postSd(k) = sqrt(myf(k, mymm1, trans(myv)));
      }
    }else{
      // via the Woodbury identity (more stable)
      arma::mat myprod = myu * tempPostOmega * trans(myu);
      arma::mat mymm1 = trans(myX)*myprod;
      for(int k=0; k<thep; k++){
        postSd(k) = sqrt(myf(k, mymm1, myX));
      }
    }
    ratio = arma::abs(postMean)/postSd;
  }
  
  //if(ii==1){
    //      Rcpp::Rcout << "expTau = " << expTau << std::endl;
    //      Rcpp::Rcout << "bRand = " << bRand << std::endl;
    //      Rcpp::Rcout << "expSig = " << expSig << std::endl;
    //      Rcpp::Rcout << "arma::as_scalar(trans(postMean)*(postMean)) = " << arma::as_scalar(trans(postMean)*(postMean)) << std::endl;
    //Rcpp::Rcout << "myprod.n_rows = " << myprod.n_rows << std::endl;
    //Rcpp::Rcout << "myprod.n_cols = " << myprod.n_cols << std::endl;
    //Rcpp::Rcout << "sum(postSd) = " << sum(postSd) << std::endl;
    //Rcpp::Rcout << "sum(postSd2) = " << sum(postSd) << std::endl;
    //Rcpp::Rcout << "equal ? = " << arma::all(postSigma==postSigma0) << std::endl;
  //}
  
  // Output
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

  return Rcpp::List::create(Rcpp::Named("postBeta") = postBeta, Rcpp::Named("L") = L, Rcpp::Named("priorRand") = priorRand, Rcpp::Named("priorSig") = priorSig, Rcpp::Named("postRand") = postRand, Rcpp::Named("postSig") = postSig);
}

arma::colvec mydigamma(arma::colvec vec){
  arma::colvec out(vec.n_elem);
  for(int k = 0; k<vec.n_elem; k++){
    out(k) = R::digamma(vec(k));
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List HiddenVarAlgo(Rcpp::List SVDs, double aRand, double bRand, int maxiter, int globalShrink, double tol){
  
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
    
    allbRandStar = allbRandStarnew;
    allbRandStarnew.zeros();
    alldSigmaStar = alldSigmaStarnew;
    alldSigmaStarnew.zeros();
    
    // Fit all models
    for(int j=0; j<SVDs.size(); j++){
      //Rcpp::Rcout << "j " << j+1 << std::endl;
      tplist = Rcpp::as<Rcpp::List>(HiddenVarRidgeiOneIter(j+1, SVDs, parTau(ct,0), parTau(ct,1), allbRandStar, alldSigmaStar, true));
      tpvals = Rcpp::as<arma::colvec>(tplist["postRand"]);
      tpvals2 = Rcpp::as<arma::colvec>(tplist["postSig"]);
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
  Rcpp::Rcout << "DONE" << std::endl;
  
  // Calculate posterior statistics
  Rcpp::Rcout << "STEP 2: Calculate summary statistics from posteriors... ";
  for(int j=0; j<SVDs.size(); j++){
    tplist = Rcpp::as<Rcpp::List>(HiddenVarRidgeiOneIter(j+1, SVDs, parTau(ct,0), parTau(ct,1), allbRandStar, alldSigmaStar, false));
    tpmat = Rcpp::as<arma::mat>(tplist["postBeta"]);
    idxs = arma::find(arma::linspace(1,SVDs.size(), SVDs.size())!=(j+1));
    tpmat2.zeros();
    tpmat2.elem(idxs) = tpmat.col(2);
    matThres.col(j) = tpmat2;
  }
  Rcpp::Rcout << "DONE" << std::endl;
  
  // Symmetrize matThres
  matThres = (matThres + matThres.t())/2;
  
  return Rcpp::List::create(Rcpp::Named("matThres") = matThres, Rcpp::Named("parTau") = parTau.rows(0,ct), Rcpp::Named("allmargs") = allmargs.rows(0,ct));
}//end varAlgo

// [[Rcpp::export]]
Rcpp::List HiddenGetSVD(int ii,  Rcpp::NumericMatrix tX){
  
  using arma::trans;
  
  // Data
  arma::mat tXbis(tX.begin(),tX.nrow(),tX.ncol(),false);
  arma::colvec myy = trans(tXbis.row(ii-1));
  arma::mat myX = trans(tXbis.rows(arma::find(arma::linspace(1,tXbis.n_rows, tXbis.n_rows)!=ii)));
  int then = myX.n_rows;
  int thep = myX.n_cols;
  
  // Fast SVD along the lines of the R package corpcor
  arma::mat u;
  arma::colvec d;
  arma::mat v;
  if(then >= thep){
    arma::mat XTX = myX.t()*myX;
    arma::svd_econ(u,d,v,XTX,"right");
    d = sqrt(d);
    u = myX*v*arma::diagmat(1/d);
  }
  else{
    arma::mat XXT = myX*myX.t();
    arma::svd_econ(u,d,v,XXT,"left");
    d = sqrt(d);
    v = myX.t()*u*arma::diagmat(1/d);
  }
  
  arma::mat myF = u*arma::diagmat(d);
  arma::mat FTF = myF.t() * myF;

  return Rcpp::List::create(Rcpp::Named("u") = u, Rcpp::Named("d") = d, Rcpp::Named("v") = v, Rcpp::Named("myF") = myF, Rcpp::Named("FTF") = FTF, Rcpp::Named("myy") = myy, Rcpp::Named("myX") = myX);
}

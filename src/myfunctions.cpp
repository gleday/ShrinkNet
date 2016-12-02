//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rcpp.h>

using namespace arma;

double VarRidgeiP(int ii, arma::mat tX, arma::mat P, double aRand, double bRand, double cSigma, double dSigma, int maxiter, arma::colvec varSel, bool intercept){
  
  //Rcpp::Rcout << "ii= " << ii << std::endl;
  
  // Indexes
  uvec uvecii = find(linspace(0,tX.n_rows-1, tX.n_rows)==ii);
  uvec uvecdiffii = find(linspace(0,tX.n_rows-1, tX.n_rows)!=ii);
  
  // Data
  colvec myy = trans(tX.rows(uvecii));
  mat tXbis = tX.rows(uvecdiffii);
  int then = tX.n_cols;

  // Initialisation
  double cSigmaStar, myL;
  
  if((sum(varSel)==0) && (!intercept)){
    cSigmaStar = cSigma + 0.5*then;
    double G1 = lgamma(cSigmaStar) -lgamma(cSigma) - (then/2)*(log(2*datum::pi)+log(dSigma));
    double G2 = cSigmaStar*log((1/(2*dSigma))*sum(square(tX.row(ii)))+1);
    myL = G1-G2;
  }else{
    // Variable subset
    uvec idxs = find(varSel == 1);
    mat myX;
    if(intercept){
      myX = ones(myy.n_elem);
      if(sum(varSel)>0){
        myX.insert_cols(myX.n_cols, trans(tXbis.rows(idxs)));
      }
    }else{
      myX = trans(tXbis.rows(idxs));
    }
    //Rcpp::Rcout << "myX.row(0) = " << myX.row(0) << std::endl;
    
    int thep = myX.n_cols;
    int ct = 1;
    colvec L = zeros(maxiter+1);
    colvec prP =  nonzeros(P);
    colvec prvalsP = sort(unique(prP));
    mat XTX;
    colvec ratesSig(maxiter+1);
    ratesSig(0) = dSigma;
    bool mybool = true;
    double expTau, expSig, Lrand, Lsig, theplus;
    float Ldiff;
    mat postSigma;
    colvec postMean;
    
    if(prvalsP.n_elem==1){
      
      double aRandStar = aRand+0.5*thep;
      double cSigmaStar = cSigma+0.5*then+0.5*thep;
      colvec ratesRand(maxiter+1);
      ratesRand(0) = bRand;
      
      if(thep<4){

        XTX = myX.t() * myX;
        
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
          XTX = myX.t()*myX;
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
        //Rcpp::Rcout << "nb iter = " << ct << std::endl;
      }
    }else{
      // Resize aRand and bRand according to selected variables
      colvec pr0 = P(uvecdiffii, uvecii);
      colvec pr = pr0.elem(idxs);
      colvec prvals = sort(unique(pr));
      int K = prvals.n_elem;
      colvec ss = zeros(K);
      for(int k=0; k<K; k++){
        ss(k) = sum(pr==prvals(k));
      }// end for
      
      // Update shape parameters
      colvec aRandStar = aRand + 0.5*ss;
      cSigmaStar = cSigma + 0.5*then + 0.5*thep;
      
      // Initialisation Algo
      mat ratesRand(maxiter+1, K);
      ratesRand.zeros();
      ratesRand.row(0) += bRand;
      double valLogDet, sign;
      colvec u = zeros(thep);
      for(int k=0; k<K; k++){
        u.elem(find(pr==prvals(k))) += aRandStar(k);
      }//end for
      colvec mybRandStarInit = zeros(thep);
      mat D;
      colvec myresid, tempval;
      XTX = myX.t() * myX;
      
      // Algo
      while(mybool){
        
        //Rcpp::Rcout << "ct = " << ct << std::endl;
        
        // Intermediate calculus
        mybRandStarInit.zeros();
        for(int k=0; k<K; k++){
          mybRandStarInit.elem(find(pr==prvals(k))) += ratesRand(ct-1,k);
        }//end for
        
        // Expectations
        D = diagmat(u/mybRandStarInit);
        expSig = cSigmaStar/ratesSig(ct-1);
        
        //Update sigma
        if(then>=thep){
          postSigma = inv(expSig*(XTX + D));
        }else{
          mat Dinv = diagmat(1/D.diag());
          postSigma = (1/expSig)*(Dinv - Dinv*myX.t()*inv(eye(then,then)+myX*Dinv*myX.t())*myX*Dinv);
        }//end if
        colvec v = postSigma.diag();
        
        // Update posterior expectation of beta
        postMean = expSig*postSigma*trans(myX)*myy;
        
        // Update posterior rate for sigma^2
        myresid = myy-myX*postMean;
        tempval = v+square(postMean);
        ratesSig(ct) = dSigma + 0.5*(as_scalar(trans(myresid)*(myresid))+trace(XTX*postSigma)) + 0.5*sum(D.diag()%tempval);
        
        // Update posterior rates for tau^2
        for(int k=0; k<K; k++){
          ratesRand(ct, k) = bRand + 0.5*expSig*sum(tempval.elem(find(pr==prvals(k))));
        }//end for
        
        // Variational lower bound on log-marginal likelihood
        Lrand = 0;
        for(int i=0; i<K; i++){
          Lrand += aRand*log(bRand)-aRandStar(i)*log(ratesRand(ct,i))+lgamma(aRandStar(i))-lgamma(aRand);
        }//end for
        Lsig = cSigma*log(dSigma)-cSigmaStar*log(ratesSig(ct))+lgamma(cSigmaStar)-lgamma(cSigma);
        theplus = 0.5*(cSigmaStar/ratesSig(ct))*sum((u/mybRandStarInit)%tempval);
        log_det(valLogDet, sign, postSigma);
        try{
          L(ct) = 0.5*thep - 0.5*then*log(2*datum::pi) + 0.5*valLogDet + Lrand + Lsig + theplus;
        }
        catch(...){
          L(ct) = datum::nan;
        }
        
        // Monitor convergence
        if(ct==maxiter){
          mybool = FALSE;
        }else{
          if(ct>2){
            Ldiff = std::abs(L(ct)-L(ct-1));
            if(Ldiff<0.001){
              mybool = FALSE;
            }//end if
          }//end if
        }//end if
        ct++;
      }//end while
    }// end if
    myL = L(ct-1);
  }
  return myL;
}//end VarRidgeiP

// [[Rcpp::export]]
arma::colvec HiddenEdgeBFprimeP(Rcpp::NumericVector idx, arma::mat themat, arma::mat tX, arma::mat P, double cSigma, double dSigma, bool intercept){
  
  //Rcpp::Rcout << "edge = " << idx << std::endl;
  
  //mat mymat(themat.begin(),themat.nrow(),themat.ncol(),false);
  
  // Kappa value
  double kappa = themat(idx(0),idx(1));
  
  // Indicator for variables to be included in regression 1 (under H0)
  colvec col1 = themat.col(idx(0));
  colvec vars1 = col1.elem(find(linspace(0,themat.n_rows-1, themat.n_rows)!=(idx(0))));
  colvec varSel1(themat.n_rows-1);
  varSel1.zeros();
  if(any(vars1>kappa)){
    varSel1.elem(find(vars1>kappa)) += 1;
  }
  
  // Indicator for variables to be included in regression 2 (under H0)
  colvec col2 = themat.col(idx(1));
  colvec vars2 = col2.elem(find(linspace(0,themat.n_rows-1, themat.n_rows)!=(idx(1))));
  colvec varSel2(themat.n_rows);
  varSel2.zeros();
  if(any(vars2>kappa)){
    varSel2.elem(find(vars2>kappa)) += 1;
  }
  
  // Compute variational lower bounds for the models under H0
  double logML01 = VarRidgeiP(idx(0), tX, P, 0.5, (double) tX.n_cols/2, cSigma, dSigma, 5000, varSel1, intercept);
  double logML02 = VarRidgeiP(idx(1), tX, P, 0.5, (double) tX.n_cols/2, cSigma, dSigma, 5000, varSel2, intercept);
  
  // Indicator for variables to be included in regression 1 (under H1)
  varSel1.zeros();
  varSel1.elem(find(vars1>=kappa)) += 1;
  
  // Indicator for variables to be included in regression 2 (under H1)
  varSel2.zeros();
  varSel2.elem(find(vars2>=kappa)) += 1;
  
  // Compute variational lower bounds for the models under H1
  double logML11 = VarRidgeiP(idx(0), tX, P, 0.5, (double) tX.n_cols/2, cSigma, dSigma, 5000, varSel1, intercept);
  double logML12 = VarRidgeiP(idx(1), tX, P, 0.5, (double) tX.n_cols/2, cSigma, dSigma, 5000, varSel2, intercept);
  
  // logarithm of Bayes factors
  colvec logBFs(2);
  logBFs(0) = logML11-logML01;
  logBFs(1) = logML12-logML02;
  
  return logBFs;
}// end HiddenEdgeBFprimeP

// [[Rcpp::export]]
Rcpp::List HiddenEstimatep0P(arma::mat themat, arma::mat tX, arma::mat P, double cSigma, double dSigma, bool intercept){
  
  //mat mymat(themat.begin(),themat.nrow(),themat.ncol(),false);
  mat mymatL = trimatl(themat);
  colvec uniquevals = unique(mymatL);
  colvec allvals = sort(uniquevals.elem(find(uniquevals!=0)), "descend");
  
  // Algo
  mat tempGraph = zeros(themat.n_rows,themat.n_cols);
  colvec allML(themat.n_rows);
  for(int j=0; j<tX.n_rows; j++){
    allML(j) = VarRidgeiP(j, tX, P, 0.5, (double) tX.n_cols/2, cSigma, dSigma, 5000, zeros(tX.n_rows-1), intercept);
  }
  //Rcpp::Rcout << "mean(allML) = " << mean(allML) << std::endl;
  double logML01, logML02, logML11, logML12;
  int maxedges = 0.5*themat.n_rows*(themat.n_rows-1);
  mat logBFs(maxedges,3);
  logBFs.zeros();
  colvec p0(maxedges);
  p0.zeros();
  uvec idx1, idx2;
  int cpt = 0;
  bool mybool = true;
  while(mybool){
    // Consider new edge
    uvec newedge = find(themat==allvals(cpt))/themat.n_rows;
    tempGraph(newedge(1), newedge(0)) = 1;
    tempGraph(newedge(0), newedge(1)) = 1;
    logML01 = allML(newedge(1));
    logML02 = allML(newedge(0));
    
    // Obtain ML for the two regression equations that include to this edge
    colvec varSel1 = tempGraph.col(newedge(1));
    colvec varSel2 = tempGraph.col(newedge(0));
    logML11 = VarRidgeiP(newedge(1), tX, P, 0.5, (double) tX.n_cols/2, cSigma, dSigma, 5000, varSel1.elem(find(linspace(0,tX.n_rows-1, tX.n_rows)!=(newedge(1)))), intercept);
    logML12 = VarRidgeiP(newedge(0), tX, P, 0.5, (double) tX.n_cols/2, cSigma, dSigma, 5000, varSel2.elem(find(linspace(0,tX.n_rows-1, tX.n_rows)!=(newedge(0)))), intercept);
    
    allML(newedge(1)) = logML11;
    allML(newedge(0)) = logML12;
    
    // logarithm of Bayes factors
    logBFs(cpt, 0) = allvals(cpt);
    logBFs(cpt, 1) = logML11-logML01;
    logBFs(cpt, 2) = logML12-logML02;
    
    // proportion of null hypothesis
    idx1 = find(logBFs.col(1)>0);
    idx2 = find(logBFs.col(2)>0);
    p0(cpt) = 1-(((double)idx1.n_elem+(double)idx2.n_elem)/((double)tX.n_rows*((double)tX.n_rows-1)));
    
    //Rcpp::Rcout << "cpt = " << cpt << " - " << p0(cpt) << " - " << std::endl;
    //Rcpp::Rcout << "cpt = " << cpt << " - " << "mean(allML)=" << mean(allML) << std::endl;
    
    // Convergence
    if(cpt==(maxedges-1)){
      mybool = false;
    }else{
      cpt++;
    }
  }
  
  // OUTPUT
  return Rcpp::List::create(Rcpp::Named("p0") = p0(cpt), Rcpp::Named("logBFs") = logBFs.rows(0,cpt));
}// end HiddenEstimatep0P

// [[Rcpp::export]]
Rcpp::List HiddenEdgeSelectionP(arma::mat themat, arma::mat tX, arma::mat P, arma::mat weights, double p0, double lfdrcut, int maxNbEdges, double cSigma, double dSigma, bool intercept){
  
  mat mymatL = trimatl(themat);
  colvec uniquevals = unique(mymatL);
  colvec allvals = sort(uniquevals.elem(find(uniquevals!=0)), "descend");
  int maxedges = allvals.n_elem;
  // Algo
  mat tempGraph = zeros(themat.n_rows,themat.n_cols);
  mat myGraph = zeros(themat.n_rows,themat.n_cols);
  mat tempAllLogMaxBF = zeros(themat.n_rows,themat.n_cols);
  colvec allML(themat.n_rows);
  for(int j=0; j<tX.n_rows; j++){
    allML(j) = VarRidgeiP(j, tX, P, 0.5, (double) tX.n_cols/2, cSigma, dSigma, 5000, zeros(tX.n_rows-1), intercept);
  }
  double logML01, logML02, logML11, logML12;
  int cpt = 0;
  int ctstop = 0;
  int nbSel = 0;
  bool mybool = true;
  double maxBF, minlfdr, edgeWeight;
  colvec logBFs(2);
  logBFs.zeros();
  while(mybool){
    
    // Consider new edge
    uvec newedge = find(themat==allvals(cpt))/themat.n_rows;
    tempGraph(newedge(1), newedge(0)) = 1;
    tempGraph(newedge(0), newedge(1)) = 1;
    logML01 = allML(newedge(1));
    logML02 = allML(newedge(0));
    
    // Obtain ML for the two regression equations that include to this edge
    colvec varSel1 = tempGraph.col(newedge(1));
    colvec varSel2 = tempGraph.col(newedge(0));
    logML11 = VarRidgeiP(newedge(1), tX, P, 0.5, (double) tX.n_cols/2, cSigma, dSigma, 5000, varSel1.elem(find(linspace(0,tX.n_rows-1, tX.n_rows)!=(newedge(1)))), intercept);
    logML12 = VarRidgeiP(newedge(0), tX, P, 0.5, (double) tX.n_cols/2, cSigma, dSigma, 5000, varSel2.elem(find(linspace(0,tX.n_rows-1, tX.n_rows)!=(newedge(0)))), intercept);

    // logarithm of Bayes factors
    logBFs(0) = logML11-logML01;
    logBFs(1) = logML12-logML02;
    maxBF = exp(max(logBFs));
    edgeWeight = weights(newedge(1), newedge(0));
    if((edgeWeight>=0) && (edgeWeight<=1)){
      minlfdr = edgeWeight/(maxBF*(1-edgeWeight)+edgeWeight);
    }else{
      minlfdr = p0/(maxBF*(1-p0)+p0);
    }
    
    // Thresholding - edge selection
    if(minlfdr<=lfdrcut){
      myGraph(newedge(1), newedge(0)) = 1;
      myGraph(newedge(0), newedge(1)) = 1;
      tempAllLogMaxBF(newedge(1), newedge(0)) = max(logBFs);
      tempAllLogMaxBF(newedge(0), newedge(1)) = max(logBFs);
      allML(newedge(1)) = logML11;
      allML(newedge(0)) = logML12;
      nbSel++;
      ctstop = 0;
    }else{
      tempGraph(newedge(1), newedge(0)) = 0;
      tempGraph(newedge(0), newedge(1)) = 0;
      ctstop++;
    }
    //Rcpp::Rcout << "cpt = " << cpt << "    nbSel = " << nbSel << "    kappa = " << allvals(cpt) << std::endl;

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
  
  //return myGraph;
  return Rcpp::List::create(Rcpp::Named("myGraph") = myGraph, Rcpp::Named("logMaxBFs") = tempAllLogMaxBF, Rcpp::Named("sumLogML") = sum(allML));
}// end HiddenEdgeSelectionP

double myf(int myk, mat mm1, mat mm2){
  rowvec myrowvec = mm1.row(myk);
  colvec mycolvec = mm2.col(myk);
  double out = as_scalar(myrowvec*mycolvec);
  return out;
}

Rcpp::List HiddenVarRidgeiOneIter(int ii,  Rcpp::List SVDs, arma::mat tX, double aRand, double bRand, arma::colvec bRandStarInit, arma::colvec dSigmaStarInit, double cSigma, double dSigma){
  
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
arma::mat HiddenVarRidgeiGetKappa(int ii,  Rcpp::List SVDs, arma::mat tX, double aRand, double bRand, arma::colvec bRandStarInit, arma::colvec dSigmaStarInit, double cSigma, double dSigma){
  
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
  colvec mypostMean(thep+1);
  mypostMean.zeros();
  mat out(thep+1, 2);
  out.zeros();
  
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
  mypostMean.elem(find(linspace(1,mypostMean.n_elem, mypostMean.n_elem)!=ii)) = postMean;
  out.col(0) = kappa;
  out.col(1) = mypostMean;
  
  return out;
}

// Class
class myVarRidgei{
private:
  double L;
  mat priorRand;
  colvec priorSig;
  mat postRand;
  colvec postSig;
  mat postBeta;
public:
  double get_L();
  mat get_priorRand();
  colvec get_priorSig();
  mat get_postRand();
  colvec get_postSig();
  mat get_postBeta();
  
  void set_L(double l);
  void set_priorRand(mat priorR);
  void set_priorSig(colvec priorS);
  void set_postRand(mat postR);
  void set_postSig(colvec postS);
  void set_postBeta(mat postB);
};

// Getters
double myVarRidgei::get_L(){
  return L;
}
mat myVarRidgei::get_priorRand(){
  return priorRand;
}
colvec myVarRidgei::get_priorSig(){
  return priorSig;
}
mat myVarRidgei::get_postRand(){
  return postRand;
}
colvec myVarRidgei::get_postSig(){
  return postSig;
}
mat myVarRidgei::get_postBeta(){
  return postBeta;
}
//Setters
void myVarRidgei::set_L(double l){
  L = l;
}
void myVarRidgei::set_priorRand(mat priorR){
  priorRand = priorR;
}
void myVarRidgei::set_priorSig(colvec priorS){
  priorSig = priorS;
}
void myVarRidgei::set_postRand(mat postR){
  postRand = postR;
}
void myVarRidgei::set_postSig(colvec postS){
  postSig = postS;
}
void myVarRidgei::set_postBeta(mat postB){
  postBeta = postB;
}

myVarRidgei VarRidgeiOneIterP(int ii, arma::mat tX, arma::mat P, arma::colvec aRand, arma::colvec bRand, arma::colvec bRandStarInit, arma::colvec dSigmaStarInit, double cSigma, double dSigma){
  
  // Data
  uvec idxs = find((linspace(1, tX.n_rows, tX.n_rows)-1)!=ii);
  uvec uvecii = find((linspace(1, tX.n_rows, tX.n_rows)-1)==ii);
  mat myX = trans(tX.rows(idxs));
  int then = myX.n_rows;
  int thep = myX.n_cols;
  colvec myy = trans(tX.row(ii));
  mat XTX = myX.t() * myX;
  
  colvec prP =  nonzeros(P);
  colvec prvalsP = sort(unique(prP));
  colvec pr = P.submat(idxs, uvecii);
  colvec prvals = sort(unique(pr));
  int K = prvals.n_elem;
  colvec ss = zeros(K);
  for(int k=0; k<K; k++){
    ss(k) = sum(pr==prvals(k));
  }
  colvec aRand2, bRand2, bRandStarInit2;
  if(K!=prvalsP.n_elem){
    aRand2 = zeros(K);
    bRand2 = zeros(K);
    bRandStarInit2 = zeros(K);
    for(int k=0; k<K; k++){
      aRand2(find(prvals==prvals(k))) = aRand(find(prvalsP==prvals(k)));
      bRand2(find(prvals==prvals(k))) = bRand(find(prvalsP==prvals(k)));
      bRandStarInit2(find(prvals==prvals(k))) = bRandStarInit(find(prvalsP==prvals(k)));
    }
  }else{
    aRand2 = aRand;
    bRand2 = bRand;
    bRandStarInit2 = bRandStarInit;
  }//end if

  // Update posterior shape parameters
  colvec aRandStar = aRand2 + 0.5*ss;
  double cSigmaStar = cSigma + 0.5*then + 0.5*thep;
  
  // Intermediate calculus
  colvec u = zeros(thep);
  colvec mybRandStarInit = zeros(thep);
  for(int k=0; k<K; k++){
    u.elem(find(pr==prvals(k))) += aRandStar(k);
    if(all(bRandStarInit2==0)){
      mybRandStarInit.elem(find(pr==prvals(k))) += bRand2(k);
    }else{
      mybRandStarInit.elem(find(pr==prvals(k))) += bRandStarInit2(k);
    }//end if
  }//end for
  
  // Expectations
  double expSig;
  mat D;
  if(all(bRandStarInit2==0)){
    D = diagmat(u/mybRandStarInit);
    expSig = cSigmaStar/dSigma;
  }else{
    D = diagmat(u/mybRandStarInit);
    expSig = cSigmaStar/dSigmaStarInit(ii);
  }//end if
  
  //Update sigma
  mat postSigma;
  if(then>=thep){
    postSigma = inv(expSig*(XTX + D));
  }else{
    mat Dinv = diagmat(1/D.diag());
    postSigma = (1/expSig)*(Dinv - Dinv*myX.t()*inv(eye(then,then)+myX*Dinv*myX.t())*myX*Dinv);
  }//end if
  colvec v = postSigma.diag();
  
  // Calculate posterior mean of beta
  colvec postMean = expSig*postSigma*trans(myX)*myy;
  
  // Update posterior rate for sigma^2
  colvec myresid = myy-myX*postMean;
  colvec tempval = v+square(postMean);
  double dSigmaStar = dSigma + 0.5*(as_scalar(trans(myresid)*(myresid))+trace(XTX*postSigma)) + 0.5*sum(D.diag()%tempval);
  
  // Update posterior rates for tau^2
  colvec bRandStar = zeros(K);
  for(int i=0; i<K; i++){
    bRandStar(i) = bRand2(i) + 0.5*expSig*sum(tempval.elem(find(pr==prvals(i))));
  }//end for
  
  // Lower bound marginal likelihood
  double Lrand = 0;
  for(int i=0; i<K; i++){
    Lrand += aRand2(i)*log(bRand2(i))-aRandStar(i)*log(bRandStar(i))+lgamma(aRandStar(i))-lgamma(aRand2(i));
  }//end for
  double Lsig = cSigma*log(dSigma)-cSigmaStar*log(dSigmaStar)+lgamma(cSigmaStar)-lgamma(cSigma);
  double theplus = 0.5*(cSigmaStar/dSigmaStar)*sum((u/mybRandStarInit)%tempval);
  double valLogDet, sign;
  log_det(valLogDet, sign, postSigma);
  double L;
  try{
    L = 0.5*thep - 0.5*then*log(2*datum::pi) + 0.5*valLogDet + Lrand + Lsig + theplus;
  }
  catch(...){
    L = datum::nan;
  }
  
  // priorRand matrix
  mat priorRand(K,2);
  priorRand.zeros();
  priorRand.col(0) += aRand2;
  priorRand.col(1) += bRand2;
  // priorSig vector
  colvec priorSig(2);
  priorSig.zeros();
  priorSig(0) = cSigma;
  priorSig(1) = dSigma;
  // postRand matrix
  mat postRand(prvalsP.n_elem,2);
  postRand.zeros();
  colvec col0 = postRand.col(0);
  colvec col1 = postRand.col(1);
  for(int k=0; k<K; k++){
    col0.elem(find(prvalsP==prvals(k))) += aRandStar(k);
    col1.elem(find(prvalsP==prvals(k))) += bRandStar(k);
  }
  postRand.col(0) += col0;
  postRand.col(1) += col1;
  // postSig vector
  colvec postSig(2);
  postSig.zeros();
  postSig(0) = cSigmaStar;
  postSig(1) = dSigmaStar;
  // postBeta matrix
  mat postBeta(thep,3);
  postBeta.zeros();
  postBeta.col(0) += postMean;
  postBeta.col(1) += sqrt(v);
  postBeta.col(2) += abs(postMean)/sqrt(v);
  
  // Output
  myVarRidgei out;
  out.set_L(L);
  out.set_priorRand(priorRand);
  out.set_priorSig(priorSig);
  out.set_postRand(postRand);
  out.set_postSig(postSig);
  out.set_postBeta(postBeta);
  
  return out;
}// end VarRidgeiOneIterP

arma::colvec mydigamma(colvec vec){
  colvec out(vec.n_elem);
  for(int k = 0; k<vec.n_elem; k++){
    out(k) = R::digamma(vec(k));
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List HiddenVarAlgo(Rcpp::List SVDs, arma::mat tX, double aRand, double bRand, double cSigma, double dSigma, int maxiter, int globalShrink, double tol, bool verbose){
  
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
      tplist = Rcpp::as<Rcpp::List>(HiddenVarRidgeiOneIter(j+1, SVDs, tX, parTau(ct,0), parTau(ct,1), allbRandStar, alldSigmaStar, cSigma, dSigma));
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

arma::colvec fixedPointIterEB(arma::colvec initab, arma::colvec myallaRandStar, arma::colvec myallbRandStar, int mymaxiter, double myeps){
  colvec myallaRandStar2 = nonzeros(myallaRandStar);
  colvec myallbRandStar2 = nonzeros(myallbRandStar);
  colvec a(mymaxiter+1);
  colvec b(mymaxiter+1);
  a.zeros();
  b.zeros();
  a(0) = initab(0);
  b(0) = initab(1);
  int cpt = 1;
  bool mybool = true;
  double tp = (myallaRandStar2.n_elem/sum(myallaRandStar2/myallbRandStar2));
  double tp2 = mean(log(myallbRandStar2) - mydigamma(myallaRandStar2))-log(tp);
  while(mybool){
    a(cpt) = a(cpt-1) + 0.5*(1/(R::digamma(a(cpt-1))-log(a(cpt-1)))) + 0.5*(1/tp2);
    b(cpt) = a(cpt)*tp;
    if((abs(a(cpt)-a(cpt-1))<myeps) && (abs(b(cpt)-b(cpt-1))<myeps)){
      mybool = false;
    }else{
      if(cpt==mymaxiter){
        mybool = false;
      }else{
        cpt++;
      }//end if
    }//end if
  }//end while
  colvec ab(2);
  ab(0) = a(cpt);
  ab(1) = b(cpt);
  
  return ab;
}

arma::colvec approximateEB(arma::colvec myallaRandStar, arma::colvec myallbRandStar){
  colvec myallaRandStar2 = nonzeros(myallaRandStar);
  colvec myallbRandStar2 = nonzeros(myallbRandStar);
  double tp = sum(myallaRandStar2/myallbRandStar2);
  double tp2 = mean(log(myallbRandStar2)-mydigamma(myallaRandStar2));
  double a = 0.5*(1/(log(tp)+tp2-log(myallaRandStar2.n_rows)));
  double b = a*myallaRandStar2.n_rows*(1/tp);
  colvec ab(2);
  ab(0) = a;
  ab(1) = b;
  return ab;
}

// [[Rcpp::export]]
Rcpp::List HiddenVarAlgoP(arma::mat tX, arma::mat P, double aRandInit, double bRandInit, int maxiter, int globalShrink, double tol, bool verbose, double cSigma, double dSigma){
  
  // Initialization
  mat allmargs(maxiter+1, tX.n_rows);
  colvec pr =  nonzeros(P);
  colvec prvals = sort(unique(pr));
  int K = prvals.n_elem;
  mat parTau(maxiter+1, 2*K);
  parTau.zeros();
  uvec aRandIdxs = regspace<uvec>(1, 2, 2*K-1) - 1;
  uvec bRandIdxs = regspace<uvec>(2, 2, 2*K) - 1;
  for(int k=0; k<K; k++){
    parTau(0,aRandIdxs(k)) += aRandInit;
    parTau(0,bRandIdxs(k)) += bRandInit;
  }//end for
  int ct = 0;
  bool mybool = true;
  mat allaRandStar(tX.n_rows, K);
  mat allbRandStar(tX.n_rows, K);
  mat allbRandStarnew(tX.n_rows, K);
  colvec alldSigmaStar(tX.n_rows);
  colvec alldSigmaStarnew(tX.n_rows);
  allaRandStar.zeros();
  allbRandStar.zeros();
  alldSigmaStar.zeros();
  allbRandStarnew.zeros();
  alldSigmaStarnew.zeros();
  myVarRidgei obj;
  mat tpvals;
  colvec tpvals2;
  colvec myab(2);
  myab.zeros();
  mat matThres(tX.n_rows, tX.n_rows);
  matThres.zeros();
  mat matBeta(tX.n_rows, tX.n_rows);
  matBeta.zeros();
  mat tpmat;
  uvec uvecct;
  
  // Algo
  while(mybool){
    if(verbose){
      Rcpp::Rcout << "iteration " << ct+1 << std::endl; 
    }//end if
    
    allbRandStar = allbRandStarnew;
    allbRandStarnew.zeros();
    alldSigmaStar = alldSigmaStarnew;
    alldSigmaStarnew.zeros();
    uvecct = ct;

    // Fit all models
    for(int j=0; j<tX.n_rows; j++){
      obj = VarRidgeiOneIterP(j, tX, P, trans(parTau.submat(uvecct,aRandIdxs)), trans(parTau.submat(uvecct,bRandIdxs)), trans(allbRandStar.row(j)), alldSigmaStar, cSigma, dSigma);
      tpmat = obj.get_postBeta();
      matThres(find(linspace(0,tX.n_rows-1, tX.n_rows)!=j), find(linspace(0,tX.n_rows-1, tX.n_rows)==j)) += tpmat.col(2);
      matBeta(find(linspace(0,tX.n_rows-1, tX.n_rows)!=j), find(linspace(0,tX.n_rows-1, tX.n_rows)==j)) += tpmat.col(0);
      tpvals = obj.get_postRand();
      tpvals2 = obj.get_postSig();
      allaRandStar.row(j) = trans(tpvals.col(0));
      allbRandStarnew.row(j) = trans(tpvals.col(1));
      alldSigmaStarnew(j) = tpvals2(1);
      allmargs(ct,j) = obj.get_L();
      //Rcpp::Rcout << "allmargs(ct,j) " << allmargs(ct,j) << std::endl; 
    }//end for
    //Rcpp::Rcout << "mean log-ML =" << mean(allmargs.row(ct)) << " sum log-ML ="<< sum(allmargs.row(ct)) << std::endl; 
    // Variational Empirical Bayes using fixed-point iteration as in Valpola and Honkela (2006)
    if(globalShrink==1){
      for(int ii=0; ii<K; ii++){
        myab(0) = parTau(ct,0+2*ii);
        myab(1) = parTau(ct,1+2*ii);
        myab = fixedPointIterEB(myab, allaRandStar.col(ii), allbRandStarnew.col(ii), 50, 1e-6);
        parTau(ct+1,0+2*ii) = myab(0);
        parTau(ct+1,1+2*ii) = myab(1);
      }// end for
    }//end if
    
    // Variational Empirical Bayes using approximate analytical solution as in Leday et al (2015)
    if(globalShrink==2){
      for(int ii=0; ii<K; ii++){
        myab = approximateEB(allaRandStar.col(ii), allbRandStarnew.col(ii));
        parTau(ct+1,0+2*ii) = myab(0);
        parTau(ct+1,1+2*ii) = myab(1);
      }// end for
    }//end if
    
    // Monitor convergence
    double maxRelDiffML, diffTotalML;
    if(ct==(maxiter-1)){
      mybool = false;
    }else{
      if(ct>20){
        // Check relative increase for each variational lower bound
        diffTotalML = sum(allmargs.row(ct))-sum(allmargs.row(ct-1));
        maxRelDiffML = max(abs((allmargs.row(ct)-allmargs.row(ct-1))/allmargs.row(ct-1)));
        if(diffTotalML>=0){
          if(maxRelDiffML<tol || diffTotalML<0){
            mybool = false;
          }else{
            ct++;
          }//end if
        }else{
          mybool = false;
          ct--;
        }//end if
      }else{
        ct++;
      }//end if
    }//end if
  }//end while
  if(verbose){
    Rcpp::Rcout << "DONE" << std::endl;
  }
  
  matThres = (matThres + trans(matThres))/2;
  return Rcpp::List::create(Rcpp::Named("parTau") = parTau.rows(0,ct), Rcpp::Named("allmargs") = allmargs.rows(0,ct), Rcpp::Named("matThres") = matThres, Rcpp::Named("matBeta") = matBeta);
}//end varAlgoP


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

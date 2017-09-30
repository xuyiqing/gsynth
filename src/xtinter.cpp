# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

/* ******************* Useful Functions  *********************** */

/* cross product */
arma::mat crossprod (arma::mat x, arma::mat y) {
  return(x.t() * y);
}

/* Three dimensional matrix inverse */
// [[Rcpp::export]]
arma::mat XXinv (arma::cube X) { 
  int p = X.n_slices ;
  arma::mat xx(p, p) ;
  for (int k = 0; k < p; k++) {
    for (int m = 0; m < p; m++) {
      xx(k, m) = trace(crossprod(X.slice(k), X.slice(m))) ;
      if (k < m) {
        xx(m, k) = xx(k, m);
      }
    }
  } 
  return(inv(xx)) ;
}  

arma::mat DM (int m) { /* Demeaning */
  arma::mat A(m, 1, arma::fill::ones)  ; 
  arma::mat MM = arma::eye<arma::mat>(m,m) - A * A.t()/m ;
  return(MM) ;
}


    /* ******************* Subsidiary Functions  *********************** */


/* Obtain OLS panel estimate */
// [[Rcpp::export]]
arma::mat panel_est (arma::cube X, arma::mat Y, arma::mat MF) {
  int p = X.n_slices ;
  arma::mat xx(p, p, arma::fill::zeros) ;
  arma::mat xy(p, 1, arma::fill::zeros);
  if (p==1) {
    xx(0,0) = trace(X.slice(0).t() * MF * X.slice(0)) ;
    xy(0,0) = trace(X.slice(0).t() * MF * Y) ;
  }
  if (p>1) {
    for (int k = 0; k < p; k++) {
      arma::mat MFX1 = MF * X.slice(k) ;
      xy(k, 0) = trace(crossprod(MFX1, Y)) ;
      for (int m = k; m < p; m++) {
        arma::mat MFX2 = MF * X.slice(m) ;
        xx(k, m) = trace(crossprod(MFX1, MFX2)) ;
        if (k < m) {
          xx(m, k) = xx(k, m) ;
        }
      }
    }
  }
  return(xx.i() * xy) ;
}

/* Obtain beta given factors and loadings */
// [[Rcpp::export]]
arma::mat panel_beta (arma::cube X, arma::mat xxinv,
                      arma::mat Y, arma::mat F, arma::mat L) {
  int p = X.n_slices ; 
  arma::mat xy(p, 1, arma::fill::zeros) ;
  for (int k = 0; k < p; k++) {
    xy(k) = trace(crossprod(X.slice(k), (Y - F * L.t()))) ;
  }
  return(xxinv * xy);
}

/* Obtain factors and loading given error */
// [[Rcpp::export]]
List panel_factor (arma::mat E, int r) {
  int T = E.n_rows ;
  int N = E.n_cols ;
  arma::mat factor(T, r) ;
  arma::mat lambda(N, r) ;
  arma::mat VNT(r, r) ;
  arma::mat U ;
  arma::vec s ;
  arma::mat V ; 
  if (T < N) { 
    arma::mat EE = E * E.t() /(N * T) ;
    arma::svd( U, s, V, EE) ;
    factor = U.head_cols(r) * sqrt(T) ;
    lambda = E.t() * factor/T ;
    VNT = diagmat(s.head_rows(r)) ;
  } 
  else {
    arma::mat EE = E.t() * E / (N * T) ;
    svd(U, s, V, EE) ;
    lambda = U.head_cols(r) * sqrt(N) ;
    factor = E * lambda / N ;
    VNT = diagmat(s.head_rows(r)) ;
  }
  List result ;
  result["lambda"] = lambda ;
  result["factor"] = factor ;
  result["VNT"] = VNT ;
  return(result) ;
  
}

/* Main iteration for beta */
// [[Rcpp::export]]
List beta_iter (arma::cube X,
                arma::mat xxinv,
                arma::mat Y,
                int r,
                double tolerate,
                arma::mat beta0) {

  /* beta.new: computed beta under iteration with error precision=tolerate
     factor: estimated factor
     lambda: estimated loadings
     V: the eigenvalues matrix
     e: estimated residuals
     niter: number of interations to achieve convergence */
  
  int p = X.n_slices ;
  double beta_norm = 1.0 ;
  arma::mat beta = beta0 ;
  arma::mat beta_old(p, 1) ;
  arma::mat VNT(r, r) ;

  /* starting value */
  arma::mat U = Y ;
  for (int k = 0; k < p; k++) {
    U = U - X.slice(k) * beta(k) ;
  }
  List pf = panel_factor(U, r)  ;
  arma::mat F = as<arma::mat>(pf["factor"]) ;
  arma::mat L = as<arma::mat>(pf["lambda"]) ;
 
  /* Loop */
  int niter = 0 ;
  while ((beta_norm > tolerate) && (niter < 500)) {
    niter++ ;
    beta = panel_beta(X, xxinv, Y, F, L) ;
    beta_norm = arma::norm(beta - beta_old, "fro") ; 
    beta_old = beta ;
    U = Y ;
    for (int k = 0; k < p; k++) {
      U = U - X.slice(k) * beta(k) ;
    }
    pf = panel_factor(U, r)  ;
    F = as<arma::mat>(pf["factor"]) ;
    L = as<arma::mat>(pf["lambda"]) ; 
  }
  VNT = as<arma::mat>(pf["VNT"]) ; 
  arma::mat e = U - F * L.t() ;

  /* Storage */
  List result ;
  result["niter"] = niter ;
  result["beta"] = beta ;
  result["e"] = e ; 
  result["lambda"] = L ;
  result["factor"] = F ;
  result["VNT"] = VNT ;
  return(result)  ;
}


/* ******************* Main Function  *********************** */

   
/* Interactive Fixed Effects */
// [[Rcpp::export]]
List inter_fe (arma::mat Y,
               arma::cube X,
               int r,
               int force,
               arma::mat beta0, 
               double tol = 1e-5,
               int trends = 0
               ) {

  
  /* Dimensions */
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = X.n_slices ;
  int obs = T * N ;
  int niter ;
  arma::mat factor ;
  arma::mat lambda ;
  arma::mat VNT ;
  arma::mat beta ; 
  arma::mat U ;
  double mu ;
  double mu_Y ;
  arma::mat mu_X(p, 1) ;
  arma::mat alpha ;
  arma::mat alpha_Y(N, 1) ;
  arma::mat alpha_X(N, p) ;
  arma::mat xi ;
  arma::mat xi_Y(T, 1) ;
  arma::mat xi_X(T, p) ;
  arma::mat Ftd(T, T) ;
  double sigma2 ;
  double IC ;

  /* duplicate data */
  arma::mat YY = Y;
  arma::cube XX = X;
   
  /* grand mean */
  mu_Y  =  accu(YY)/obs ;
  YY  =  YY - mu_Y ;
  for (int i = 0; i < p; i++) {
    mu_X(i,0)  =  accu(XX.slice(i))/obs ;
    XX.slice(i) =  XX.slice(i)- mu_X(i,0) ;
  }
  
  /* unit fixed effects */
  if (force ==1 || force ==3 ){
    alpha_Y  =  mean(YY, 0).t() ; // colMeans, (N * 1) matrix
    YY  =  YY - repmat(alpha_Y.t(), T, 1) ; // (T * N) matrix 
    if (p > 0) {
      for (int i = 0; i < p; i++) {
        alpha_X.col(i)  = mean(XX.slice(i), 0).t(); // colMeans 
        XX.slice(i) = XX.slice(i) - repmat(alpha_X.col(i).t(), T, 1) ;
      }
    } 
  }
  
  /* time fixed effects  */
  if ( force == 2 || force == 3 ) {
    xi_Y  =  mean(YY, 1); //rowMeans, (N * 1) matrix
    YY  =  YY - repmat(xi_Y, 1, N);
    if (p > 0) {
      for (int i = 0; i < p; i++) {
        xi_X.col(i)  =  mean(XX.slice(i), 1) ; //rowMeans
        XX.slice(i) =  XX.slice(i) - repmat(xi_X.col(i), 1, N) ;
      }
    }  
  }
  
  /* trends and unit FEs: needs debugging */
  if (trends == 1) {
    arma::mat A(T, 1, arma::fill::ones) ; A = cumsum(A) ; 
    Ftd = DM(T) * A ; 
  } else if (trends == 2) {
    arma::mat A(T, 1, arma::fill::ones) ; A = cumsum(A) ;
    arma::mat AA = join_rows(A, A % A);
    Ftd  =  DM(T) * AA ;
  } else if (trends == 3) {
    arma::mat A(T, 1, arma::fill::ones) ; A = cumsum(A) ;
    arma::mat AA = join_rows(A, A % A);
    arma::mat AAA = join_rows(AA, A % A % A) ;
    Ftd  =  DM(T) * AAA ;
  }
  if (trends == 1 || trends ==2 || trends ==3) {
    arma::mat M_td  =  arma::eye<arma::mat>(T,T) -
      Ftd * inv(crossprod(Ftd, Ftd)) * Ftd.t() ;
    Y  =  M_td * Y ;
    if (p > 0) {
      for (int i = 0; i < p; i++) {
        X.slice(i) =  M_td * X.slice(i) ;
      }
    } 
  }
      
  /* Main Algorithm */ 
  if (p == 0) {
    if (r > 0) {
      List pf = panel_factor(YY, r)  ;
      factor = as<arma::mat>(pf["factor"]) ;
      lambda = as<arma::mat>(pf["lambda"]) ;
      VNT = as<arma::mat>(pf["VNT"]) ;
      U  =  YY - factor * lambda.t() ;
    } else {
      U = YY ;
    } 
  } else { 
    /* starting value:  the OLS/LSDV estimator */
    if (accu(abs(beta0))< 1e-10 || r==0) {  //
      beta0 = panel_est(XX, YY, arma::eye<arma::mat>(T,T)) ; //
    }
    if (r==0) {
      beta  =  beta0 ;
      for (int k = 0; k < p; k++) {
        U =  YY - XX.slice(k) * beta0(k,0);
      }
    } else if (r > 0) {  
      arma::mat invXX =  XXinv(XX) ;  // compute (X'X)^{-1}, outside beta iteration      
      List out  =  beta_iter(XX, invXX, YY, r, tol, beta0) ;
      beta  = as<arma::mat>(out["beta"]) ;
      factor  =  as<arma::mat>(out["factor"]) ;
      lambda  =  as<arma::mat>(out["lambda"]) ;
      VNT  =  as<arma::mat>(out["VNT"]) ;
      U  =  as<arma::mat>(out["e"]);
      niter = as<int>(out["niter"]) ;
    }
  } 
    
  /* save fixed effects */
  if (p == 0) {
    mu =  mu_Y;
    if (force ==1 || force ==3) {
      alpha =  alpha_Y ;
    }
    if (force == 2 || force == 3) {
      xi =  xi_Y;
    }   
  } else { // with covariates
    mu  =  mu_Y - crossprod(mu_X, beta)(0,0) ;
    if (force ==1 || force == 3) {
      alpha  =  alpha_Y - alpha_X * beta ; 
    }
    if (force == 2 || force == 3) {
      xi  =  xi_Y - xi_X * beta  ;
    } 
  }

  /* sigma2 and IC */
  sigma2 = trace(U * U.t())/ (N * T - r * (N + T) + pow(r,2) - p - trends * N) ;
  
  IC = log(sigma2) + (r * ( N + T ) - pow(r,2)) * log ( N * T ) / ( N * T ) ;
    
  //-------------------------------#
  // Storage
  //-------------------------------# 

  List output ;
  output["mu"] = mu ;
  if (p > 0)  {
    output["beta"] = beta ;
  }
  if (r > 0) {
    output["factor"] = factor ;
    output["lambda"] = lambda ;
    output["VNT"] = VNT ;
  }
  if (p > 0 && r > 0) {
    output["niter"] = niter ;
  }
  if (force ==1 || force == 3) {
    output["alpha"] = alpha ;
  }
  if (force ==2 || force == 3) {
    output["xi"] = xi ;
  }
  output["residuals"] = U ;
  output["sigma2"] = sigma2 ;
  output["IC"] = IC; 
  return(output);

 
}

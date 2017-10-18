# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp ;

/* ******************* Useful Functions  *********************** */

/* cross product */
arma::mat crossprod (arma::mat x, arma::mat y) {
  return(x.t() * y);
}

/* Expectation :E if Iij==0, Wij=FEij */
arma::mat E_adj (arma::mat E, arma::mat FE,
                 arma::mat I) {
  int T = E.n_rows ;
  int N = E.n_cols ;
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      if (I(i, j) == 0) {
        E(i, j) = FE(i, j) ;
      }
    }
  }
  return(E) ;
}

/* reset FEij=0 if Iij==0 */
arma::mat FE_adj (arma::mat FE, arma::mat I) {
  int T = FE.n_rows ;
  int N = FE.n_cols ;
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      if (I(i, j) == 0) {
        FE(i, j) = 0 ;
      }
    }
  }
  return(FE) ;
}

/* estimate standard deviation for a matrix */
//double std_dev (arma::mat X) {
//  int T = X.n_rows;
//  int N = X.n_cols;
//  double mu = arma::accu(X)/(N*T) ;
//  arma::mat xx(T, N);
//  for(int i=0; i<T; i++) {
//    for(int j=0; j<N; j++) {
//      xx(i,j)=pow((X(i,j)-mu),2);
//    }
//  }
//  double sd=arma::accu(xx)/(N*T-1);
//  return(sqrt(sd));
//}
 
/* estimate standard deviation for a matrix for ub data */
//double std_dev_ub (arma::mat X, arma::mat I) {
//  int T = X.n_rows;
//  int N = X.n_cols;
//  double mu = arma::accu(X)/arma::accu(I) ;
//  arma::mat xx(T, N, arma::fill::zeros);
//  for(int i=0; i<T; i++) {
//    for(int j=0; j<N; j++) {
//      if(I(i,j)!=0){
//        xx(i,j)=pow((X(i,j)-mu),2);
//      }
//    }
//  }
//  double sd=arma::accu(xx)/(arma::accu(I)-1);
//  return(sqrt(sd));
//}

/* calculate mean square of a matrix */
//double MS (arma::mat X) {
//  int T = X.n_rows;
//  int N = X.n_cols;
//  arma::mat xx(T, N);
//  for(int i=0; i<T; i++) {
//    for(int j=0; j<N; j++) {
//      xx(i,j)=pow((X(i,j)),2);
//    }
//  }
//  double mse = arma::accu(xx)/(N*T);
//  return(mse);
//}

/* calculate mean of square of a matrix for ub data */
//double MS_ub (arma::mat X, arma::mat I) {
//  int T = X.n_rows;
//  int N = X.n_cols;
//  arma::mat xx(T, N, arma::fill::zeros);
//  for(int i=0; i<T; i++) {
//    for(int j=0; j<N; j++) {
//      if(I(i,j)!=0){
//        xx(i,j)=pow((X(i,j)),2);  
//      }
//    }
//  }
//  double mse = arma::accu(xx)/arma::accu(I);
//  return(mse);
//}

/* adjust unbalanced data */
// [[Rcpp::export]]
arma::mat data_ub_adj (arma::mat I_data, arma::mat data) {
  int count = I_data.n_rows ;
  //int total = data.n_rows ;
  int nov = data.n_cols ;
  arma::mat data_adj(count,nov) ;
  data_adj.fill(arma::datum::nan) ;
  int j = 0;
  for(int i=0; i<count; i++){
    if(I_data(i,0)==1){
      data_adj.row(i) = data.row(j);
      j++;
    }
  }
  return(data_adj);
}

/* Three dimensional matrix inverse */
// [[Rcpp::export]]
arma::mat XXinv (arma::cube X) { 
  int p = X.n_slices ;
  arma::mat xx(p, p) ;
  for (int k = 0; k < p; k++) {
    for (int m = 0; m < p; m++) {
      if (k > m) {
        xx(k, m) = xx(m, k);
      }
      else {
        xx(k, m) = trace(crossprod(X.slice(k), X.slice(m))) ;
      }
    }
  } 
  return(inv(xx)) ;
}

/* ******************* Subsidiary Functions  *********************** */

/* Obtain OLS panel estimate */
// [[Rcpp::export]]
arma::mat panel_est (arma::cube X, arma::mat Y, arma::mat MF) {
  int p = X.n_slices ;
  arma::mat xx(p, p, arma::fill::zeros);
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

/* Obtain beta given interactive fe */
// [[Rcpp::export]]
arma::mat panel_beta (arma::cube X, arma::mat xxinv,
                      arma::mat Y, arma::mat FE) {
  int p = X.n_slices ; 
  arma::mat xy(p, 1, arma::fill::zeros) ;
  for (int k = 0; k < p; k++) {
    xy(k) = trace(crossprod(X.slice(k), (Y - FE))) ;
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
  arma::mat FE (T, N, arma::fill::zeros) ;
  arma::mat VNT(r, r) ;
  arma::mat U ;
  arma::vec s ;
  arma::mat V ; 
  if (T < N) { 
    arma::mat EE = E * E.t() /(N * T) ;
    arma::svd( U, s, V, EE) ;
    factor = U.head_cols(r) * sqrt(double(T)) ;
    lambda = E.t() * factor/T ;
    VNT = diagmat(s.head_rows(r)) ;
  } 
  else {
    arma::mat EE = E.t() * E / (N * T) ;
    svd(U, s, V, EE) ;
    lambda = U.head_cols(r) * sqrt(double(N)) ;
    factor = E * lambda / N ;
    VNT = diagmat(s.head_rows(r)) ;
  }
  FE = factor * lambda.t() ;
  List result ;
  result["lambda"] = lambda ;
  result["factor"] = factor ;
  result["VNT"] = VNT ;
  result["FE"] = FE ;
  return(result) ;
  
}

/* Obtain factors and loading given error for ub data */
// [[Rcpp::export]]
List panel_factor_ub (arma::mat E, arma::mat I, int r, double tolerate) {
  int T = E.n_rows ;
  int N = E.n_cols ;
  int niter = 0;
  double L_norm = 1.0 ;
  arma::mat F(T, r) ;
  arma::mat L(N, r) ;
  arma::mat L_old(N, r, arma::fill::zeros) ;
  arma::mat FE_0(T, N) ;
  arma::mat FE ;
  arma::mat E_use(T, N) ;
  arma::mat VNT(r, r) ;

  List pf = panel_factor(E, r)  ;
  F = as<arma::mat>(pf["factor"]) ;
  L = as<arma::mat>(pf["lambda"]) ;

  while ( (niter<500) && (L_norm>tolerate) ) {
    niter++ ;
    FE_0 = F * L.t() ;
    E_use = E_adj(E, FE_0, I) ;
    pf = panel_factor(E_use, r)  ;
    F = as<arma::mat>(pf["factor"]) ;
    L = as<arma::mat>(pf["lambda"]) ;
    VNT = as<arma::mat>(pf["VNT"]) ;
    L_norm = arma::norm(L - L_old, "fro") ;
    L_old = L ;
  }
  FE = FE_adj(F*L.t(), I) ;

  List result ;
  result["niter"] = niter ;
  result["lambda"] = L ;
  result["factor"] = F ;
  result["VNT"] = VNT ;
  result["FE"] = FE ;
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
  int b_r = beta0.n_rows ; 
  double beta_norm = 1.0 ;
  arma::mat beta ;
  if (b_r != p) {
      beta.zeros(p, 1) ;
  } else {
      beta = beta0 ;
  }
  arma::mat beta_old = beta ;
  arma::mat VNT(r, r) ;
  arma::mat FE ;

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
    FE = F * L.t() ;
    beta = panel_beta(X, xxinv, Y, FE) ;
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

/* Main iteration for beta: unbalanced */
// [[Rcpp::export]]
List beta_iter_ub (arma::cube X,
                   arma::mat xxinv,
                   arma::mat Y,
                   arma::mat I,
                   int r,
                   double tolerate,
                   arma::mat beta0) { 
  int p = X.n_slices ;
  int b_r = beta0.n_rows ; 
  double beta_norm = 1.0 ;
  arma::mat beta ;
  if (b_r != p) {
      beta.zeros(p, 1) ;
  } else {
      beta = beta0 ;
  }
  arma::mat beta_old = beta ;  
  arma::mat VNT(r, r) ;
  arma::mat FE ;
  arma::mat FE_use ;
  arma::mat U_use ;

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
    FE = F * L.t() ;
    // set missing value = 0
    FE_use = FE_adj(FE = FE, I = I) ;
    beta = panel_beta(X, xxinv, Y, FE_use) ;
    beta_norm = arma::norm(beta - beta_old, "fro") ; 
    beta_old = beta ;
    U = Y ;
    for (int k = 0; k < p; k++) {
      U = U - X.slice(k) * beta(k) ;
    }
    // Expectation for missing value
    U_use = E_adj(U, FE, I) ;
    pf = panel_factor(U_use, r)  ;
    F = as<arma::mat>(pf["factor"]) ;
    L = as<arma::mat>(pf["lambda"]) ; 
  }
  VNT = as<arma::mat>(pf["VNT"]) ; 
  FE = F * L.t() ;
  FE_use = FE_adj(FE = FE, I = I) ;
  arma::mat e = U - FE_use ;

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

/* Interactive Fixed Effects */
// [[Rcpp::export]]
List inter_fe (arma::mat Y,
               arma::cube X,
               int r,
               int force,
               arma::mat beta0, 
               double tol = 1e-5
               ) { 
  /* Dimensions */
  int b_r = beta0.n_rows ; 
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
  double sigma2 ;
  double IC ;
  //arma::mat FE(T, N, arma::fill::zeros) ;
  arma::mat invXX ;

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

  /* check if XX has enough variation */
  int p1=p; 
  arma::mat X_invar(p, 1, arma::fill::zeros); // =1 if invar

  int j = 0;
  for(int i=0; i<p1; i++){
    if (arma::accu(abs(XX.slice(i))) < 1e-5) {
      XX.shed_slice(i);
      mu_X.shed_row(i);
      alpha_X.shed_col(i);
      xi_X.shed_col(i);
      X_invar(j,0) = 1;
      i--;
      p1--;
    }
    j++;
  }

  int validX = 1;
  if(p1==0){
    validX = 0;
  }
  else {
    invXX =  XXinv(XX) ;
  }
 
  /* Main Algorithm */ 
  if (p1 == 0) {
    if (r > 0) {
      List pf = panel_factor(YY, r)  ;
      factor = as<arma::mat>(pf["factor"]) ;
      lambda = as<arma::mat>(pf["lambda"]) ;
      VNT = as<arma::mat>(pf["VNT"]) ;
      U  =  YY - factor * lambda.t() ;
    } 
    else {
      U = YY ;
    } 
  } 
  else {
    /* starting value:  the OLS/LSDV estimator */
    if (accu(abs(beta0))< 1e-10 || r==0 || b_r != p1 ) {  //
      beta0 = panel_beta(XX, invXX, YY, arma::zeros<arma::mat>(T,N)); //
    }
    if (r==0) {
      U = YY;
      beta  =  beta0 ;
      for (int k = 0; k < p1; k++) {
        U =  U - XX.slice(k) * beta(k,0);
      }
    } 
    else if (r > 0) {  
      // arma::mat invXX =  XXinv(XX) ;  // compute (X'X)^{-1}, outside beta iteration      
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
  if (p1 == 0) {
    
    mu =  mu_Y;
    
    if (force ==1 || force ==3) {
      alpha =  alpha_Y ;
    }
    if (force == 2 || force == 3) {
      xi =  xi_Y;
    }   
  } else { // with valid covariates
    
    mu  =  mu_Y - crossprod(mu_X, beta)(0,0) ;
    
    if (force ==1 || force == 3) {
      alpha  =  alpha_Y - alpha_X * beta ; 
    }
    if (force == 2 || force == 3) {
      xi  =  xi_Y - xi_X * beta  ;
    } 
  }

  /* sigma2 and IC */
  sigma2 = trace(U * U.t())/ (N * T - r * (N + T) + pow(double(r),2) - p1 ) ;
  
  IC = log(sigma2) + (r * ( N + T ) - pow(double(r),2) + p1) * log ( double(N * T) ) / ( N * T ) ;
    
  //-------------------------------#
  // Storage
  //-------------------------------# 

  List output ;
  
  output["mu"] = mu ;
  // output["p1"] = p1 ;  
  
  if(p>0){
    // output["beta_valid"] = beta ;
    arma::mat beta_total(p,1);
    // arma::mat beta_tot(p,1);

    if(p>p1) {
      int j4= 0;
      for(int i=0; i<p; i++) {
        if(X_invar(i,0)==1) {
          beta_total(i,0) = arma::datum::nan;
          // beta_tot(i,0) = 0;
        }
        else {
          beta_total(i,0) = beta(j4,0);
          // beta_tot(i,0) = beta(j4,0);
          j4++;
        }
      }
    }
    else {
      beta_total = beta;
      // beta_tot = beta;
    }
    output["beta"] = beta_total;
    // output["beta_tot"] = beta_tot;
  }   

  if (r > 0) {
    output["factor"] = factor ;
    output["lambda"] = lambda ;
    output["VNT"] = VNT ;
    //FE = factor * lambda.t() ;
    //output["FE"] = FE ;
  }
  if ((p1 > 0) && (r > 0)) {
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
  output["IC"] = IC ;
  output["validX"] = validX ;
  return(output);
}

/* Interactive Fixed Effects: ub */
// [[Rcpp::export]]
List inter_fe_ub (arma::mat Y,
                  arma::cube X,
                  arma::mat I,
                  int r,
                  int force,
                  arma::mat beta0, 
                  double tol = 1e-5
                  ) {
  
  /* Dimensions */
  int b_r = beta0.n_rows ; 
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = X.n_slices ;
  double obs = accu(I) ;
  int niter ;
  arma::mat factor ;
  arma::mat lambda ;
  arma::mat FE_0 ;
  //arma::mat FE(T, N, arma::fill::zeros) ;
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
  double sigma2 ;
  double IC ;

  /* duplicate data */
  arma::mat subX(T, N) ;
  arma::mat YY = Y;
  arma::cube XX = X;
   
  /* grand mean */
  mu_Y  =  accu(YY)/obs ;
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      if (I(i, j) == 1) {
        YY(i, j) = YY(i, j) - mu_Y;
      }
    }
  }

  for (int k = 0; k < p; k++) {
    mu_X(k,0)  =  accu(XX.slice(k))/obs ;
    subX = XX.slice(k) ;
    for (int i = 0; i < T; i++) {
      for (int j = 0; j < N; j++) {
        if (I(i, j) == 1) {
          subX(i, j) = subX(i, j) - mu_X(k,0);
        }
      }
    }
    XX.slice(k) =  subX ;
  }
  
  /* unit fixed effects */
  if (force ==1 || force ==3 ) {
    alpha_Y  =  sum(YY, 0).t() / sum(I, 0).t() ; // colMeans, (N * 1) matrix
    for (int i = 0; i < T; i++) {
      for (int j = 0; j < N; j++) {
        if (I(i, j) == 1) {
          YY(i, j) = YY(i, j) - alpha_Y(j, 0);
        }
      }
    }    
    if (p > 0) {
      for (int k = 0; k < p; k++) {
        alpha_X.col(k)  = sum(XX.slice(k), 0).t() / sum(I, 0).t() ; // colMeans 
        subX = XX.slice(k) ; 
        for (int i = 0; i < T; i++) {
          for (int j = 0; j < N; j++) {
            if (I(i, j) == 1) {
              subX(i, j) = subX(i, j) - alpha_X(j, k);
            }
          }
        }
        XX.slice(k) = subX ;        
      }
    } 
  }
  
  /* time fixed effects  */
  if ( force == 2 || force == 3 ) {
    xi_Y  =  sum(YY, 1) / sum(I, 1); //rowMeans, (T * 1) matrix
    for (int i = 0; i < T; i++) {
      for (int j = 0; j < N; j++) {
        if (I(i, j) == 1) {
          YY(i, j) = YY(i, j) - xi_Y(i, 0);
        }
      }
    }   
    if (p > 0) {
      for (int k = 0; k < p; k++) {
        xi_X.col(k)  =  sum(XX.slice(k), 1) / sum(I, 1); //rowMeans
        subX = XX.slice(k) ; 
        for (int i = 0; i < T; i++) {
          for (int j = 0; j < N; j++) {
            if (I(i, j) == 1) {
              subX(i, j) = subX(i, j) - xi_X(i, k);
            }
          }
        }
        XX.slice(k) = subX ;          
      }
    }  
  }

  /* check if XX has enough variation */
  int p1 = p; 
  arma::mat X_invar(p, 1, arma::fill::zeros); // =1 if invar

  int j = 0;

  for(int i=0; i<p1; i++){
    if (arma::accu(abs(XX.slice(i))) < 1e-5) {
      XX.shed_slice(i);
      mu_X.shed_row(i);
      alpha_X.shed_col(i);
      xi_X.shed_col(i);
      X_invar(j,0) = 1;
      i--;
      p1--;
    }
    j++;
  }

  int validX = 1;
  if(p1==0){
    validX = 0;
  }

  /* Main Algorithm */ 
  if (p1 == 0) {
    if (r > 0) {
      List pf = panel_factor_ub(YY, I, r, tol)  ;
      factor = as<arma::mat>(pf["factor"]) ;
      lambda = as<arma::mat>(pf["lambda"]) ;
      VNT = as<arma::mat>(pf["VNT"]) ;
      FE_0 = FE_adj(factor * lambda.t(), I) ;
      U  =  YY - FE_0 ;
    } 
    else {
      U = YY ;
    } 
  } 
  else {
    /* starting value:  the OLS estimator */
    arma::mat invXX = XXinv(XX) ; // compute (X'X)^{-1}, outside beta iteration 
    if (accu(abs(beta0))< 1e-10 || r==0 || b_r != p1) {  //
      beta0 = panel_beta(XX, invXX, YY, arma::zeros<arma::mat>(T,N)) ; //
    }
    if (r==0) {
      beta  =  beta0 ;
      U = YY;
      for (int k = 0; k < p1; k++) {
        U =  U - XX.slice(k) * beta(k,0);
      }
    } 
    else if (r > 0) {       
      List out  =  beta_iter_ub(XX, invXX, YY, I, r, tol, beta0) ;
      beta  =  as<arma::mat>(out["beta"]) ;
      factor  =  as<arma::mat>(out["factor"]) ;
      lambda  =  as<arma::mat>(out["lambda"]) ;
      VNT  =  as<arma::mat>(out["VNT"]) ;
      /* U = YY;
      for (int k = 0; k < p1; k++) {
        U =  U - XX.slice(k) * beta(k,0);
      }
      List pf = panel_factor_ub(U, I, r, tol)  ;
      factor = as<arma::mat>(pf["factor"]) ;
      lambda = as<arma::mat>(pf["lambda"]) ;
      VNT = as<arma::mat>(pf["VNT"]) ;
      FE = FE_adj(factor * lambda.t(), I) ;
      U  =  U - FE ; */
      U  =  as<arma::mat>(out["e"]) ;
      niter = as<int>(out["niter"]) ;
    }
  } 
    
  /* save fixed effects */
  if (p1 == 0) {
    
    mu =  mu_Y;
    
    if (force ==1 || force ==3) {
      alpha =  alpha_Y ;
    }
    if (force == 2 || force == 3) {
      xi =  xi_Y;
    }   
  } 
  else { // with valid covariates
    
    mu  =  mu_Y - crossprod(mu_X, beta)(0,0) ;
    
    if (force ==1 || force == 3) {
      alpha  =  alpha_Y - alpha_X * beta ; 
    }
    if (force == 2 || force == 3) {
      xi  =  xi_Y - xi_X * beta  ;
    } 
  }

  /* sigma2 and IC */
  sigma2 = trace(U * U.t())/ (obs - r * (N + T) + pow(double(r),2) - p1 ) ;
  
  IC = log(sigma2) + (r * ( N + T ) - pow(double(r),2) + p1)
   * log ( obs ) / ( obs ) ;
    
  //-------------------------------#
  // Storage
  //-------------------------------# 

  List output ;
  output["mu"] = mu ;
  // output["p1"] = p1 ;
  // output["beta_valid"] = beta ;
  
  
  if(p>0){
    arma::mat beta_total(p,1);
    if(p>p1) {
      int j4= 0;
      for(int i=0; i<p; i++) {
        if(X_invar(i,0)==1) {
          beta_total(i,0) = arma::datum::nan;
        }
        else {
          beta_total(i,0) = beta(j4,0);
          j4++;
        }
      }
    }
    else {
      beta_total = beta;
    }
    output["beta"] = beta_total;
  }      

  if (r > 0) {
    output["factor"] = factor ;
    output["lambda"] = lambda ;
    //FE = FE_adj(factor * lambda.t(), I) ;
    //output["FE"] = FE;
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
  output["validX"] = validX;
  return(output);
 
}


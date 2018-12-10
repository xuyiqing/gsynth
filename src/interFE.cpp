# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp ;

/* ******************* Useful Functions  *********************** */

/* cross product */
arma::mat crossprod (arma::mat x, arma::mat y) {
  return(x.t() * y);
}

/* Expectation :E if Iij==0, Eij=FEij */
arma::mat E_adj (arma::mat E, arma::mat FE,
                 arma::mat I) {
  int T = E.n_rows ;
  int N = E.n_cols ;
  arma::mat EE = E ;
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      if (I(i, j) == 0) {
        EE(i, j) = FE(i, j) ;
      }
    }
  }
  return(EE) ;
}

/* reset FEij=0 if Iij==0 , for residuals or IC*/
arma::mat FE_adj (arma::mat FE, arma::mat I) {
  int T = FE.n_rows ;
  int N = FE.n_cols ;
  arma::mat FEE = FE ;
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      if (I(i, j) == 0) {
        FEE(i, j) = 0 ;
      }
    }
  }
  return(FEE) ;
}

/* drop values if Iij == 1 */
arma::mat FE_missing (arma::mat FE, arma::mat I) {
  int T = FE.n_rows ;
  int N = FE.n_cols ;
  arma::mat FEE = FE ;
  for (int i = 0; i < T; i++) {
    for (int j = 0; j < N; j++) {
      if (I(i, j) == 1) {
        FEE(i, j) = 0 ;
      }
    }
  }
  return(FEE) ;
}

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
  arma::mat xx(p, p, arma::fill::zeros) ;
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

/* unbalanced panel: response demean function */
// [[Rcpp::export]]
List Y_demean (arma::mat Y, int force) {
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  double mu_Y = 0 ;
  arma::mat alpha_Y(N, 1, arma::fill::zeros) ; 
  arma::mat xi_Y(T, 1, arma::fill::zeros) ;
  arma::mat YY = Y ;

  mu_Y  =  accu(YY)/(N*T) ;
  if (force == 0) {
    YY = YY - mu_Y ;
  } 
  
  /* unit fixed effects */
  if (force == 1) {
    alpha_Y  =  mean(YY, 0).t() ; // colMeans, (N * 1) matrix
    YY  =  YY - repmat(alpha_Y.t(), T, 1) ; // (T * N) matrix  
  }
  
  /* time fixed effects  */
  if ( force == 2) {
    xi_Y  =  mean(YY, 1) ; //rowMeans, (N * 1) matrix
    YY  =  YY - repmat(xi_Y, 1, N);  
  }

  if (force == 3) {
    alpha_Y  =  mean(YY, 0).t() ;
    xi_Y  =  mean(YY, 1) ;
    YY  =  YY - repmat(alpha_Y.t(), T, 1) - repmat(xi_Y, 1, N) + mu_Y ;
  }

  List result ;
  result["mu_Y"] = mu_Y ;
  result["YY"] = YY ;
  if (force==1 || force==3) {
    result["alpha_Y"] = alpha_Y ;
  }
  if (force==2 || force==3) {
    result["xi_Y"] = xi_Y ;  
  }

  return(result) ;
}

/* estimate additive fe for unbalanced panel, without covariates */
// [[Rcpp::export]]
List fe_add (arma::mat alpha_Y,
             arma::mat xi_Y,
             double mu_Y,
             int T,
             int N,
             int force) {
  arma::mat FE_ad(T, N, arma::fill::zeros) ;
  double mu = 0 ;
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;

  mu =  mu_Y;
  if (force ==1 || force ==3) {
    alpha =  alpha_Y - mu_Y ;
  }
  if (force == 2 || force == 3) {
    xi =  xi_Y - mu_Y ;
  }   

  FE_ad = FE_ad + mu ;

  if (force ==1 || force == 3) {
    FE_ad = FE_ad + repmat(alpha.t(), T, 1) ;
  }
  if (force == 2 || force == 3) {
    FE_ad = FE_ad + repmat(xi, 1, N) ;
  }

  List result ;
  result["mu"] = mu ;
  result["FE_ad"] = FE_ad ;
  if (force ==1 || force == 3) {
    result["alpha"] = alpha ;
  }
  if (force == 2 || force == 3) {
    result["xi"] = xi ;
  }

  return(result) ;
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
  arma::mat factor(T, r, arma::fill::zeros) ;
  arma::mat lambda(N, r, arma::fill::zeros) ;
  arma::mat FE (T, N, arma::fill::zeros) ;
  arma::mat VNT(r, r, arma::fill::zeros) ;
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

/* Obtain interactive fe directly */
// [[Rcpp::export]]
arma::mat panel_FE (arma::mat E, double lambda) {
  int T = E.n_rows ;
  int N = E.n_cols ;
  int r = T ;
  if (T >= N) {
    r = N ;
  }

  arma::mat FE (T, N, arma::fill::zeros) ;
  arma::mat D(r, r, arma::fill::zeros) ;
  arma::mat U ;
  arma::vec s ;
  arma::mat V ;
  arma::svd( U, s, V, E) ;
  
  for (int i = 0; i < r; i++) {
    if (s(i) > lambda) {
      D(i, i) = s(i) - lambda ;
    } else {
      D(i, i) = 0 ;
    }
  }
  if (T >= N) {
    arma::mat UU = U.cols(0, r-1) ; 
    FE = UU * D * V.t() ;
  }
  else {
    arma::mat VV = V.cols(0, r-1) ;
    FE = U * D * VV.t() ;
  }
  return(FE) ;
}

/* factor analysis: mu add ife*/
// [[Rcpp::export]]
List ife (arma::mat E, 
          int force,
          int mc, // whether pac or mc method
          int r,
          double lambda
         ) {

  int T = E.n_rows ; 
  int N = E.n_cols ;

  arma::mat VNT(r, r) ;
  arma::mat EE(T, N, arma::fill::zeros) ;
  arma::mat FE_add_use(T, N, arma::fill::zeros) ;
  arma::mat FE_inter_use(T, N, arma::fill::zeros) ;
  arma::mat FE(T, N, arma::fill::zeros) ;

  arma::mat F(T, r, arma::fill::zeros) ;
  arma::mat L(N, r, arma::fill::zeros) ;

  double mu_E = 0 ;
  arma::mat alpha_E(N, 1, arma::fill::zeros) ;
  arma::mat xi_E(T, 1, arma::fill::zeros) ;

  double mu = 0 ;
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;

  List E_ad ;
  List E_fe_ad ;
  List pf ;

  E_ad = Y_demean(E, force) ;
  EE = as<arma::mat>(E_ad["YY"]) ;

  mu_E = as<double>(E_ad["mu_Y"]) ;
  if (force==1||force==3) {
    alpha_E = as<arma::mat>(E_ad["alpha_Y"]) ;
  }
  if (force==2||force==3) {
    xi_E = as<arma::mat>(E_ad["xi_Y"]) ;
  }

  E_fe_ad = fe_add(alpha_E, xi_E, mu_E,
                   T, N, force) ;

  FE_add_use = as<arma::mat>(E_fe_ad["FE_ad"]) ; // additive fe

  if (r > 0) {
    if (mc == 0) {
      pf = panel_factor(EE, r)  ;
      F = as<arma::mat>(pf["factor"]) ;
      L = as<arma::mat>(pf["lambda"]) ;
      VNT = as<arma::mat>(pf["VNT"]) ;
      FE_inter_use = F * L.t() ; // interactive fe
    }
    else {
      FE_inter_use = panel_FE(EE, lambda) ;
    }
  }

  FE = FE_add_use + FE_inter_use ;

  List result ;
    
  mu = as<double>(E_fe_ad["mu"]) ;
  result["mu"] = mu ;
  result["FE"] = FE ;
  if (force==1||force==3) {
    alpha = as<arma::mat>(E_fe_ad["alpha"]) ;
    result["alpha"] = alpha ;
  }
  if (force==2||force==3) {
    xi = as<arma::mat>(E_fe_ad["xi"]) ;
    result["xi"] = xi ;
  }
  if (r > 0) {
    if (mc == 0) {
      result["lambda"] = L ;
      result["factor"] = F ;
      result["VNT"] = VNT ;
    }
    result["FE_inter_use"] = FE_inter_use ;
  }
  return(result) ;
}

/* Obtain additive fe for ub data; assume r=0, without covar */
// [[Rcpp::export]]
List fe_ad_iter (arma::mat Y,
                 arma::mat Y0,
                 arma::mat I,
                 int force,
                 double tolerate) {
  
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  double mu = 0 ;
  double dif = 1.0 ;
  int niter = 0 ;

  arma::mat fit = Y0 ; // initial value
  arma::mat fit_old = Y0 ;
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;
  double mu_Y = 0 ;
  arma::mat alpha_Y(N, 1, arma::fill::zeros) ;
  arma::mat xi_Y(T, 1, arma::fill::zeros) ;

  arma::mat e(T, N, arma::fill::zeros) ; // residual

  arma::mat YY = Y ;

  List Y_ad ;
  List Y_fe_ad ;

  while (dif > tolerate && niter <= 500) {

    YY =  E_adj (Y, fit, I) ; // e step: expeactation
    
    Y_ad = Y_demean(YY, force) ; 
    mu_Y = as<double>(Y_ad["mu_Y"]) ;
    if (force==1||force==3) {     
      alpha_Y = as<arma::mat>(Y_ad["alpha_Y"]) ;
    }
    if (force==2||force==3) {
      xi_Y = as<arma::mat>(Y_ad["xi_Y"]) ;
    }
    Y_fe_ad = fe_add(alpha_Y, xi_Y, mu_Y, T, N, force) ; // m step: estimate fe
    
    fit = as<arma::mat>(Y_fe_ad["FE_ad"]) ;
    
    dif = arma::norm(fit - fit_old, "fro")/arma::norm(fit_old, "fro") ;
    fit_old = fit ;

    niter = niter + 1 ;
  }

  e = FE_adj(YY - fit, I) ;

  List result;
  mu = as<double>(Y_fe_ad["mu"]) ;
  result["mu"] = mu ;
  result["fit"] = fit ;
  result["niter"] = niter ;
  result["e"] = e ;

  if (force==1||force==3) {
    alpha = as<arma::mat>(Y_fe_ad["alpha"]) ;
    result["alpha"] = alpha ;
  }
  if (force==2||force==3) {
    xi = as<arma::mat>(Y_fe_ad["xi"]) ;
    result["xi"] = xi ;
  }
  return(result) ;
}


/* Obtain additive fe for ub data; assume r=0, with covariates */
// [[Rcpp::export]]
List fe_ad_covar_iter (arma::cube XX,
                       arma::mat xxinv,
                       arma::mat Y,
                       arma::mat Y0,
                       arma::mat I,
                       arma::mat beta0,
                       int force,
                       double tolerate) {
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = XX.n_slices ;
  double dif = 1.0 ;
  int niter = 0 ;
  
  arma::mat beta(p, 1, arma::fill::zeros) ;
  if (beta0.n_rows == p) {
    beta = beta0 ;
  }

  double mu = 0 ;

  arma::mat fit(T, N, arma::fill::zeros) ;
  arma::mat fit_old(T, N, arma::fill::ones) ;

  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;

  arma::mat YY = Y ;

  arma::mat U(T, N, arma::fill::zeros) ; 
  arma::mat e(T, N, arma::fill::zeros) ; // residual

  List ife_inner ;

  arma::mat covar_fit(T, N, arma::fill::zeros) ;
  for (int i = 0; i < p; i++) {
    covar_fit = covar_fit + XX.slice(i) * beta(i) ;
  }

  fit = Y0 ;
  fit_old = fit ;

  arma::mat FE = fit - covar_fit ; // initial fixed effects

  while (dif > tolerate && niter <= 500) {

    YY =  E_adj (Y, fit, I) ; // e-step: expectation
    
    // m1: estimate beta
    beta = panel_beta(XX, xxinv, YY, FE) ;
    covar_fit.zeros() ;
    for (int i = 0; i < p; i++) {
      covar_fit = covar_fit + XX.slice(i) * beta(i) ;
    }

    // m2: estimate interactive fe, additive fe, and mu
    U = E_adj (YY - covar_fit, FE, I) ;
    ife_inner = ife(U, force, 0, 0, 0) ;
    FE = as<arma::mat>(ife_inner["FE"]) ; // overall fe 

    fit = covar_fit + FE ;

    dif = arma::norm(fit - fit_old, "fro")/arma::norm(fit_old, "fro") ;
    fit_old = fit ;

    niter = niter + 1 ;
  }

  e = FE_adj(YY - fit, I) ;

  List result ;
  mu = as<double>(ife_inner["mu"]) ;
  result["mu"] = mu ;
  result["fit"] = fit ;
  result["niter"] = niter ;
  result["e"] = e ;
  if (p>0) {
    result["beta"] = beta ;
  }
  if (force==1||force==3) {
    alpha = as<arma::mat>(ife_inner["alpha"]) ;
    result["alpha"] = alpha ;
  }
  if (force==2||force==3) {
    xi = as<arma::mat>(ife_inner["xi"]) ;
    result["xi"] = xi ;
  }
  return(result) ;
}

/* Obtain additive fe for ub data; assume r>0 but p=0*/
// [[Rcpp::export]]
List fe_ad_inter_iter (arma::mat Y,
                       arma::mat Y0,
                       arma::mat I,
                       int force,
                       int mc, // whether pac or mc method
                       int r,
                       double lambda,
                       double tolerate
                       ) {
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  double mu = 0 ;
  double dif = 1.0 ;
  int niter = 0 ;
  int validF = 1 ; // whether has a factor structure

  arma::mat VNT(r, r) ;
  arma::mat FE_inter_use(T, N, arma::fill::zeros) ;
  arma::mat fit(T, N, arma::fill::zeros) ;
  arma::mat fit_old(T, N, arma::fill::ones) ;
  arma::mat U = FE_adj(Y - Y0, I) ;
  arma::mat e(T, N, arma::fill::zeros) ; // residual
  
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;

  arma::mat F(T, r, arma::fill::zeros) ;
  arma::mat L(N, r, arma::fill::zeros) ;


  arma::mat YY = Y ;

  List pf ;
  List ife_inner ;

  // initial value for ife
  if (mc == 0) {
    pf = panel_factor(U, r)  ;
    F = as<arma::mat>(pf["factor"]) ;
    L = as<arma::mat>(pf["lambda"]) ;
    FE_inter_use = F * L.t() ; // interactive fe
  }
  else {
    FE_inter_use = panel_FE(U, lambda) ;
  }
  fit = Y0 + FE_inter_use ;
  fit_old = fit ;

  while (dif > tolerate && niter <= 500) {
    
    YY = E_adj(Y, fit, I) ; // e-step: expectation
    
    if (mc == 0) {
      ife_inner = ife(YY, force, 0, r, 0) ;
    }
    else {
      ife_inner = ife(YY, force, 1, 1, lambda) ;
    }
    fit = as<arma::mat>(ife_inner["FE"]) ; // new overall fe

    dif = arma::norm(fit - fit_old, "fro")/arma::norm(fit_old, "fro") ;
    fit_old = fit ;

    niter = niter + 1 ;
  }
  e = FE_adj(YY - fit, I) ;
  FE_inter_use = as<arma::mat>(ife_inner["FE_inter_use"]) ;
  if (arma::accu(abs(FE_inter_use)) < 1e-10) {
    validF = 0 ;
  }

  List result;
  mu = as<double>(ife_inner["mu"]) ;
  result["mu"] = mu ;
  result["niter"] = niter ;
  result["fit"] = fit ;
  result["e"] = e ;
  result["validF"] = validF ;
  if (force==1||force==3) {
    alpha = as<arma::mat>(ife_inner["alpha"]) ;
    result["alpha"] = alpha ;
  }
  if (force==2||force==3) {
    xi = as<arma::mat>(ife_inner["xi"]) ;
    result["xi"] = xi ;
  }
  if (mc == 0) {
    result["lambda"] = L ;
    result["factor"] = F ;
    result["VNT"] = VNT ;
  }
  return(result) ;
}

/* Obtain additive fe for ub data; assume r>0 p>0*/
// [[Rcpp::export]]
List fe_ad_inter_covar_iter (arma::cube XX,
                             arma::mat xxinv,
                             arma::mat Y,
                             arma::mat Y0,
                             arma::mat I,
                             arma::mat beta0, 
                             int force,
                             int mc, // whether pac or mc method
                             int r,
                             double lambda,
                             double tolerate
                             ) {
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = XX.n_slices ;
  double dif = 1.0 ;
  int niter = 0 ;
  int validF = 1 ;

  arma::mat beta(p, 1, arma::fill::zeros) ;
  if (beta0.n_rows == p) {
    beta = beta0 ;
  }

  double mu = 0 ;
  arma::mat VNT(r, r) ;
  arma::mat FE_inter_use(T, N, arma::fill::zeros) ; // ife

  arma::mat fit(T, N, arma::fill::zeros) ;
  arma::mat fit_old(T, N, arma::fill::ones) ;
  arma::mat U = FE_adj(Y - Y0, I) ;

  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;

  arma::mat YY = Y ;

  arma::mat e(T, N, arma::fill::zeros) ; // residual

  List ife_inner ;
  List pf ;

  arma::mat F ;
  arma::mat L ;

  // initial value for ife
  if (mc == 0) {
    pf = panel_factor(U, r)  ;
    F = as<arma::mat>(pf["factor"]) ;
    L = as<arma::mat>(pf["lambda"]) ;
    FE_inter_use = F * L.t() ; // interactive fe
  }
  else {
    FE_inter_use = panel_FE(U, lambda) ;
  }

  fit = Y0 + FE_inter_use ;
  fit_old = fit ;

  arma::mat covar_fit(T, N, arma::fill::zeros) ;
  for (int i = 0; i < p; i++) {
    covar_fit = covar_fit + XX.slice(i) * beta(i) ;
  }

  arma::mat FE = fit - covar_fit ;


  while (dif > tolerate && niter <= 500) {
    
    YY =  E_adj (Y, fit, I) ; // e-step: expectation
    
    // m1: estimate beta
    beta = panel_beta(XX, xxinv, YY, FE) ;
    
    covar_fit.zeros() ;
    for (int i = 0; i < p; i++) {
      covar_fit = covar_fit + XX.slice(i) * beta(i) ;
    }
    
    // m2: estimate interactive fe, additive fe, and mu
    U = E_adj (YY - covar_fit, FE, I) ;
    ife_inner = ife(U, force, mc, r, lambda) ;
    FE = as<arma::mat>(ife_inner["FE"]) ;

    fit = covar_fit + FE ; // overall fe */

    dif = arma::norm(fit - fit_old, "fro")/arma::norm(fit_old, "fro") ;
    fit_old = fit ;

    niter = niter + 1 ;
  }
  e = FE_adj(Y - fit, I) ;

  FE_inter_use = as<arma::mat>(ife_inner["FE_inter_use"]) ;
  if (arma::accu(abs(FE_inter_use)) < 1e-10) {
    validF = 0 ;
  }

  List result;
  mu = as<double>(ife_inner["mu"]) ;
  result["mu"] = mu ;
  result["niter"] = niter ;
  result["e"] = e ;
  result["beta"] = beta ;
  result["fit"] = fit ;
  result["validF"] = validF ;

  if (force==1||force==3) {
    alpha = as<arma::mat>(ife_inner["alpha"]) ;
    result["alpha"] = alpha ;
  }
  if (force==2||force==3) {
    xi = as<arma::mat>(ife_inner["xi"]) ;
    result["xi"] = xi ;
  }
  if (mc == 0) {
    L = as<arma::mat>(ife_inner["lambda"]) ;
    F = as<arma::mat>(ife_inner["factor"]) ;
    VNT = as<arma::mat>(ife_inner["VNT"]) ;
    result["lambda"] = L ;
    result["factor"] = F ;
    result["VNT"] = VNT ;
  }
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
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = X.n_slices ;
  int b_r = beta0.n_rows ; 
  double beta_norm = 1.0 ;
  arma::mat beta(p, 1, arma::fill::zeros) ;
  if (b_r == p) {
      beta = beta0 ;
  } // beta should have the same dimension as X, if not it will be reset to 0
  arma::mat beta_old = beta ;
  arma::mat VNT(r, r, arma::fill::zeros) ;
  arma::mat FE(T, N, arma::fill::zeros) ;

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
  int niter = 0 ;
  arma::mat factor ;
  arma::mat lambda ;
  arma::mat VNT ;
  arma::mat beta ; 
  arma::mat U ;
  double mu = 0 ;
  double mu_Y = 0 ;
  arma::mat mu_X(p, 1) ;
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat alpha_Y(N, 1) ;
  arma::mat alpha_X(N, p) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;
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
  if (p > 0) {
    for (int i = 0; i < p; i++) {
      mu_X(i,0)  =  accu(XX.slice(i))/obs ;
      XX.slice(i) =  XX.slice(i)- mu_X(i,0) ;
    }
  }
  
  /* unit fixed effects */
  if (force ==1 || force ==3 ) {
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

  // number of estimated parameters  
  // force = 0
  double np = r * (N + T) - pow(double(r),2) + p1 + 1;
  if (force == 1) {
    np = np + (N - 1) - r ;
  }
  else if (force == 2) {
    np = np + (T - 1) - r ;
  } 
  else if (force == 3) {
    np = np + (N - 1) + (T - 1) - 2 * r ;
  }

  sigma2 = trace(U * U.t())/ (N * T - np) ;
  
  IC = log(sigma2) + np * log ( double(N * T) ) / ( N * T ) ;
  
  // PC criterion in Li 2018
  // mean squared error
  double mse = trace(U * U.t()) / (N * T) ;

  double m1 = 0 ;
  if (N < 60) {
    m1 = 60 - N ;
  }
  double m2 = 0 ;
  if (T < 60) {
    m2 = 60 - T ;
  }

  double C = (N + m1) * (T + m2) / (N * T) ;
  double PC = mse + r * sigma2 * C * (N + T) / (N * T) * log (double(N * T)/double(N + T)) ;
    
  //-------------------------------#
  // Storage
  //-------------------------------# 

  List output ;
  
  output["mu"] = mu ;
  // output["p1"] = p1 ; 
 
  if(p>0) {
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
  output["PC"] = PC ;
  output["validX"] = validX ;
  return(output);
}

/* Interactive Fixed Effects: ub */
// [[Rcpp::export]]
List inter_fe_ub (arma::mat Y,
                  arma::mat Y0,
                  arma::cube X,
                  arma::mat I,
                  arma::mat beta0, 
                  int r, // r > 0, the outcome has a factor-type fixed effect; r = 0 else
                  int force,
                  double tol = 1e-5
                  ) {
  
  /* Dimensions */
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = X.n_slices ;
  double obs = accu(I) ;
  int niter = 0 ;
  arma::mat factor ;
  arma::mat lambda ;
  arma::mat VNT ;
  arma::mat beta ; 
  arma::mat U ;
  double mu = 0 ;
  double mu_Y = 0 ;
  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;
  arma::mat fit(T, N, arma::fill::zeros) ;
  double sigma2 = 0 ;
  double IC = 0 ;

  arma::mat invXX ;

  /* duplicate data */
  arma::mat YY = Y;
  arma::cube XX = X;

  /* check if XX has enough variation */
  int p1 = p; 
  arma::mat X_invar(p, 1, arma::fill::zeros); // =1 if invar

  int j = 0;

  for(int i=0; i<p1; i++){
    if (arma::accu(abs(XX.slice(i))) < 1e-5) {
      XX.shed_slice(i);
      X_invar(j,0) = 1;
      i--;
      p1--;
    }
    j++;
  }

  int validX = 1 ;
  if(p1==0){
    validX = 0 ;
    if (force == 0 && r == 0) { // no covariate and force == 0 and r == 0 
      mu_Y = accu(YY)/obs ;
      mu = mu_Y ;
      YY = FE_adj(YY - mu_Y, I) ;
    }
  }

  /* Main Algorithm */ 
  if (p1 == 0) {
    if (r > 0) {
      // add fe ; inter fe ; iteration
      List fe_ad_inter = fe_ad_inter_iter(YY, Y0, I, force, 0, r, 0, tol) ;
      mu = as<double>(fe_ad_inter["mu"]) ;
      U = as<arma::mat>(fe_ad_inter["e"]) ;
      fit = as<arma::mat>(fe_ad_inter["fit"]) ;

      factor = as<arma::mat>(fe_ad_inter["factor"]) ;
      lambda = as<arma::mat>(fe_ad_inter["lambda"]) ;
      VNT = as<arma::mat>(fe_ad_inter["VNT"]) ;

      if (force==1||force==3) {
        alpha = as<arma::mat>(fe_ad_inter["alpha"]) ;
      }
      if (force==2||force==3) {
        xi = as<arma::mat>(fe_ad_inter["xi"]) ;
      }
      niter = as<int>(fe_ad_inter["niter"]) ;
    } 
    else {
      if (force==0) {
        U = YY ;
        fit.fill(mu) ;
      } else {
        // add fe; iteration
        List fe_ad = fe_ad_iter(YY, Y0, I, force, tol) ;
        mu = as<double>(fe_ad["mu"]) ;
        U = as<arma::mat>(fe_ad["e"]) ;
        fit = as<arma::mat>(fe_ad["fit"]) ;
        if (force==1||force==3) {
          alpha = as<arma::mat>(fe_ad["alpha"]) ;
        }
        if (force==2||force==3) {
          xi = as<arma::mat>(fe_ad["xi"]) ;
        }
        niter = as<int>(fe_ad["niter"]) ;
      }
    } 
  } 
  else {
    /* starting value:  the OLS estimator */
    invXX = XXinv(XX) ; // compute (X'X)^{-1}, outside beta iteration 
    if (r==0) {
      // add fe, covar; iteration
      List fe_ad = fe_ad_covar_iter(XX, invXX, 
                                    YY, Y0, I, beta0, force, tol) ;
      mu = as<double>(fe_ad["mu"]) ;
      beta = as<arma::mat>(fe_ad["beta"]) ;
      U = as<arma::mat>(fe_ad["e"]) ;
      fit = as<arma::mat>(fe_ad["fit"]) ;
      if (force==1||force==3) {
        alpha = as<arma::mat>(fe_ad["alpha"]) ;
      }
      if (force==2||force==3) {
        xi = as<arma::mat>(fe_ad["xi"]) ;
      }
      niter = as<int>(fe_ad["niter"]) ;
    } 
    else if (r > 0) {       
      // add, covar, interactive, iteration
      List fe_ad_inter_covar = fe_ad_inter_covar_iter(XX, invXX,
               YY, Y0, I, beta0, force, 0, r, 0, tol) ;
      mu = as<double>(fe_ad_inter_covar["mu"]) ;
      beta = as<arma::mat>(fe_ad_inter_covar["beta"]) ;
      U = as<arma::mat>(fe_ad_inter_covar["e"]) ;
      fit = as<arma::mat>(fe_ad_inter_covar["fit"]) ;

      factor = as<arma::mat>(fe_ad_inter_covar["factor"]) ;
      lambda = as<arma::mat>(fe_ad_inter_covar["lambda"]) ;
      VNT = as<arma::mat>(fe_ad_inter_covar["VNT"]) ;

      if (force==1||force==3) {
        alpha = as<arma::mat>(fe_ad_inter_covar["alpha"]) ;
      }
      if (force==2||force==3) {
        xi = as<arma::mat>(fe_ad_inter_covar["xi"]) ;
      }
      niter = as<int>(fe_ad_inter_covar["niter"]) ;
    }
  } 
    
  /* sigma2 and IC */
  // number of estimated parameters
  double np = r * (N + T) - pow(double(r),2) + p1 + 1 ;
  if (force == 1) {
    np = np + (N - 1) - r ;
  }
  else if (force == 2) {
    np = np + (T - 1) - r ;
  } 
  else if (force == 3) {
    np = np + (N - 1) + (T - 1) - 2 * r ;
  }

  sigma2 = trace(U * U.t())/ (obs - np) ;

  IC = log(sigma2) + np * log ( obs ) / ( obs ) ;

  // PC criterion in Li 2018
  // mean squared error
  double mse = trace(U * U.t()) / obs ;

  double m1 = 0 ;
  if (N < 60) {
    m1 = 60 - N ;
  }
  double m2 = 0 ;
  if (T < 60) {
    m2 = 60 - T ;
  }

  double C = (N + m1) * (T + m2) / obs ;
  double PC = mse + r * sigma2 * C * (N + T) / obs * log (double(obs)/double(N + T)) ;
    
  //-------------------------------#
  // Storage
  //-------------------------------# 

  List output ;

  if(p>0) {
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

  output["mu"] = mu ;   
  output["fit"] = fit ;  

  if ( !(force == 0 && r == 0 && p1 == 0) ) {
    output["niter"] = niter ;
  }
  if (force ==1 || force == 3) {
    output["alpha"] = alpha ;
  }
  if (force ==2 || force == 3) {
    output["xi"] = xi ;
  }
  if (r > 0) {
    output["factor"] = factor ;
    output["lambda"] = lambda ;
    output["VNT"] = VNT ;
  }
  output["residuals"] = U ;
  output["sigma2"] = sigma2 ;
  output["IC"] = IC ;
  output["PC"] = PC ;
  output["validX"] = validX ;
  return(output);
 
}


/* Interactive Fixed Effects: matrix completion */
// [[Rcpp::export]]
List inter_fe_mc (arma::mat Y,
                  arma::mat Y0,
                  arma::cube X,
                  arma::mat I,
                  arma::mat beta0,
                  int r, // r > 0, the outcome has a factor-type fixed effect; r = 0 else
                  double lambda,
                  int force,
                  double tol = 1e-5
                  ) {
  
  /* Dimensions */
  int T = Y.n_rows ;
  int N = Y.n_cols ;
  int p = X.n_slices ;
  double obs = accu(I) ;
  int niter = 0 ;
  int validF = 1 ;
  arma::mat beta ; 
  arma::mat U ;
  double mu = 0 ;
  double mu_Y = 0 ;

  arma::mat alpha(N, 1, arma::fill::zeros) ;
  arma::mat xi(T, 1, arma::fill::zeros) ;
  arma::mat fit(T, N, arma::fill::zeros) ;
  // double sigma2 = 0;


  arma::mat invXX ;

  /* duplicate data */
  arma::mat YY = Y;
  arma::cube XX = X;

  
  /* check if XX has enough variation */
  int p1 = p; 
  arma::mat X_invar(p, 1, arma::fill::zeros); // =1 if invar

  int j = 0;

  for(int i=0; i<p1; i++){
    if (arma::accu(abs(XX.slice(i))) < 1e-5) {
      XX.shed_slice(i);
      X_invar(j,0) = 1;
      i--;
      p1--;
    }
    j++;
  }

  int validX = 1 ;
  if(p1==0){
    validX = 0 ;
    if (force == 0 && r == 0) { // no covariate and force == 0 and r == 0 
      mu_Y = accu(YY)/obs ;
      mu = mu_Y ;
      YY = FE_adj(YY - mu_Y, I) ;
    }
  }

  /* Main Algorithm */ 
  if (p1 == 0) {
    if (r > 0) {
      // add fe ; inter fe ; iteration
      List fe_ad_inter = fe_ad_inter_iter(YY, Y0, I, force, 1, 1, lambda, tol) ;
      mu = as<double>(fe_ad_inter["mu"]) ;
      U = as<arma::mat>(fe_ad_inter["e"]) ;
      fit = as<arma::mat>(fe_ad_inter["fit"]) ;
      if (force==1||force==3) {
        alpha = as<arma::mat>(fe_ad_inter["alpha"]) ;
      }
      if (force==2||force==3) {
        xi = as<arma::mat>(fe_ad_inter["xi"]) ;
      }
      niter = as<int>(fe_ad_inter["niter"]) ;
      validF = as<int>(fe_ad_inter["validF"]) ;
    } 
    else {
      if (force==0) {
        U = YY ;
        fit.fill(mu) ;
        validF = 0 ;
      } else {
        // add fe; iteration
        List fe_ad = fe_ad_iter(YY, Y0, I, force, tol) ;
        mu = as<double>(fe_ad["mu"]) ;
        U = as<arma::mat>(fe_ad["e"]) ;
        fit = as<arma::mat>(fe_ad["fit"]) ;
        if (force==1||force==3) {
          alpha = as<arma::mat>(fe_ad["alpha"]) ;
        }
        if (force==2||force==3) {
          xi = as<arma::mat>(fe_ad["xi"]) ;
        }
        niter = as<int>(fe_ad["niter"]) ;
        validF = 0 ;
      }
    } 
  } 
  else {
    /* starting value:  the OLS estimator */
    invXX = XXinv(XX) ; // compute (X'X)^{-1}, outside beta iteration 
    if (r==0) {
      // add fe, covar; iteration
      List fe_ad = fe_ad_covar_iter(XX, invXX,
                                    YY, Y0, I, beta0, force, tol) ;
      mu = as<double>(fe_ad["mu"]) ;
      beta = as<arma::mat>(fe_ad["beta"]) ;
      U = as<arma::mat>(fe_ad["e"]) ;
      fit = as<arma::mat>(fe_ad["fit"]) ;
      if (force==1||force==3) {
        alpha = as<arma::mat>(fe_ad["alpha"]) ;
      }
      if (force==2||force==3) {
        xi = as<arma::mat>(fe_ad["xi"]) ;
      }
      niter = as<int>(fe_ad["niter"]) ;
      validF = 0 ;
    } 
    else if (r > 0) {       
      // add, covar, interactive, iteration
      List fe_ad_inter_covar = fe_ad_inter_covar_iter(XX, invXX,
               YY, Y0, I, beta0, 
               force, 1, 1, lambda, tol) ;
      mu = as<double>(fe_ad_inter_covar["mu"]) ;
      beta = as<arma::mat>(fe_ad_inter_covar["beta"]) ;
      U = as<arma::mat>(fe_ad_inter_covar["e"]) ;
      fit = as<arma::mat>(fe_ad_inter_covar["fit"]) ;
      if (force==1||force==3) {
        alpha = as<arma::mat>(fe_ad_inter_covar["alpha"]) ;
      }
      if (force==2||force==3) {
        xi = as<arma::mat>(fe_ad_inter_covar["xi"]) ;
      }
      niter = as<int>(fe_ad_inter_covar["niter"]) ;
      validF = as<int>(fe_ad_inter_covar["validF"]) ;
    }
  } 
    
  //-------------------------------#
  // Storage
  //-------------------------------# 

  List output ;

  if(p>0) {
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

  output["mu"] = mu ;   
  output["fit"] = fit ; 
  output["validF"] = validF ; 

  if ( !(force == 0 && r == 0 && p1 == 0) ) {
    output["niter"] = niter ;
  }
  if (force ==1 || force == 3) {
    output["alpha"] = alpha ;
  }
  if (force ==2 || force == 3) {
    output["xi"] = xi ;
  }
  output["residuals"] = U ;
  output["validX"] = validX;
  return(output);
 
}

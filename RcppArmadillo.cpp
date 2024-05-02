#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


double log_lik(arma::mat X, arma::vec y, double alpha, arma::vec beta){
  double end =y.n_elem-1;
  arma::vec indices =arma::linspace<arma::vec>(0,end,y.n_elem);
  indices.transform([&](double val) { unsigned int val1=val; double s=alpha+arma::as_scalar(X.row(val1) * beta);  return ( ((y(val1)-1.0)*s) - std::log(1+std::exp(-s))         );} );      
  return( arma::sum(indices) );
  
  //arma::uword n=y.n_elem;
  //arma::vec log_lik_all(n);
  //for(int i=0;i< n;i++){
    //double s=alpha+arma::as_scalar(X.row(i) * beta);
    //log_lik_all(i)=((y(i)-1.0)*s) - std::log(1+std::exp(-s)) ;
  //}
  
  //return(arma::sum(log_lik_all));
}

double log_lik_sel(arma::mat X, arma::vec y, double alpha, arma::vec beta, arma::vec gamma){
  double end =y.n_elem-1;
  arma::uvec zero_ind=arma::find(gamma==0.0);
  beta.shed_rows(zero_ind);
  X.shed_cols(zero_ind);
  arma::vec indices =arma::linspace<arma::vec>(0,end,y.n_elem);
  indices.transform([&](double val) { unsigned int val1=val; double s=alpha+arma::as_scalar(X.row(val1) * beta);  return ( ((y(val1)-1.0)*s) - std::log(1+std::exp(-s))         );} );      
  return( arma::sum(indices) );
  
  
}




// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat X, arma::vec y, double alpha_init, arma::vec beta_init, int M, arma::vec tuning_param, bool gibbs){
  int dim=1+beta_init.n_elem;
  double r;
  double u;
  arma::vec log_density(M,arma::fill::zeros);
  arma::mat diagonal(dim-1, dim-1, arma::fill::eye);
  arma::mat samples(dim,M,arma::fill::zeros);
  arma::mat acceptance(1,tuning_param.n_elem,arma::fill::zeros);
  double alpha_old=alpha_init;
  arma::vec beta_old =beta_init;
  double log_lik_old=log_lik(X,y,alpha_old,beta_old);
  samples.col(0)=arma::join_cols(alpha_old*arma::ones(1),beta_old);
  log_density(0)=log_lik_old;
  
  double alpha_current;
  double log_lik_current;
  arma::vec beta_current;
  
  double beta_prime;
  
  for(int t=1;t<M;t++){
    alpha_current=arma::randn<double>(arma::distr_param(alpha_old,tuning_param(0)));
    log_lik_current=log_lik(X,y,alpha_current,beta_old);
    r=std::min(1.0,std::exp( log_lik_current-log_lik_old  ));
    u=arma::randu();
    if(u<=r){
      alpha_old=alpha_current;
      log_lik_old=log_lik_current;
      acceptance[0]=acceptance[0]+1;
    }
  
  if(gibbs==0){
  beta_current=arma::mvnrnd(beta_old,tuning_param(1)*diagonal);
  
  log_lik_current=log_lik(X,y,alpha_old,beta_current);
  r=std::min(1.0,std::exp( log_lik_current-log_lik_old  ));
  u=arma::randu();
  if(u<=r){
    beta_old=beta_current;
    log_lik_old=log_lik_current;
    acceptance[1]=acceptance[1]+1;
  }
  } else{
    for(int j=0;j<dim-1;j++){
     beta_prime=arma::randn<double>(arma::distr_param(beta_old(j),tuning_param(1+j)));
     beta_current=beta_old;
     beta_current(j)=beta_prime;
     log_lik_current=log_lik(X,y,alpha_old,beta_current);
     r=std::min(1.0,std::exp( log_lik_current-log_lik_old  ));
     u=arma::randu();
     if(u<=r){
       beta_old=beta_current;
       log_lik_old=log_lik_current;
       acceptance[1+j]=acceptance[1+j]+1;
     }
    }
  }
  samples.col(t)=arma::join_cols(alpha_old*arma::ones(1),beta_old);
  log_density(t)=log_lik_old;
  }
  
 arma::mat t_samples=samples.t();
 Rcpp::List L(3);
  L[0]=t_samples;
  L[1]=acceptance/M;
  L[2]=log_density;
  return L; 
  
}

// [[Rcpp::export]]
arma::mat prediction(arma::mat X_new, arma::mat mcmc_samples){
  mcmc_samples=mcmc_samples.t();
  arma::uword p=X_new.n_cols;
  int n1=X_new.n_rows;
  double end =n1-1;
  arma::mat pred(n1,mcmc_samples.n_cols,arma::fill::zeros);
  for(int i=0;i<mcmc_samples.n_cols;i++){
    arma::vec indices =arma::linspace<arma::vec>(0,end,n1);
    indices.transform([&](double val) { unsigned int val1=val; double s=mcmc_samples.col(i)(0)+arma::as_scalar(X_new.row(val1) * mcmc_samples.col(i).subvec(1,p));  return ( 1.0/(1+std::exp(-s))         );} ); 
    pred.col(i)=indices;
  }
 
 return(pred.t()); 
}



// [[Rcpp::export]]
Rcpp::List run_mcmc_sel(arma::mat X, arma::vec y, double alpha_init, arma::vec beta_init, arma::vec gamma_init, int M, arma::vec tuning_param){
  int dim_var=beta_init.n_elem;
  int dim=1+2*dim_var;
  
  double r;
  double u;
  arma::uword k;
  arma::vec log_density(M,arma::fill::zeros);
  arma::mat samples(dim,M,arma::fill::zeros);
  arma::mat acceptance(1,tuning_param.n_elem,arma::fill::zeros);
  double alpha_old=alpha_init;
  arma::vec beta_old =beta_init;
  arma::vec gamma_old=gamma_init;
  double log_lik_old=log_lik_sel(X,y,alpha_old,beta_old,gamma_old);
  samples.col(0)=arma::join_cols(alpha_old*arma::ones(1),beta_old,gamma_old);
  log_density(0)=log_lik_old;
  
  double alpha_current;
  double log_lik_current;
  arma::vec beta_current;
  arma::vec gamma_current;
  
  double beta_prime;
  
  for(int t=1;t<M;t++){
    alpha_current=arma::randn<double>(arma::distr_param(alpha_old,tuning_param(0)));
    log_lik_current=log_lik_sel(X,y,alpha_current,beta_old,gamma_old);
    r=std::min(1.0,std::exp(arma::log_normpdf(alpha_current, 0.0, 100.0)-arma::log_normpdf(alpha_old, 0.0, 100.0) +log_lik_current-log_lik_old  ));
    u=arma::randu();
    if(u<=r){
      alpha_old=alpha_current;
      log_lik_old=log_lik_current;
      acceptance[0]=acceptance[0]+1;
    }
    
      k=arma::randi<arma::uword>( arma::distr_param(0,dim_var-1) );
      if(gamma_old(k)==0.0){
       gamma_current=gamma_old;
       gamma_current(k)=1.0;
       beta_current=beta_old;
       beta_current(k)=arma::randn<double>(arma::distr_param(beta_old(k),tuning_param(1)));
       log_lik_current=log_lik_sel(X,y,alpha_old,beta_current,gamma_current);
       r=std::min(1.0,std::exp(arma::log_normpdf(beta_current(k), 0.0, 100.0)+ log_lik_current-log_lik_old  ));
      
        
      } else{
        gamma_current=gamma_old;
        gamma_current(k)=0.0;
        beta_current=beta_old;
        beta_current(k)=0.0;
        log_lik_current=log_lik_sel(X,y,alpha_old,beta_current,gamma_current);
        r=std::min(1.0,std::exp(log_lik_current-log_lik_old-arma::log_normpdf(beta_old(k), 0.0, 100.0)  ));
      }
      
      u=arma::randu();
      if(u<=r){
        beta_old=beta_current;
        gamma_old=gamma_current;
        log_lik_old=log_lik_current;
        acceptance[1]=acceptance[1]+1;
      }
     for(int j=0;j<dim_var;j++){
       if(gamma_old(j)==1.0){
         beta_current=beta_old;
         beta_current(j)=arma::randn<double>(arma::distr_param(beta_old(j),tuning_param(1)));
         log_lik_current=log_lik_sel(X,y,alpha_old,beta_current,gamma_old);
         r=std::min(1.0,std::exp(arma::log_normpdf(beta_current(j), 0.0, 100.0)-arma::log_normpdf(beta_old(j), 0.0, 100.0)+ log_lik_current-log_lik_old  ));
         u=arma::randu();
         if(u<=r){
           beta_old=beta_current;
           
           log_lik_old=log_lik_current;
           //acceptance[1]=acceptance[1]+1;
         }
       }
     }
     
     
    samples.col(t)=arma::join_cols(alpha_old*arma::ones(1),beta_old,gamma_old);
    log_density(t)=log_lik_old;
  }
  
  arma::mat t_samples=samples.t();
  Rcpp::List L(3);
  L[0]=t_samples;
  L[1]=acceptance/M;
  L[2]=log_density;
  return L; 
  
}



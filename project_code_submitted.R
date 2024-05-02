
# source the c++ file 
Rcpp::sourceCpp("RcppArmadillo.cpp")

#Load Data
data=read.csv("https://raw.githubusercontent.com/jbrownlee/Datasets/master/pima-indians-diabetes.csv")

#data=read.csv("F:/Biostat/Bayesian/diabetes.csv")

# Define Predictor Matrix
X=as.matrix(data[,-9])
X=scale(X)

#Binary Outcome
y=as.vector(data[,9])

n=nrow(data) # sample size
p=ncol(X)   #no. of predictors


# Run the MCMC algorithm
M=50000   #No of MCMC samples
alpha_old=0  #initial values of alpha
beta_old=rep(0,p)  #initial value of betas
tuning=c(0.25,rep(0.25,p)) #tuning params of proposal

#mcmc samples
samples=run_mcmc(X,y,alpha_old,beta_old,M,tuning,gibbs=T)

#Posterior estimates of the coefficients
round(apply(samples[[1]][(burn_in+1):M,],2,mean),3)[1:9]

#GLM estimates of the coefficients
glm(y~X,family=binomial(link="logit"), control=list('maxit'=100))$coefficients

quant=function(x){
  quantile(x,c(0.025,0.975))
}

#get credible intervals
apply(samples[[1]][(burn_in+1):M,],2,quant)

#log-density trace plot
plot(1:M,samples[[3]],type='l',col='black', xlab="Index",ylab='log_density')
abline(v=burn_in,col="green")

# Define a function to run three parallel chains
run_mcmc_i=function(i){
  set.seed(100+i)
  alpha_old=rnorm(1)
  beta_old=rnorm(p)
  run_mcmc(X,y,alpha_old,beta_old,M,tuning, gibbs=T)
}


library(parallel)

cl <- parallel::makeCluster(4)

parallel::clusterEvalQ(cl, {library(Rcpp);library(RcppArmadillo);Rcpp::sourceCpp("RcppArmadillo.cpp")})
clusterExport(cl,list('X','y','M','p','tuning'))

#t1=Sys.time()
chains=parLapply(cl,1:3,run_mcmc_i) # save the three chains
#t2=Sys.time()
stopCluster(cl)

#MCMC draws for beta8 for the three chains---superimposing trace plots
var=9
plot(1:M,chains[[1]][[1]][,var],type='l',col='black', xlab="Index",ylab='MCMC Draws', main="Draws for beta_8")
lines(1:M,chains[[2]][[1]][,var],type='l',col='green')
lines(1:M,chains[[3]][[1]][,var],type='l',col='red')

burn_in=5000

#Gelman Rubin---for computing scale reduction factors
library(coda)
l=list(mcmc(chains[[1]][[1]][(burn_in+1):M,]), mcmc(chains[[2]][[1]][(burn_in+1):M,]   ), mcmc(chains[[3]][[1]][(burn_in+1):M,]   )  )
l=mcmc.list(l)
diagnostics=gelman.diag(l)
diagnostics


# Prediction 

n_train=600
n_test=n-n_train
set.seed(100+1)
train_rows=sample(1:n, size=n_train, replace=F)

X_train=X[train_rows,] #train set
y_train=y[train_rows]
test_rows=setdiff(1:n,train_rows)
X_test=X[test_rows,] # test set
y_test=y[test_rows]

set.seed(100+1)
alpha_old=rnorm(1)
beta_old=rnorm(p)
tuning=c(0.25,rep(0.25,p))
samples=run_mcmc(X_train,y_train,alpha_old,beta_old,M,tuning, gibbs=T)

pred=prediction(X_test,samples[[1]])

# Draws of posterior probability for a particular test observation
obs=94
plot(1:M,pred[,obs],type='l',col='black',xlab='Indices',ylab="Draws of p")
abline(v=burn_in,col="green")


mean_prob=apply(pred,2,mean)

# predicted lables
y_pred=ifelse(mean_prob>=0.5,1,0)

#confusion matrix
table(y_test,y_pred)


# Bayesian Variable Selection

M=50000

tuning=c(0.25,0.1)
run_mcmc_i=function(i){
  set.seed(100+i)
  alpha_old=rnorm(1)
  beta_old=rnorm(p)
  gamma_old=sample(c(0,1),size=p,replace=T)
  beta_old[which(gamma_old==0)]=0
  run_mcmc_sel(X,y,alpha_old,beta_old,gamma_old,M,tuning)
  
}


library(parallel)

cl <- parallel::makeCluster(4)

parallel::clusterEvalQ(cl, {library(Rcpp);library(RcppArmadillo);Rcpp::sourceCpp("RcppArmadillo.cpp")})
clusterExport(cl,list('X','y','M','p','tuning'))


chains=parLapply(cl,1:3,run_mcmc_i)
stopCluster(cl)


burn_in=5000
library(coda)
l=list(mcmc(chains[[1]][[1]][(burn_in+1):M,1:9]), mcmc(chains[[2]][[1]][(burn_in+1):M,1:9]   ), mcmc(chains[[3]][[1]][(burn_in+1):M,1:9]   )  )
l=mcmc.list(l)
diagnostics=gelman.diag(l)
diagnostics

#get posterior estimates of the coefficients
round(apply(chains[[1]][[1]][(burn_in+1):M,],2,mean),3)[2:9]

#get posterior probability of inclusion
round(apply(chains[[1]][[1]][(burn_in+1):M,],2,mean),3)[10:17]

#trace plot of log-likelihood
plot((burn_in+1):M,chains[[1]][[3]][(burn_in+1):M],type='l',col='black', xlab="Index",ylab='log_likelihood', main="log_likelihood")
lines((burn_in+1):M,chains[[2]][[3]][(burn_in+1):M],type='l',col='green')
lines((burn_in+1):M,chains[[3]][[3]][(burn_in+1):M],type='l',col='red')

#prediction under variable selection setup
n_train=600
n_test=n-n_train
set.seed(100+1)
train_rows=sample(1:n, size=n_train, replace=F)

X_train=X[train_rows,]
y_train=y[train_rows]
test_rows=setdiff(1:n,train_rows)
X_test=X[test_rows,]
y_test=y[test_rows]

set.seed(100+1)
alpha_old=rnorm(1)
beta_old=rnorm(p)
gamma_old=sample(c(0,1),size=p,replace=T)
beta_old[which(gamma_old==0)]=0
tuning=c(0.25,0.5)


burn_in=5000

samples=run_mcmc_sel(X_train,y_train,alpha_old,beta_old,gamma_old,M,tuning)

pred=prediction(X_test,samples[[1]])


mean_prob=apply(pred[(1+burn_in):M,],2,mean)

y_pred=ifelse(mean_prob>=0.5,1,0)

table(y_test,y_pred)


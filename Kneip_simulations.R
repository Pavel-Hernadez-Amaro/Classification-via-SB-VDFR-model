# Packages, Sources, Data, ... --------------------------------------------
library(ReconstPoFD)
library(pROC)
library(Epi)
library(MASS)
library(fpca)
library(fda)  
library(mvtnorm)
library(refund)
library(foreach)
library(doParallel)
#setwd("/afs/sta.uniroma1.it/user/m/mstefanu")      # <--- change this
#setwd("C:/Users/Marco/Desktop/Dottorato - III anno/Paper - PCA--based discrimination of partially observed functional data with an application to Aneurisk65 dataset/Codes/Simulations/Functions")
source("fpca-modifica.R")
source("myfunctions.R")
source("myfunctions-cv.R")
#source("C:/Users/Marco/Desktop/Dottorato - II anno/ANEURISK-final/R/Codes/simulation XIII - cross validation/myfunctions-cv.R")
#setwd("C:/Users/Marco/Desktop/Dottorato - II anno/ANEURISK-final/R/Codes")
options(error=NULL)

# no_cores <- 4
# registerDoParallel(no_cores)

# parameters --------------------------------------------------------------

# SETTINGS THAT CAME WITH THE CODE

set.seed(1327765)
S      <- 50
maxi   <- 5
p      <- 150
n1     <- 50
n2     <- 50
N=n      <- n1+n2
maxl   <- 9 # 9
k      <- seq(10, 40, by=5)
lambda.sm   <- 100
lambda.buja <- 10^seq(-2, 4, length.out = 12)
grid   <- seq(0, 1, length=p)
err1   <- 0.1
err2   <- 0.1
unif.miss=TRUE

nknots <- 18
par.knots <- nknots/3

knots =  c(sort(runif(par.knots, grid[1], grid[p]/3)),
           sort(runif(par.knots, grid[p]/3, 2*grid[p]/3)),
           sort(runif(par.knots, 2*grid[p]/3, grid[p])))

nbasis <- length(knots) +2

# mean1  <- c(0, 0, 0, 0, 1, 2, 1, 0,-1, 2, 2,-1, 0, 0.5, 1, 0.5, 0, 0, 0, 0)
# mean2  <- rev(mean1)
mean1  <- rnorm(nbasis)
mean2  <- rev(mean1)
coef   <- rep(0.3,nbasis)

c_kneip=18
R=50


X_list=t_list=list()

AUC=res_clas_12=res_clas_34=prob_pred=optimal_cut_12=optimal_cut_34=err_opt_12=err_opt_34=array(dim = N)
miss_classification_error=array(dim = R)
prob_pred_Kneip=res_clas_final=optimal_cut=array(dim = c(R,N))


for (j in 1:R) {
  
print(c("j =", j))

x1       <- data.gen.splines.uneven(n1, p, grid, knots, mean1, err1, coef)
x2       <- data.gen.splines.uneven(n2, p, grid, knots, mean2, err2, coef)
xbind    <- rbind(x1$data, x2$data)
x.smbind <- rbind(x1$smooth.data, x2$smooth.data)
y        <- c(rep(0,n1), rep(1,n2))  
# nu=abs(rowSums(x$x.miss,na.rm=TRUE)/100) # POSSIBLE NU THAT IS DATA RELATED (BAD RESULTS IN 5 ITERATIONS)
# y=rbinom(n1+n2,1,nu) 
if(unif.miss==TRUE) {
  
  x=data.agg.miss(xbind)
  
}else {
  print("Beta")
  x=data.agg.miss.betadecay(xbind,1) 
}

############## Reconstruction of curves

for (i in 1:dim(x$x.miss)[1]) {
  
  X_list[[i]]=x$x.miss[i,x$domain[i,1]:x$domain[i,2]]
  t_list[[i]]=grid[x$domain[i,1]:x$domain[i,2]]
}

reconst_result <- reconstructKneipLiebl(Ly = X_list, Lu = t_list, method = 'Error>0_AlignYES')

X_reconst_mat  <- t(matrix(unlist(reconst_result[['Y_reconst_list']]), ncol=100))

for (i in 1:length(y)) {
  
  print(c("i =", i))
  
  X_aux=X_reconst_mat[-i,]
  y_aux=y[-i]
  
  fit_Kneip <- pfr( y_aux~ lf(X_aux, argvals = grid, bs="ps",k=18),family=binomial())
  
  aux=ROC(fit_Kneip$fitted.values,y_aux,plot = NULL)
  
  best_pos_34=which.max(aux$res[,3]+aux$res[,4]) # which.max(aux$res[,1]+aux$res[,2])
  optimal_cut_34[i]=aux$res[best_pos_34,5]
  
  best_pos_12=which.max(aux$res[,1]+aux$res[,2]) # which.max(aux$res[,1]+aux$res[,2])
  optimal_cut_12[i]=aux$res[best_pos_12,5]
  
  AUC[i]=1-aux$AUC
  
  prob_pred[i]=predict(fit_Kneip, newdata = list(X_aux=t(as.matrix(X_reconst_mat[i,]))),type='response')
  
  if (prob_pred[i]<optimal_cut_34[i]) {
    
    res_clas_34[i]=0
  }else{
    res_clas_34[i]=1
  }
  
  err_opt_34[i]  = as.double(res_clas_34[i]!=y[i])
  
  if (prob_pred[i]<optimal_cut_12[i]) {
    
    res_clas_12[i]=0
  }else{
    res_clas_12[i]=1
  }
  
  err_opt_12[i]  = as.double(res_clas_12[i]!=y[i])
  
  
}

aux=c(sum(err_opt_12),sum(err_opt_34))
aux_2=cbind(res_clas_12,res_clas_34)

min_track=which.min(aux)

aux_3=cbind(optimal_cut_12,optimal_cut_34)

res_clas_final[j,]=aux_2[,min_track]
miss_classification_error[j]=aux[min_track]
optimal_cut[j,]=aux_3[,min_track]
prob_pred_Kneip[j,]=prob_pred

}




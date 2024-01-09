# STUDY OF SIMULATION USING A LOGISTIC REGRESSION MODEL TO SIMULATE THE CURVES AND THE RESPONSE VARIABLE (USING A FUNCTIONAL COEFFICIENT)
#
# NAMED SIMULATION_2

# Packages, Sources, Data, ... --------------------------------------------
library(pROC)
library(Epi)
library(MASS)
library(fpca)
library(fda)  
library(mvtnorm)
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

# parameters --------------------------------------------------------------

# set.seed(1327765)
# S      <- 50
# maxi   <- 5
# p      <- 150
# n1     <- 50
# n2     <- 50
# n      <- n1+n2
# maxl   <- 9
# k      <- seq(10, 40, by=5)
# lambda.sm   <- 100
# lambda.buja <- 10^seq(-2, 4, length.out = 12)
# grid   <- seq(0, 1, length=p)
# err1   <- .4
# err2   <- .4
# 
# 
# unif.miss=T
# 
# nknots <- 27
# par.knots <- nknots/3
# knots = c(sort(runif(par.knots, grid[1], grid[p]/3)),
#           sort(runif(par.knots, grid[p]/3, 2*grid[p]/3)),
#           sort(runif(par.knots, 2*grid[p]/3, grid[p])))
# nbasis <- length(knots) +2
# 
# mean1  <- rnorm(nbasis)
# mean2  <- mean1 + .1*c(1:10,5,rep(0,7),5, 10:1)
# coef <- .2*c( sort(seq(1,30, length.out = 10), decreasing = T) ,rep(5,9), seq(1,30, length.out = 10))


set.seed(1327765)
S      <- 100
maxi   <- 5
N=100
p      <- 150
maxl   <- 9
k      <- seq(10, 40, by=5)
# lambda.sm   <- 100
lambda.buja <- 10^seq(-2, 4, length.out = 12)
case=1

### MY ADDITIONS

probs=array(dim=c(N,maxl+1,S)) # probabilities of the logistic model

Beta_estimated=array(dim=c(p,maxl+1,S)) # ESTIMATED FUNCTIONAL COEFFICIENT FOR THE SB_VDFR MODEL

err_t=err_james=err_yao=err_buja=err_t_lda=err_james_lda=err_yao_lda=err_buja_lda=err_SB=array(dim=c(S,maxl+1)) # classifiactions errors

eigenfun_t=eigenfun_james=eigenfun_yao=eigenfun_buja=eigenfun_t_lda=eigenfun_james_lda=eigenfun_yao_lda=eigenfun_buja_lda=array(dim = c(p,maxi,maxl+1,S)) # eigenfunctions of all method

scores_t=scores_james=scores_yao=scores_buja=scores_t_lda=scores_james_lda=scores_yao_lda=scores_buja_lda=array(dim = c(N,maxi,maxl+1,S)) # scores of all methods

best_dim_t=best_dim_james=best_dim_yao=best_dim_buja=best_dim_t_lda=best_dim_james_lda=best_dim_yao_lda=best_dim_buja_lda=array(dim=c(S,maxl+1)) # optimal number of scores 
#(the one that had achieved the lowest classification error)

optimal_cut=array(dim=c(S,maxl+1))

domain=array(dim = c(N,2))

# c1=25 # number of basis for the SB_VDFR model
# c2=25

# computations ------------------------------------------------------------

# options(error=recover)
# no_cores <- detectCores() - 4
# registerDoParallel(no_cores)
# 
# 
# # PARALLEL / PLAIN
# ao <- foreach(s = 1:S, .packages=c("MASS", "fda", "fpca", "mvtnorm"), .export="myfpca.mle", .combine = rbind)  %dopar% 
#   all.together.splines.last(knots, mean1, err1, mean2, err2, n1, n2, p, grid, c1, c2, maxi,
#                             k, lambda.sm, lambda.buja, coef, unif.miss=T,
#                             rate=NA)

for (s in 1:S) {
  
  set.seed(1000+s)
  
  print(c("s=", s))
  
  ao=all.together.splines.last.totalcv_sim_2(N, p, case, maxi, maxl, k, lambda.buja,unif.miss=T, rate=NA)
  
  # ao_SB_CV=all.together.splines.last_solo_SB_CV(knots, mean1, err1, mean2, err2, n1, n2, p, grid, maxi, maxl,
  #                           k, lambda.sm, lambda.buja, coef, unif.miss=T,
  #                           rate=NA)
  
  
  err_t[s,]=ao$err.cv.t
  err_james[s,]=ao$err.cv.james
  err_yao[s,]=ao$err.cv.yao
  err_buja[s,]=ao$err.cv.buja
  # 
  err_SB[s,]=ao$err_SB_VDFR
  probs[,,s]=ao$prob_SB_VDFR
  optimal_cut[s,]=ao$optimal_cut_SB_VDFR
  Beta_estimated[,,s]=ao$Beta_estimated
  # 
  eigenfun_t[,,,s]=ao$eigenfun.t
  eigenfun_james[,,,s]=ao$eigenfun.james
  eigenfun_yao[,,,s]=ao$eigenfun.yao
  eigenfun_buja[,,,s]=ao$eigenfun.buja
  
  scores_t[,,,s]=ao$scores.t
  scores_james[,,,s]=ao$scores.james
  scores_yao[,,,s]=ao$scores.yao
  scores_buja[,,,s]=ao$scores.buja
  
  best_dim_t[s,]=ao$best.t
  best_dim_james[s,]=ao$best.james
  best_dim_yao[s,]=ao$best.yao
  best_dim_buja[s,]=ao$best.buja
  
  
  ########
  
  err_t_lda[s,]=ao$err.cv.t_lda
  err_james_lda[s,]=ao$err.cv.james_lda
  err_yao_lda[s,]=ao$err.cv.yao_lda
  err_buja_lda[s,]=ao$err.cv.buja_lda
  
  eigenfun_t_lda[,,,s]=ao$eigenfun.t_lda
  eigenfun_james_lda[,,,s]=ao$eigenfun.james_lda
  eigenfun_yao_lda[,,,s]=ao$eigenfun.yao_lda
  eigenfun_buja_lda[,,,s]=ao$eigenfun.buja_lda
  
  scores_t_lda[,,,s]=ao$scores.t_lda
  scores_james_lda[,,,s]=ao$scores.james_lda
  scores_yao_lda[,,,s]=ao$scores.yao_lda
  scores_buja_lda[,,,s]=ao$scores.buja_lda
  
  best_dim_t_lda[s,]=ao$best.t_lda
  best_dim_james_lda[s,]=ao$best.james_lda
  best_dim_yao_lda[s,]=ao$best.yao_lda
  best_dim_buja_lda[s,]=ao$best.buja_lda
  
  
}

### THESE STATISTICS SUMMARIZE ALL THE DOMAIN GAPS IN ONE NUMBER AND DO NOT MAKE SUMMARIES FOR EVERY DOMAIN   

err_t_min=apply(err_t, 1, min)
err_james_min=apply(err_james, 1, min)
err_yao_min=apply(err_yao, 1, min)
err_buja_min=apply(err_buja, 1, min)
err_SB_min=apply(err_SB, 1, min)

err_t_min_lda=apply(err_t_lda, 1, min)
err_james_min_lda=apply(err_james_lda, 1, min)
err_yao_min_lda=apply(err_yao_lda, 1, min)
err_buja_min_lda=apply(err_buja_lda, 1, min)

where_err_t_min=apply(err_t, 1, which.min)
where_err_james_min=apply(err_james, 1, which.min)
where_err_yao_min=apply(err_yao, 1, which.min)
where_err_buja_min=apply(err_buja, 1, which.min)
where_err_SB_min=apply(err_SB, 1, which.min)

where_err_t_min_lda=apply(err_t_lda, 1, which.min)
where_err_james_min_lda=apply(err_james_lda, 1, which.min)
where_err_yao_min_lda=apply(err_yao_lda, 1, which.min)
where_err_buja_min_lda=apply(err_buja_lda, 1, which.min)

## these graphs tell where is the best domain and how many times is it the best 
barplot(table(where_err_t_min))
barplot(table(where_err_james_min))
barplot(table(where_err_yao_min))
barplot(table(where_err_buja_min))
barplot(table(where_err_SB_min))
barplot(table(where_err_t_min_lda))
barplot(table(where_err_james_min_lda))
barplot(table(where_err_yao_min_lda))
barplot(table(where_err_buja_min_lda))
#

where_class_error=rbind(where_err_t_min,where_err_james_min,where_err_yao_min,where_err_buja_min,where_err_SB_min)

class_error=rbind(err_t_min,err_james_min,err_yao_min,err_buja_min,err_SB_min)

rowSums(class_error)/S

boxplot(t(class_error))

#

where_class_error_lda=rbind(where_err_t_min_lda,where_err_james_min_lda,where_err_yao_min_lda,where_err_buja_min_lda,where_err_SB_min)

class_error_lda=rbind(err_t_min_lda,err_james_min_lda,err_yao_min_lda,err_buja_min_lda,err_SB_min)

rowSums(class_error_lda)/S

boxplot(t(class_error_lda))

# using the mean error in every domain gap

err_t_min_domains_median=apply(err_t, 2, median)
err_james_min_domains_median=apply(err_james, 2, median)
err_yao_min_domains_median=apply(err_yao, 2, median)
err_buja_min_domains_median=apply(err_buja, 2, median)
err_SB_min_domains_median=apply(err_SB, 2, median)

err_t_min_domains=apply(err_t, 2, mean)
err_james_min_domains=apply(err_james, 2, mean)
err_yao_min_domains=apply(err_yao, 2, mean)
err_buja_min_domains=apply(err_buja, 2, mean)
err_SB_min_domains=apply(err_SB, 2, mean)

err_t_min_lda_domains=apply(err_t_lda, 2, mean)
err_james_min_lda_domains=apply(err_james_lda, 2, mean)
err_yao_min_lda_domains=apply(err_yao_lda, 2, mean)
err_buja_min_lda_domains=apply(err_buja_lda, 2, mean)

err_t_min_lda_domains_median=apply(err_t_lda, 2, median)
err_james_min_lda_domains_median=apply(err_james_lda, 2, median)
err_yao_min_lda_domains_median=apply(err_yao_lda, 2, median)
err_buja_min_lda_domains_median=apply(err_buja_lda, 2, median)

class_error_domains=rbind(err_t_min_domains,err_james_min_domains,err_yao_min_domains,err_buja_min_domains,err_SB_min_domains)

class_error_domains_median=rbind(err_t_min_domains_median,err_james_min_domains_median,err_yao_min_domains_median,err_buja_min_domains_median,err_SB_min_domains_median)

rowSums(class_error_domains)/10

class_error_domains_lda=rbind(err_t_min_lda_domains,err_james_min_lda_domains,err_yao_min_lda_domains,err_buja_min_lda_domains,err_SB_min_domains)

# rowSums(class_error_domains_lda)/10

indx=10 # 1:10

boxplot(err_t[,indx],err_james[,indx],err_yao[,indx],err_buja[,indx],err_SB[,indx],xlab="")
axis(1, at=c(1,2,3,4,5),gap.axis = 0.2, labels=c("Without-Missing ","James","Yao","Buja","SB_VDFR")) #c("Without Missing", "James et al.", "Yao et al.", "Huang et al.", "SB_VDFR"))
legend("top",c("indx=",indx),ncol = 2)

boxplot(err_t_lda[,indx],err_james_lda[,indx],err_yao_lda[,indx],err_buja_lda[,indx],err_SB[,indx],xlab="")
axis(1, at=c(1,2,3,4,5),gap.axis = 0.75, labels=c("Without Missing","James","Yao","Buja","SB_VDFR")) #c("Without Missing", "James et al.", "Yao et al.", "Huang et al.", "SB_VDFR"))
legend("top",c("indx=",indx),ncol = 2)

#

########## experiments

# l=1
# t.par  = round(p/3+1-(l+1)*(p/3)/(maxl+1)):round(2*p/3+(l+1)*(p/3)/(maxl+1))
# 

indx_aux=39

curve=c(rep(NA,(domain[indx_aux][1]-1)),SB_VDFR_Class$x_h[[indx_aux]],rep(NA,dim(x)[2]-domain[indx_aux,][2]))

plot(x[indx_aux,],type="l")
lines(curve,col=2,lwd=2)
lines(curve,col=3,lwd=2)
legend("top",c("5 basis", "30 basis"),col=c(2,3),lwd=2)


curve_h=Beta_estimated[,10,50]

plot(c(1:150),curve_h)
lines(Beta_estimated[,10,50],col=2,lwd=2)

plot(c(1:150),Beta[1,])


# 
# 
# aux=roc(y,probs[,1])
# 
# is_naive_bad=0
# where=NULL
# err_aux=rep(0,S)
# 
# for (sim in 1:S) {
# 
#   aux=ROC(probs[,sim],y,plot = NULL)
#   
#   best_pos=which.max(aux$res[,1]+aux$res[,2])
#   
#   best_cut=aux$res[best_pos,5]
#   
#   res_clas_aux=rep(1,length(y))
#   res_clas_aux[which(probs[,sim]<best_cut)]=0
#   
#   err_aux[sim]=sum(res_clas_aux!=y)
#   
#   aux=err_SB[sim]<=err_aux[sim]
#   
#   if (aux==FALSE) {
#     
#     where=c(where,sim)
#     
#   }
#   
#   is_naive_bad=sum(is_naive_bad,aux)
#   
# }
# 
# res_clas_aux[where]
# 
# mix_err_SB=err_SB
# 
# mix_err_SB[where]=err_aux[where]
# 
# class_error_mix=rbind(class_error,mix_err_SB)
# 
# rowSums(class_error_mix)/S
# 
# boxplot(t(class_error_mix))
# 
# 
# # saving ------------------------------------------------------------------
# 
# #save.image("xxx.RData")





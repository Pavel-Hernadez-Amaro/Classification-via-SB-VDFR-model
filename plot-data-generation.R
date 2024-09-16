# Prelim ------------------------------------------------------------------

#setwd("C:/Users/Marco/Desktop/Dottorato - III anno/Paper - PCA--based discrimination of partially observed functional data with an application to Aneurisk65 dataset/Codes/Simulations/Functions")
source("myfunctions.R")


# Parameters --------------------------------------------------------------

# set.seed(1327765)
# n1    <- 20 # number of curves of type 1
# n2    <- 20 # number of curves of type 2
# p     <- 150
# grid  <- seq(0, 1, length.out = p)
# nknots <- 27
# par.knots <- nknots/3
# knots  <- c(sort(runif(par.knots, grid[1], grid[p]/3)),
#            sort(runif(par.knots, grid[p]/3, 2*grid[p]/3)),
#            sort(runif(par.knots, 2*grid[p]/3, grid[p])))
# nbasis <- length(knots) +2
# mean1  <- rnorm(nbasis)
# mean2  <- mean1 + .1*c(1:10,5,rep(0,7),5, 10:1)
# err1   <- .4
# err2   <- .4
# coef   <- .2*c( sort(seq(1,30, length.out = 10), decreasing = T) ,rep(5,9), seq(1,30, length.out = 10))
# unif.miss <- T #this decide the distributions of the Na's according to an uniform distribution (True) or Beta (False)

# CREATING THE SAME DATA AS THE PAPER

set.seed(1327765)
n1    <- 50 # number of curves of type 1
n2    <- 50 # number of curves of type 2
p     <- 150
grid  <- seq(0, 1, length.out = p)
nknots <- 18
par.knots <- nknots/3
knots  <- c(sort(runif(par.knots, grid[1], grid[p]/3)),
            sort(runif(par.knots, grid[p]/3, 2*grid[p]/3)),
            sort(runif(par.knots, 2*grid[p]/3, grid[p])))
nbasis <- length(knots) +2
mean1  <- rnorm(nbasis)#c(0, 0, 0, 0, 1, 2, 1, 0,-1, 2, 2,-1, 0, 0.5, 1, 0.5, 0, 0, 0, 0)
mean2  <- rev(mean1)
err1   <- .1
err2   <- .1
coef   <- rep(0.6,nbasis)
unif.miss <- TRUE #this decide the distributions of the Na's according to an uniform distribution (True) or Beta (False)

##
# plot(mean1, type="l",lwd=2)
# lines(mean2, type="l",col=2,lwd=2)

# Data --------------------------------------------------------------------

x1       <- data.gen.splines.uneven(n1, p, grid, knots, mean1, err1, coef)
x2       <- data.gen.splines.uneven(n2, p, grid, knots, mean2, err2, coef)
xbind    <- rbind(x1$data, x2$data)
x.smbind <- rbind(x1$smooth.data, x2$smooth.data)
y        <- c(rep(0,n1), rep(1,n2))  
x        <- if(unif.miss==TRUE) data.agg.miss(xbind) else data.agg.miss.betadecay(xbind) 


# Take a look
names(x)
dim(x$x.miss)
# matplot(t(x$x.miss), type= "l", lty = 1, col = gray(.8))
matplot(t(x$x.miss), type= "l", lty = 1,col = y+1)#c("#ff9900", "#0868ac") )

# matplot(x[1:50,], type= "l", lty = 1, col = y)
# matplot(x[51:100,], type= "l", lty = 1, col = y+1)
# 
# plot(x[2,],type="l",col=2,lwd=2)
# lines(x[100,],col=3,lwd=2)

# Useful fun --------------------------------------------------------------

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Common domain  ----------------------------------------------------------

foo <- x$x.miss
0==sum(foo == 0, na.rm = T) # THIS HAVE TO BE ZERO
foo[is.na(foo)] <- 0
dim(foo)
jnk <- apply(foo, 2, prod)
length(jnk)
ext <- range(which(jnk != 0))
ext

#the common domain is obtained more simply by range(x$timepoints)

gr <- seq(0,1,length.out = 150) # this is the same as grid

# Plot --------------------------------------------------------------------
dat <- x$x.miss

# Actual plot
#pdf("simul-gen2.pdf", width = 15, height = 8)


X11()
# (rosa, blue) = (LN, U) = (otherwise, aneurysm at the terminal bifurcation)
colo   <- c("#ff9900", "#0868ac") 
colo.t <- c(rgb(215/255,181/255,216/255,.25), rgb(8/255,104/255,172/255,.25))


matplot(gr[6:145],t(xbind[,6:145]), type = "l", main = "Partially Observed Data with Error",
        lty = 2, col = rgb(0,0,0,.2), 
        axes = TRUE, xlab = "", ylab = "", bty="n")
matplot(gr[6:145],t(dat[,6:145]), type = "l", lty = 1, 
        col = colo[y+1], 
        axes = F, ylab = "", lwd = .85, add = TRUE)
matplot(gr[6:145],t(dat[,6:145]), type = "l",
        lty = 1, col = colo[y+1], 
        axes = F, ylab = "", lwd = .85, add = TRUE); 
rect(gr[ext[1]], -10, gr[ext[2]], +10, col = rgb(0,0,0,.05), border = NA)
legend("topleft", paste("Group", c("A", "B")), col = colo, bty = "n", lty = 1, lwd = 4, cex = .8)
for (i in 1:nrow(dat)){
  jnk <- dat[i,]
  foo <- range(which(!is.na(jnk)))
  xjnk <- gr[range(which(!is.na(jnk)))]
  if(xjnk[1]>gr[5] & xjnk[2]<gr[146]) {
  points(xjnk[1], jnk[foo[1]], pch = 21, col = rgb(0,0,0,.3), 
         bg = rgb(0,0,0,.1),
         cex = .85)
  points(xjnk[2], jnk[foo[2]], pch = 21, col = rgb(0,0,0,.3), 
         bg = rgb(0,0,0,.1),
         cex = .85)
  }
}
#dev.off()


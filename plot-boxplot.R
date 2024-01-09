# Prelim ------------------------------------------------------------------

library(ggplot2)
library(ggthemes)
library(ggjoy)
library(reshape2)

# Data --------------------------------------------------------------------

#load("simulation.RData")


# Prepare -----------------------------------------------------------------

tutti <- list()
tipi <- c("Without Missing", "James et al.", "Yao et al.", "Huang et al.")

for (j in 1:4){
  mydati <- t(sapply(ao[,j], unlist)) # THIS SELECT THE MIN CLASSIFICATION ERROR OF ALL THE SIMULATIONS FOR METHOD j
  nomi <- c()
  for(i in 1:(dim(mydati)[2])){
    nomi[i] <- 10^i
  }
  colnames(mydati) <- nomi
  dat <- stack(as.data.frame(mydati))  
  dat$Technique <- tipi[j]
  tutti[[j]] <- dat
}

# Put all together
dat <- rbind(tutti[[1]], tutti[[2]], tutti[[3]], tutti[[4]])
str(dat)
dat$Technique <- ordered(dat$Technique, 
                         levels = c("James et al.","Yao et al.", "Huang et al.",
                                     "Without Missing"))



### my changes

# bad=rbind(c(30,1),c(35,1),c(45,3),c(47,3),c(52,1),c(56,10),c(70,2),c(76,3),c(87,2),c(89,2),c(100,1))
# 
# for (i in 1:dim(bad)[1]) {
# 
#     err_SB[bad[i,1],bad[i,2]]=NA
#     
# }

colMeans(err_t)
colMeans(err_james)
colMeans(err_yao)
colMeans(err_buja)
colMeans(err_SB)

errors=as.data.frame( x = rbind(err_t, err_james, err_yao, err_buja, err_SB))

# errors=errors/S

errors$Technique =  as.factor(c(rep("Without Missing",S),rep("James et al.",S),rep("Yao et al.",S), rep("Huang et al.",S),
                                       rep("SB_VDFR",S)))

errors$Technique =  ordered(errors$Technique,levels = c("James et al.","Yao et al.", "Huang et al.","Without Missing", "SB_VDFR"))

dat=melt(errors, id= "Technique")
  
# Plot --------------------------------------------------------------------

bp <- ggplot(dat) + 
  geom_boxplot(aes(x = variable, y = value, fill = Technique), 
               outlier.shape = NA, alpha = .6 ) + 
  labs(x = "Domain Extention", y = "Leave-One-Out CV error (100 replications)",
       caption = " ") + ylim(0,50) + 
  theme_hc() + theme(axis.ticks.x = element_blank(), legend.position="top") +
  scale_x_discrete(labels = c("Common Domain",rep("",maxl-1), "Full Domain", "")) +
  scale_fill_manual(values = c("#4a1486", "#807dba", "#bcbddc", "#efedf5", "#bcbbcbdc")) 
# +
#   stat_summary(fun=mean, geom="point", size=2, color="black")



### LDA PLOTS

errors_lda=as.data.frame( x = rbind(err_t_lda, err_james_lda, err_yao_lda, err_buja_lda, err_SB))

errors_lda$Technique =  as.factor(c(rep("Without Missing",S),rep("James et al.",S),rep("Yao et al.",S), rep("Huang et al.",S),
                                rep("SB_VDFR",S)))

errors_lda$Technique =  ordered(errors_lda$Technique,levels = c("James et al.","Yao et al.", "Huang et al.","Without Missing", "SB_VDFR"))

dat_lda=melt(errors_lda, id= "Technique")

# Plot --------------------------------------------------------------------

bp_lda <- ggplot(dat_lda) + 
  geom_boxplot(aes(x = variable, y = value, fill = Technique), 
               outlier.shape = NA, alpha = .6 ) + 
  labs(x = "Domain Extention", y = "Leave-One-Out CV error (100 replications)",
       caption = " ") + ylim(0,50) + 
  theme_hc() + theme(axis.ticks.x = element_blank(), legend.position="top") +
  scale_x_discrete(labels = c("Common Domain",rep("",maxl-1), "Full Domain", "")) +
  scale_fill_manual(values = c("#4a1486", "#807dba", "#bcbddc", "#efedf5", "#bcbbcbdc")) 
# +
#   stat_summary(fun=mean, geom="point", size=2, color="black")


bp <- ggplot(dat) + 
  geom_boxplot(aes(x = ind, y = values, fill = Technique), 
               outlier.shape = NA, alpha = .6 ) + 
  labs(x = "Domain Extention", y = "Leave-One-Out CV error (20 replications)",
       caption = " ") + ylim(5,40) + 
  theme_hc() + theme(axis.ticks.x = element_blank(), legend.position="top") +
  scale_x_discrete(labels = c("Common Domain",rep("",maxl-1), "Full Domain", "")) +
  scale_fill_manual(values = c("#4a1486", "#807dba", "#bcbddc", "#efedf5"))


# Show
#pdf("allbox.pdf", height = 8, width = 11)
X11()
bp
#dev.off()

# Prelim ------------------------------------------------------------------

library(ggplot2)
library(ggthemes)
library(ggjoy)
library(reshape2)

# Data --------------------------------------------------------------------

aux_2=load("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/2do paper/codigo/Paper 2 Classification/Classification-via-SB-VDFR-model/Final simulations/50_iter_13_07_2023_0,6 and 0,1_Kneip - copia.RData")
# aux_2=load("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/2do paper/codigo/Paper 2 Classification/Classification-via-SB-VDFR-model/Final simulations/50_iter_13_07_2023_0,6 and 0,4_Kneip - copia.RData")

err_kneip=matrix(nrow=50,ncol=10)

err_kneip[,10]=get(aux_2[which(aux_2=="miss_classification_error")])

load("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/2do paper/codigo/Paper 2 Classification/Classification-via-SB-VDFR-model/Final simulations/50_iter_06_06_2023_good_cutting_point_unif.RData")
# load("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/2do paper/codigo/Paper 2 Classification/Classification-via-SB-VDFR-model/Final simulations/50_iter_26_06_2023_good_cutting_point 0,6 and 0,4.RData")

errors=as.data.frame( x = rbind(err_t, err_james, err_yao, err_buja, err_SB))

aux=load("C:/Users/user/Desktop/Trabajo/Escuela/Doctorado/Pavel/Tesis/2do paper/codigo/Paper 2 Classification/Classification-via-SB-VDFR-model/50_iter_08_05_2023.RData")

  err_t=get(aux[which(aux=="err_t")])
  err_james=get(aux[which(aux=="err_james")])
  err_yao=get(aux[which(aux=="err_yao")])
  err_buja=get(aux[which(aux=="err_buja")])

errors[1:(4*S),]=rbind(err_t, err_james, err_yao, err_buja)

errors=rbind(errors, err_kneip)

errors$Technique =  as.factor(c(rep("Without Missing",S),rep("James et al.",S),rep("Yao et al.",S), rep("Huang et al.",S),
                                rep("POFRM",S),rep("Kneip",S)))

errors$Technique =  ordered(errors$Technique,levels = c("James et al.","Yao et al.", "Huang et al.","Without Missing", "Kneip", "POFRM"))

dat=melt(errors, id= "Technique")

dat$value=log(dat$value+1)

# Plot --------------------------------------------------------------------

bp <- ggplot(dat) + 
  geom_violin(aes(x = variable, y = value, fill = Technique), alpha = .6, scale="width") + 
  labs(x = "Domain Extention", y = "Leave-One-Out CV error (50 replications)",
       caption = " ") + ylim(0,4) + 
  theme_hc() + theme(axis.ticks.x = element_blank(), legend.position="top") +
  scale_x_discrete(labels = c("Common Domain",rep("",maxl-1), "Full Domain", "")) 

# +
#   scale_fill_manual(values = c("#bcbbcbdc", "#bcbddc", "#807dba", "#efedf5", "lightgreen"))
# +
#   scale_fill_manual(values = c("#4a1486", "#807dba", "#bcbddc", "#efedf5", "#bcbbcbdc"))

# bp <- ggplot(dat) + 
#   geom_boxplot(aes(x = variable, y = value, fill = Technique), 
#                outlier.shape = NA, alpha = .6, outlier.color="red", outlier.size=2 ) + 
#   labs(x = "Domain Extention", y = "Leave-One-Out CV error (100 replications)",
#        caption = " ") + ylim(0,40) + 
#   theme_hc() + theme(axis.ticks.x = element_blank(), legend.position="top") +
#   scale_x_discrete(labels = c("Common Domain",rep("",maxl-1), "Full Domain", "")) +
#   scale_fill_manual(values = c("#4a1486", "#807dba", "#bcbddc", "#efedf5", "#bcbbcbdc")) 
# +
#   stat_summary(fun=mean, geom="point", size=2, color="black")

err_SB=errors[((4*S)+1):(5*S),1:10]

colMeans(err_t)
colMeans(err_james)
colMeans(err_yao)
colMeans(err_buja)
colMeans(err_SB)
colMeans(err_kneip)

apply(err_t, 2, var)
apply(err_james, 2, var)
apply(err_yao, 2, var)
apply(err_buja, 2, var)
apply(err_SB, 2, var)
apply(err_kneip, 2, var)

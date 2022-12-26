#install.packages("psych", dependencies=TRUE)
library(psych)
library(FactoMineR)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

select <- dplyr::select
rename <- dplyr::rename
library(car)
library(GPArotation)
library(nFactors)




###################################
### FUNCTION TO REMOVE OUTLIERS ###
###################################

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}




#############################################
### DIRECTORIES TO SAVE IMAGES AND TABLES ###
#############################################

## We set the base DIR
base_dir =dirname(rstudioapi::getSourceEditorContext()$path)
setwd(base_dir)
img_dir = paste(base_dir,  "/R_output", sep="")


#################
### LOAD FILE ###
#################

my_data <-read.csv("data/2.1_output_frequencies_zscore.csv", header=T)


######################
###  TRIM THE DATA ###
######################

#Data is z-scored and outliers -3>Z>3 are already removed

dat1 <- my_data %>%
  mutate(
 
    hope = as.numeric(hope_out),
    rage = as.numeric(rage_out),
    ) %>%
  dplyr::select(hope,rage)


psych::describe(dat1)
boxplot(dat1)

dat2= dat1

#build a correlation matrix
dat3 = cor(dat2,use="pairwise.complete.obs")
dat3


# Determine Number of Factors to Extract
ev <- eigen(dat3) # get eigenvalues
ap <- parallel(subject=nrow(dat3),var=ncol(dat3),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)


# Factor Analysis
fit = fa(dat3,nfactors=2, rotate="oblimin", scores=TRUE, residuals=FALSE, SMC=TRUE, covar=FALSE,missing=TRUE,impute="median")
fit

load <- fit$loadings[,1:2] 

df_out <- as.data.frame(load)
theme<-theme(plot.title = element_text(hjust = 0.5,size=18),legend.position = "none",axis.title=element_text(size=18), panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black", size=18),axis.text.y=element_text(colour="black",size=18),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

pdf(file = paste(img_dir, "/Factor Analysis.pdf", sep=""))

ggplot(df_out,aes(x=MR1,y=MR2,color=row.names(df_out),label=row.names(df_out)))+
  geom_point()+theme + geom_text(size=5, nudge_y =0.05, nudge_x =0.1) + xlim(-0.3,0.8)+
  geom_segment(df_out,mapping=aes(x=0,y=0,xend = MR1, yend =MR2), size  =1.2,arrow = arrow())+
  geom_hline(yintercept=0, linetype="dashed", color = "blue")+
  ggtitle('Exploratory Factor Analysis France')

dev.off()

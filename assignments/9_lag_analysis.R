#################
### LIBRARIES ###
#################

library(tseries)
library(robustbase)
library(lmtest)
library(astsa)
library(robust)
library(rr2)
library(sjPlot)
library(chron)
library(smooth)
library(timeSeries)
library(Mcomp)
library(lme4)
library(glmulti)
library(LMERConvenienceFunctions)
library(MASS)
library(psych)
library(Hmisc)
library(ggplot2)
library(rcompanion)
require(cowplot)
library(blme)
library(emmeans)
library(lsmeans)

library(dplyr)
select <- dplyr::select
rename <- dplyr::rename
library(stargazer)
library(car)
library(reshape2)
library(nlme)


#############################################
### DIRECTORIES TO SAVE IMAGES AND TABLES ###
#############################################

## We set the base DIR
base_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(base_dir)
img_dir = paste(base_dir,  "/R_output", sep = "")
tab_dir = paste(base_dir,  "/R_output", sep = "")

my_data <-
  read.csv("data/3.3_output_complete_data_lags.csv", header = T)

################
###  1. GDP  ###
################


## Initialize  all variables and transform them

dat1 <- my_data %>%
  mutate(
    year2 = year,
    year = as.numeric(scale(year)),
    year = as.numeric(scale(year)),
    GDP = as.numeric(scale(GDPpc)),
    GDPm1 = as.numeric(scale(GDPpc.1)),
    GDPm2 = as.numeric(scale(GDPpc.2)),
    GDPm3 = as.numeric(scale(GDPpc.3)),
    GDPm4 = as.numeric(scale(GDPpc.4)),
    GDPm5 = as.numeric(scale(GDPpc.5)),
    GDPm6 = as.numeric(scale(GDPpc.6)),
    GDPm7 = as.numeric(scale(GDPpc.7)),
    GDPm8 = as.numeric(scale(GDPpc.8)),
    GDPm9 = as.numeric(scale(GDPpc.9)),
    GDPm10 = as.numeric(scale(GDPpc.10)),
    GDPm11 = as.numeric(scale(GDPpc.11)),
    GDPm12 = as.numeric(scale(GDPpc.12)),
    GDPm13 = as.numeric(scale(GDPpc.13)),
    GDPm14 = as.numeric(scale(GDPpc.14)),
    GDPm15 = as.numeric(scale(GDPpc.15)),
    GDPm16 = as.numeric(scale(GDPpc.16)),
    GDPm17 = as.numeric(scale(GDPpc.17)),
    GDPm18 = as.numeric(scale(GDPpc.18)),
    GDPm19 = as.numeric(scale(GDPpc.19)),
    GDPm20 = as.numeric(scale(GDPpc.20)),
    
    GDPp1 = as.numeric(scale(GDPpc1)),
    GDPp2 = as.numeric(scale(GDPpc2)),
    GDPp3 = as.numeric(scale(GDPpc3)),
    GDPp4 = as.numeric(scale(GDPpc4)),
    GDPp5 = as.numeric(scale(GDPpc5)),
    GDPp6 = as.numeric(scale(GDPpc6)),
    GDPp7 = as.numeric(scale(GDPpc7)),
    GDPp8 = as.numeric(scale(GDPpc8)),
    GDPp9 = as.numeric(scale(GDPpc9)),
    GDPp10 = as.numeric(scale(GDPpc10)),
    GDPp11 = as.numeric(scale(GDPpc11)),
    GDPp12 = as.numeric(scale(GDPpc12)),
    GDPp13 = as.numeric(scale(GDPpc13)),
    GDPp14 = as.numeric(scale(GDPpc14)),
    GDPp15 = as.numeric(scale(GDPpc15)),
    GDPp16 = as.numeric(scale(GDPpc16)),
    GDPp17 = as.numeric(scale(GDPpc17)),
    GDPp18 = as.numeric(scale(GDPpc18)),
    GDPp19 = as.numeric(scale(GDPpc19)),
    GDPp20 = as.numeric(scale(GDPpc20)),
    Love = as.numeric(ratio_feelings2_Out)
  ) %>%
  
  dplyr::select(
    year,
    year2,
    GDP ,
    Love,
    text_per_year,
    GDPm1,
    GDPp1,
    GDPm2,
    GDPp2,
    GDPm3,
    GDPp3,
    GDPm4,
    GDPp4,
    GDPm5,
    GDPp5,
    GDPm6,
    GDPp6,
    GDPm7,
    GDPp7,
    GDPm8,
    GDPp8,
    GDPm9,
    GDPp9,
    GDPm10,
    GDPp10,
    GDPm11,
    GDPp11,
    GDPm12,
    GDPp12,
    GDPm13,
    GDPp13,
    GDPm14,
    GDPp14,
    GDPm15,
    GDPp15,
    GDPm16,
    GDPp16,
    GDPm17,
    GDPp17,
    GDPm18,
    GDPp18,
    GDPm19,
    GDPp19,
    GDPm20,
    GDPp20
  )




## Eliminate NAs
dat3 = na.omit(dat1)



############################################
### 1.2. MODELS WITH THE LAGS T-20, T+20 ###
############################################


###################################
### transform data using BoxCox ###
###################################

dat3$Love2 = dat3$Love + 4

Box = boxcox(
  Love2 ~ year + GDP + GDPm1 + GDPp1 + GDPm2 + GDPp2 + GDPm3 + GDPp3 + GDPm4 + GDPp4 + GDPm5 + GDPp5 + GDPm6 + GDPp6 + GDPm7 + GDPp7 + GDPm8 + GDPp8 + GDPm9 + GDPp9 + GDPm10 + GDPp10,
  #+ GDPm11+GDPp11+GDPm12+GDPp12+GDPm13+GDPp13+GDPm14+GDPp14+GDPm15+GDPp15+GDPm16+GDPp16+GDPm17+GDPp17+GDPm18+GDPp18+GDPm19+GDPp19+GDPm20+GDPp20,
  data = dat3,
  lambda = seq(-6, 6, 0.1)
)

Cox = data.frame(Box$x, Box$y)

Cox2 = Cox[with(Cox, order(-Cox$Box.y)), ]

Cox2[1, ]

lambda = Cox2[1, "Box.x"]
lambda ## optimal lambda 1

dat3$Love3 = scale(((dat3$Love2 ** lambda) - 1) / lambda)


####################################################################
### MODEL controlling for autocorrelations in the time dimension ###
####################################################################


## First we calculate 5 models with lags separated by 5 years intervals to avoid co-linearities


full_model_love_GDP = gls(
  Love3 ~ year + GDP + GDPm5 + GDPp5 +  GDPm10 + GDPp10 +
    GDPm15 + GDPp15 + GDPm20 + GDPp20 + GDPm1 + GDPp1 + GDPm6 + GDPp6 +
    GDPm11 + GDPp11 + GDPm16 + GDPp16 + GDPm2 + GDPp2 +  GDPm7 + GDPp7 +
    GDPm12 + GDPp12 + GDPm17 + GDPp17 +  GDPm3 + GDPp3 +  GDPm8 + GDPp8 +
    GDPm13 + GDPp13 + GDPm18 + GDPp18 + GDPm4 + GDPp4 + GDPm9 + GDPp9 +
    GDPm14 + GDPp14 + GDPm19 + GDPp19
  ,
  data = dat3,
  corCAR1(form =  ~ year),
  method = 'ML'
)


## Then we do model selection with stepAIC to find the best model of each above

best_model_love_GDP <-
  stepAIC(
    full_model_love_GDP,
    scope = list(lower = ~ year),
    direction = "both",
    k = log(nrow(dat3)),
    trace = TRUE
  )



Anova(best_model_love_GDP, type = 'II')
summary(best_model_love_GDP, type = 'II')
confint(best_model_love_GDP)


BIC(best_model_love_GDP)

res_best_model_love_GDP = residuals(best_model_love_GDP)
shapiro.test(res_best_model_love_GDP)

## Export models##
setwd(tab_dir)

stargazer(
  best_model_love_GDP,
  title = "Lags",
  align = TRUE,
  type = 'html',
  out = paste(img_dir, "/lag_analysis.htm", sep = ""),
  flip = TRUE
)



## Export diagnostics ##

setwd(img_dir)
png(
  file = paste(img_dir, "/lags_diagnostics.png", sep = ""),
  width = 12,
  pointsize = 16,
  height = 7,
  units = "in",
  res = 300
)
par(mfrow = c(1, 2))
acf(res_best_model_love_GDP, main = "GDP-Love residual autocorrelation")
plotNormalHistogram(res_best_model_love_GDP, main = "GDP-Love Residuals")

dev.off()

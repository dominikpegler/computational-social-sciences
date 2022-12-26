

############################
### Import all libraries ###
############################

library (lme4)
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
library(ggthemes)
library(lmerTest)
library(robustlmm)
library(ggplot2)
library(tseries)
library(nlme)
library(gridExtra)
library(robustlmm)
library(dplyr)
select <- dplyr::select
rename <- dplyr::rename

library(car)
library(reshape2)
library(glmulti)
library(sjPlot)
library(stargazer)


###################################
### FUNCTION TO REMOVE OUTLIERS ###
###################################

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
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
base_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(base_dir)
img_dir = paste(base_dir,  "/R_output", sep = "")
tab_dir = paste(base_dir,  "/R_output", sep = "")




#################
### LOAD FILE ###
#################

my_data <- read.csv("data/3.1_output_complete_dataset.csv", header = T)

######################
###  TRIM THE DATA ###
######################

## Note: The ratios ranged from -1 to 1, but for statistics we made greater or equal to 0,
## so that we could apply power transformations as suggested by box cox


dat1 <- my_data %>%
  mutate(
    movies = as.factor(movies),
    year2 = as.numeric(year),
    year = as.numeric(scale(year)),
    
    
    hope = as.numeric(hope_out),
    rage = as.numeric(rage_out),
    hope.to.rage.ratio = as.numeric(ratio_hope_rage2_out),

    GDPpc = as.numeric(scale(GDPpc))
  ) %>%
  select(
    movies,
    year,
    hope,
    rage,
    hope.to.rage.ratio,
    GDPpc
  )


##############################################################################
#################              HOPE  SECTION            ######################
##############################################################################

#############################
### GET DESCRIPTIVE STATS ###
#############################

psych::describe(dat1)

plot(dat1$hope.to.rage.ratio ~ dat1$year)
plotNormalHistogram(dat1$hope.to.rage.ratio)


##############################################################################
#################              HOPE  SECTION            ######################
##############################################################################


cor.test(dat1$hope.to.rage.ratio, dat1$GDPpc)

#############################
### SIMPLE TIME+GDP MODEL ###
#############################

###################################
### transform data using BoxCox ###
###################################


Box = boxcox(hope.to.rage.ratio ~ year + GDPpc,
             data = dat1,
             lambda = seq(-6, 6, 0.1))

Cox = data.frame(Box$x, Box$y)

Cox2 = Cox[with(Cox, order(-Cox$Box.y)), ]

Cox2[1, ]

lambda = Cox2[1, "Box.x"]
lambda ## optimal lambda 0.7

dat1$hope = scale(((dat1$hope.to.rage.ratio ** lambda) - 1) / lambda)


#######################
### ECONOMIC MODELS ###
#######################


### GDP model without time ###

GDP_model = lmer(hope ~ GDPpc + (1 | movies), data = dat1)
Anova(GDP_model, type = '2')
summary(GDP_model)
confint(GDP_model)


res_GDP_model = residuals(GDP_model)
shapiro.test(res_GDP_model)

acf(res_GDP_model)
BIC(GDP_model)



### GDP model with time ###

GDP_time_model = lmer(hope ~ year + GDPpc + (1 |
                                               movies), data = dat1)
Anova(GDP_time_model, type = '2')
summary(GDP_time_model)
confint(GDP_time_model)


res_GDP_time_model = residuals(GDP_time_model)
shapiro.test(res_GDP_time_model)


vif(GDP_time_model)
acf(res_GDP_time_model)
BIC(GDP_time_model)

class(GDP_model) <- "lmerMod"
class(GDP_time_model) <- "lmerMod"

stargazer(
  GDP_model,
  GDP_time_model,
  title = "Lags",
  align = TRUE,
  type = 'html',
  out = paste(img_dir, "/lmm.htm", sep = ""),
  flip = TRUE
)

## Plot model diagnostics

setwd(img_dir)
png(
  file = paste(img_dir, "/lmm_diagnostics.png", sep = ""),
  width = 9,
  height = 3.5,
  units = "in",
  res = 300
)
par(mfrow = c(1, 3))
acf(res_GDP_time_model, main = "GDP-Hope residual autocorrelation")
pacf(res_GDP_time_model, main = "GDP-Hope residual Partial autocorrelation")
plotNormalHistogram(res_GDP_time_model, main = "GDP-Hope Residuals")
dev.off()

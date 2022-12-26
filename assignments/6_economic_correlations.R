

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

library(dplyr)
select <- dplyr::select
rename <- dplyr::rename

library(car)
library(reshape2)
library(glmulti)

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
    hope = as.numeric(hope_out),
    rage = as.numeric(rage_out),
    hope.to.rage.ratio = as.numeric(ratio_hope_rage2_out),
    GDPpc = as.numeric(scale(GDPpc))
  ) %>%
  dplyr::select(hope.to.rage.ratio, year, GDPpc)

############################
### GET RAW CORRELATIONS ###
############################

mydata = dat1


cormat <-
  round(cor(mydata, method = c("pearson"), use = "pairwise.complete.obs"), 2)
melted_cormat <- melt(cormat)

mycor <- rcorr(as.matrix(mydata), type = "pearson")

###########################
### GET  R and P-VALUES ###
###########################

round(mycor$r, 3)
round(mycor$P, 4)
round(mycor$n, 0)

library(corrplot)

res1 <- cor.mtest(mydata, conf.level = .95)

theme <-
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.title = element_text(size = 20),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "line")
  )


pdf(file = paste(img_dir, "/SocioDemographic_Correlations.pdf", sep = ""))


col <-
  colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(
  cormat,
  method = "color",
  col = col(200),
  type = "upper",
  order = "original",
  number.cex = 1,
  addCoef.col = "black",
  # Add coefficient of correlation
  tl.col = "black",
  tl.srt = 45,
  tl.cex = 1.2,
  # Text label color and rotation
  # Combine with significance
  p.mat = mycor$P,
  sig.level = 0.05,
  insig = "blank",
  # hide correlation coefficient on the principal diagonal
  diag = FALSE
)

dev.off()

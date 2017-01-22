#'Calculate Durbin's h-statistic
#'
#'Calculates Durbin's h-statistic for autoregressive models.
#'@param model The model being assessed
#'@return Numeric Durbin's h statistic 
#'@details
#'Using the Durbin-Watson (DW) test for autoregressive models (like ECM) is inappropriate because the 
#'DW test itself tests for first order autocorrelation. Since an ECM model inherently uses the lag term
#'of the target variable as one of the predictors, Durbin's h-statistic should be used to test for 
#'autocorrelation. If Durbin's h-statistic is greater than 1.96, it is likely that autocorrelation exists.
#'
#'@seealso \code{lm}
#'@examples
#'#Use ecm to predict Wilshire 5000 index based on corporate profits, 
#'#Federal Reserve funds rate, and unemployment rate
#'data(Wilshire)
#'
#'#Use 2014-12-01 and earlier data to build models
#'trn <- Wilshire[Wilshire$date<='2014-12-01',]
#'
#'#Assume all predictors are needed in the equilibrium and transient terms of ecm
#'xeq <- xtr <- trn[c('CorpProfits', 'FedFundsRate', 'UnempRate')]
#'model1 <- ecm(trn$Wilshire5000, xeq, xtr)
#'
#'#Check Durbin's h-statistic on model1
#'durbinH(model1)
#'#The h-statistic is 4.55, which means there is likely autocorrelation in the data.
#'
#'@export
#'@importFrom car durbinWatsonTest
durbinH <- function(model){
  d <- car::durbinWatsonTest(model)
  n <- length(model$fitted.values) + 1
  v <- summary(model)$coef[nrow(summary(model)$coef),2]^2
  durbinH <- (1 - 0.5 * d$dw) * sqrt(n / (1 - n*v))
  return(durbinH)
}
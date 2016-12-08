#'Backwards selection to build an error correction model
#'
#'Much like the ecm function, this builds an error correction model.
#'However, it uses backwards selection to select the optimal predictors based on lowest AIC or BIC, rather than using all predictors.
#'ecmback has the same parameters and output as ecm.
#'@param y The target variable
#'@param xeq The variables to be used in the equilibrium term of the error correction model
#'@param xtr The variables to be used in the transient term of the error correction model
#'@param criterion Whether AIC (default) or BIC should be used to select variables 
#'@return an lm object representing an error correction model using backwards selection
#'@seealso \code{lm}
#'@examples
#'#Use ecm to predict Wilshire 5000 index based on corporate profits, 
#'#Federal Reserve funds rate, and unemployment rate
#'data(Wilshire)
#'
#'#Use 2014-12-01 and earlier data to build models
#'trn <- Wilshire[Wilshire$date<='2014-12-01',]
#'
#'#Use backwards selection to choose which predictors are needed 
#'xeq <- xtr <- trn[c('CorpProfits', 'FedFundsRate', 'UnempRate')]
#'modelback <- ecmback(trn$Wilshire5000, xeq, xtr)
#'print(modelback)
#'#Backwards selection chose CorpProfits in the equilibrium term, 
#'#CorpProfits and UnempRate in the transient term.
#'
#'@export
#'@importFrom stats lm complete.cases step
ecmback <- function (y, xeq, xtr, criterion="AIC") 
{
  if (missing(xeq)) {
    xtrnames <- names(xtr)
  }
  else {
    xeqnames <- names(xeq)
    if (class(xeq) != "data.frame") {
      xeqnames <- deparse(substitute(xeq))
      xeqnames <- substr(xeqnames, regexpr("\\$", xeqnames) + 
                           1, nchar(xeqnames))
    }
    if (missing(xtr)) {
      xtrnames <- xeqnames
      xtr <- xeq
    }
    else {
      xtrnames <- names(xtr)
    }
    xeqnames <- paste0(xeqnames, "Lag1")
    xeq <- as.data.frame(xeq)
    ifelse(ncol(xeq) > 1, xeq <- rbind(rep(NA, ncol(xeq)), xeq[1:(nrow(xeq) - 1), ]), xeq <- data.frame(c(NA, xeq[1:(nrow(xeq) - 1), ])))
  }
  
  xtrnames <- paste0("delta", xtrnames)
  xtr <- as.data.frame(xtr)
  xtr <- data.frame(apply(xtr, 2, diff, 1))
  dy <- diff(y, 1)
  if (!missing(xeq)) {
    yLag1 <- y[1:(length(y) - 1)]
    x <- cbind(xtr, xeq[complete.cases(xeq), ])
    x <- cbind(x, yLag1)
    names(x) <- c(xtrnames, xeqnames, "yLag1")
  }
  else {
    x <- xtr
    names(x) <- xtrnames
  }
  
  full <- lm(dy ~ ., data = x)
  null <- lm(dy ~ yLag1, data = x)
  if(criterion=='AIC'){
    k=2
  } else if(criterion=='BIC'){
    k=log(nrow(x))
  }
  
  ecm <- step(full, data = x, scope = list(upper = full, lower = null), direction = "backward", k=k, trace=0)
  return(ecm)
}

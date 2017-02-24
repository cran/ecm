#'Predict using an ecm object
#'
#'Takes an ecm object and uses it to predict based on new data.
#'@param ecm ecm object used to make predictions
#'@param newdata Data frame to on which to predict
#'@param init Initial value for prediction
#'@return Numeric predictions on new data based ecm object
#'@details 
#'Since error correction models only model the change in the target variable, an initial value must be specified.
#'@examples
#'data(Wilshire)
#'
#'#Rebuilding model1 from ecm example
#'trn <- Wilshire[Wilshire$date<='2014-12-01',]
#'xeq <- xtr <- trn[c('CorpProfits', 'FedFundsRate', 'UnempRate')]
#'model1 <- ecm(trn$Wilshire5000, xeq, xtr)
#'
#'#Use 2014-12-01 and onwards data as test data to predict
#'tst <- Wilshire[Wilshire$date>='2014-12-01',]
#'
#'#predict on tst using model1 and initial FedFundsRate
#'tst$model1Pred <- ecmpredict(model1, tst, tst$Wilshire5000[1])
#'
#'@export
#'@importFrom stats predict
ecmpredict <- function (ecm, newdata, init) {
  if (sum(grepl('^delta', names(ecm$coefficients))) >= 1) {
    form <- names(ecm$coefficients)
    xtrnames <- form[grep("^delta", form)]
    xtrnames <- substr(xtrnames, 6, max(nchar(xtrnames)))
    xtr <- newdata[which(names(newdata) %in% xtrnames)]
    xtr <- data.frame(apply(xtr, 2, diff, 1))
    names(xtr) <- paste0("delta", names(xtr))
    xtrnames <- names(xtr)
  }
  
  if (sum(grepl('Lag1$', names(ecm$coefficients))) > 1) {
    form <- names(ecm$coefficients)
    xeqnames <- form[grep("^(?!delta).*", form, perl = T)]
    xeqnames <- xeqnames[-c(1, length(xeqnames))]
    xeqnames <- substr(xeqnames, 1, unlist(lapply(gregexpr("Lag", xeqnames), function(x) x[length(x)])) - 1)
    xeq <- newdata[which(names(newdata) %in% xeqnames)]
    names(xeq) <- paste0(names(xeq), "Lag1")
    xeqnames <- names(xeq)
  }
  
  if (exists('xeq')) {
    if (ncol(xeq) > 1) {
      xeq <- rbind(rep(NA, ncol(xeq)), xeq[1:(nrow(xeq) - 1), ])
    } else if (ncol(xeq) == 1) {
      xeq <- data.frame(c(NA, xeq[1:(nrow(xeq) - 1), ]))
    }
  }   
  
  if (exists('xeq') & exists('xtr')) {
    if ((is.na(xeq))/nrow(xeq) == 1) {
      x <- xtr
    } else {
      x <- cbind(xtr, xeq[complete.cases(xeq), ])
    }
    x$yLag1 <- init
    names(x) <- c(xtrnames, xeqnames, "yLag1")
  } else if (!exists('xeq') & exists('xtr')) {
    x <- xtr 
    x$yLag1 <- init
    names(x) <- c(xtrnames, "yLag1")
  } else if (exists('xeq') & !exists('xtr')) {
    x <- xeq[complete.cases(xeq),]
    x$yLag1 <- init
    names(x) <- c(xeqnames, "yLag1")
  }
  
  ecmpred <- predict(ecm, x[1, ])
  for (i in 2:nrow(x)) {
    x$yLag1[i] <- x$yLag1[i - 1] + ecmpred
    ecmpred <- predict(ecm, x[i, ])
  }
  ecmpred <- predict(ecm, x)
  ecmpred <- cumsum(c(init, ecmpred))
  return(ecmpred)
}




cv <- function(x){
  sd(x)/abs(mean(x))
}

synchrony <- function (x) {
  var(rowSums(x)) / (sum(apply(x, 2, sd)) ^ 2)
}

synchrony_detrended <- function (x) {
  var(residuals(lm(rowSums(x)~c(1:nrow(x))))) / (sum(apply(apply(x, 2, function(y) residuals(lm(y~c(1:length(y))))), 2, sd)) ^ 2)
}


pe_mv <- function(x, type = c("linear", "linear_robust", "quadratic",
                              "linear_quad_avg",  "linear_detrended", "loess_detrended"), ci =
                    FALSE, na.rm = FALSE) {
  
  type <- type[1]
  
  if(!type %in% c("linear", "linear_robust", "quadratic",
                  "linear_quad_avg", "linear_detrended", "loess_detrended")){
    stop("not a valid type")
  }
  
  if(!type %in% c("linear", "linear_robust", "linear_detrended",
                  "loess_detrended")){
    if(ci == TRUE){
      warning("Confidence intervals aren't supported for this type of
              mean-variance model. Setting ci = FALSE.")
    }
    ci <- FALSE
    }
  
  if(na.rm) x <- na.omit(x)
  
  total_nas <- sum(is.na(x))
  return_na <- ifelse(!na.rm & total_nas > 0, TRUE, FALSE)
  
  ## first get the means:
  m <- apply(x, 2, mean)
  single_asset_mean <- mean(rowSums(x))
  
  cv_portfolio <- cv(rowSums(x))
  
  ## now detrend if desired:
  if(type == "linear_detrended") {
    ## first get cv of detrended portfolio abundance:
    sd_portfolio <- sd(residuals(lm(rowSums(x)~c(1:nrow(x)))))
    mean_portfolio <- mean(rowSums(x))
    cv_portfolio <- sd_portfolio / mean_portfolio
    ## now detrend:
    x <- apply(x, 2, function(y) residuals(lm(y~c(1:length(y)))))
  }
  if(type == "loess_detrended") {
    ## first get CV of detrended portfolio abundance:
    sd_portfolio <- sd(residuals(loess(rowSums(x)~c(1:nrow(x)))))
    mean_portfolio <- mean(rowSums(x))
    cv_portfolio <- sd_portfolio / mean_portfolio
    ## now detrend:
    x <- apply(x, 2, function(y) residuals(loess(y~c(1:length(y)))))
  }
  
  ## and get the variances for the assets:
  v <- apply(x, 2, var)
  
  log.m <- log(m)
  log.m[which(log.m=="-Inf")]<-NA
  log.v <- log(v)
  log.v[which(log.v=="-Inf")]<-NA

  d <- data.frame(log.m = log.m, log.v = log.v, m = m, v = v)
  
  taylor_fit <- switch(type[1], 
                       linear = {
                         lm(log.v ~ log.m, data = d)
                       },
                       linear_robust = {
                         robustbase::lmrob(log.v ~ log.m, data = d)
                       },
                       quadratic = {
                         nls(log.v ~ B0 + B1 * log.m + B2 * I(log.m ^ 2),
                             data = d, start = list(B0 = 0, B1 = 2, B2 = 0), lower =
                               list(B0 = -1e9, B1 = 0, B2 = 0), algorithm = "port")
                       },
                       linear_detrended = {
                         lm(log.v ~ log.m, data = d)
                       }, 
                       loess_detrended = {
                         lm(log.v ~ log.m, data = d)
                       },
                       linear_quad_avg = {
                         linear <- nls(log.v ~ B0 + B1 * log.m, data = d, start = list(B0 =
                                                                                         0, B1 = 2), lower = list(B0 = -1e9, B1 = 0), algorithm =
                                         "port")
                         quadratic <- nls(log.v ~ B0 + B1 * log.m + B2 * I(log.m ^ 2), data
                                          = d, start = list(B0 = 0, B1 = 2, B2 = 0), lower = list(B0 =
                                                                                                    -1e9, B1 = 0, B2 = 0), algorithm = "port")
                         MuMIn::model.avg(list(linear=linear, quad=quadratic), rank = MuMIn::AICc)
                       }
  )
  
  if(ci) {
    single_asset_variance_predict <- predict(taylor_fit, newdata =
                                               data.frame(log.m = log(single_asset_mean)), se = TRUE)
    single_asset_variance <- exp(single_asset_variance_predict$fit)
  } else {
    single_asset_variance_predict <- predict(taylor_fit, newdata =
                                               data.frame(log.m = log(single_asset_mean)), se = FALSE)
    single_asset_variance <- exp(single_asset_variance_predict)
  }
  
  cv_single_asset <- sqrt(single_asset_variance) / single_asset_mean
  pe <- as.numeric(cv_single_asset / cv_portfolio)
  
  if(ci == TRUE) {
    single_asset_variance_ci <- exp(single_asset_variance_predict$fit
                                    + c(-1.96, 1.96) * single_asset_variance_predict$se.fit)
    cv_single_asset_ci <- sqrt(single_asset_variance_ci) / single_asset_mean
    pe_ci <- as.numeric(cv_single_asset_ci / cv_portfolio)
    pe_ci <- pe_ci[order(pe_ci)] # make sure the lower value is first
    out <- list(pe = pe, ci = pe_ci, cv_single_asset=cv_single_asset, cv_portfolio=cv_portfolio)
  } else {
    out <- list(pe = pe, cv_single_asset=cv_single_asset, cv_portfolio=cv_portfolio)
  }
  
  if(return_na) out <- NA
  
  out
}

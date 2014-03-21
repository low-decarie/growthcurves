#' logistic_growth_nls
#'
#' Least-squares estimates of the parameters of a logistic growth curve
#'
#' @param time Numeric vector of time since start of measurement (currently does not use date format)
#' @param abundance Abundance (eg. Optical density, counts etc).
#' @param do.interval Logical: should a 95\% confidence interval be produced
#' @keywords growth
#' @examples
#' data(navicula_growth)
#' fitted.readings <- with(navicula_growth,
#'                        logistic_growth_nls(time=Time,
#'                                            abundance=ABS))
#' head(fitted.readings)
#' @export


logistic_growth_nls<-function(time, abundance, do.interval=T){

  start.values=c("K"=max(abundance, na.rm=T), "r"=1, "N0"=min(abundance, na.rm=T))
  
  
  culture.model <- try(nls(formula=abundance ~(K*N0*exp(r*time) ) / (K + N0 * (exp(r*time)-1)),
                           start=start.values,
                           na.action=na.exclude,
                           algorithm="port",
                           lower=c(0,0,0),
                           upper=c(2*max(abundance, na.rm=T),Inf, Inf)))
  
  fitted.readings <- data.frame(time,abundance)
  
  if(class(culture.model)=="try-error"){
    
    fitted.readings$logistic.nls.N0<-NA
    fitted.readings$logistic.nls.K<-NA
    fitted.readings$logistic.nls.K.lower<-NA
    fitted.readings$logistic.nls.K.upper<-NA
    fitted.readings$logistic.nls.r<-NA
    fitted.readings$logistic.nls.r.lower<-NA
    fitted.readings$logistic.nls.r.upper<-NA
    fitted.readings$logistic.nls.predicted<-NA
    
    
  }else{
    
    #extract the parameters from the model
    parameters<-coef(culture.model)    
    fitted.readings$logistic.nls.N0<-parameters["N0"]
    fitted.readings$logistic.nls.K<-parameters["K"]
    fitted.readings$logistic.nls.r<-parameters["r"]
    fitted.readings$logistic.nls.predicted<-predict(culture.model)
    
    
    #add confidence interval
    interval<-try(confint(culture.model))   
    if(any(class(interval)=="try-error", !do.interval)){
      fitted.readings$logistic.nls.K.lower<-NA
      fitted.readings$logistic.nls.K.upper<-NA
      fitted.readings$logistic.nls.r.lower<-NA
      fitted.readings$logistic.nls.r.upper<-NA}else{
      fitted.readings$logistic.nls.r.lower<-interval["r", "2.5%"]
      fitted.readings$logistic.nls.r.upper<-interval["r", "97.5%"]
      fitted.readings$logistic.nls.K.lower<-interval["K", "2.5%"]
      fitted.readings$logistic.nls.K.upper<-interval["K", "97.5%"]}
    
  }
  
  return(fitted.readings)
}

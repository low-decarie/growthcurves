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


logistic_growth_nls<-function(time, abundance, do.interval=F){

  start.values=c("K"=max(abundance, na.rm=T), "r"=1, "N0"=min(abundance, na.rm=T))
  
  
  culture.model <- try(nls(formula=abundance ~(K*N0*exp(r*time) ) / (K + N0 * (exp(r*time)-1)),
                           start=start.values,
                           na.action=na.exclude,
                           algorithm="port",
                           lower=c(0,0,0),
                           upper=c(2*max(abundance, na.rm=T),Inf, Inf)))
  
  fitted.readings <- data.frame(time,abundance)
  
  if(class(culture.model)=="try-error"){
    
    fitted.readings$N0<-NA
    fitted.readings$K<-NA
    fitted.readings$K.lower<-NA
    fitted.readings$K.upper<-NA
    fitted.readings$r<-NA
    fitted.readings$r.lower<-NA
    fitted.readings$r.upper<-NA
    fitted.readings$predicted<-NA
    
    
  }else{
    
    #extract the parameters from the model
    parameters<-coef(culture.model)    
    fitted.readings$N0<-parameters["N0"]
    fitted.readings$K<-parameters["K"]
    fitted.readings$r<-parameters["r"]
    fitted.readings$predicted<-predict(culture.model)
    
    
    #add confidence interval
    interval<-try(confint(culture.model))   
    if(all(class(interval)=="try-error", do.interval)){
      fitted.readings$K.lower<-NA
      fitted.readings$K.upper<-NA
      fitted.readings$r.lower<-NA
      fitted.readings$r.upper<-NA}
    
    if(all(!class(interval)=="try-error", do.interval)){
      fitted.readings$r.lower<-interval["r", "2.5%"]
      fitted.readings$r.upper<-interval["r", "97.5%"]
      fitted.readings$K.lower<-interval["K", "2.5%"]
      fitted.readings$K.upper<-interval["K", "97.5%"]}
    
  }
  
  return(fitted.readings)
}

#' exponential_decline_nls
#'
#' Least-squares estimates of the parameters of a exponential growth curve
#'
#' @param time Numeric vector of time since start of measurement (currently does not use date format)
#' @param abundance Abundance (eg. Optical density, counts etc).
#' @param do.interval Logical: should a 95\% confidence interval be produced
#' @keywords growth
#' @examples
#' data(navicula_growth)
#' fitted.readings <- with(navicula_growth,
#'                        exponential_growth_nls(time=Time,
#'                                            abundance=ABS))
#' head(fitted.readings)
#' @export

exponential_decline_nls<-function(time, abundance, do.interval=F){
  
  
  model<-try(nls(formula =abundance ~ abundance0 * 2 ^ (growth.rate * time),
                 start = list(abundance0 = max(abundance,na.rm=T),
                              growth.rate = -1),
                 na.action=na.exclude,
                 algorithm="port",
                 lower=c(0,-Inf),
                 upper=c(Inf,0)))
  
  fitted.readings <- data.frame(time,abundance)
  
  if(class(model)=="try-error"){
    
    fitted.readings$N0<-NA
    fitted.readings$r<-NA
    fitted.readings$predicted<-NA
    
  }else{
    
    parameters<-t(data.frame(coef(model,matrix=T)))
    
    fitted.readings$N0<-parameters[1]
    fitted.readings$r<-parameters[2]
    fitted.readings$predicted<-predict(model)
    
    interval<-try(confint(model)) 
    if(all(class(interval)=="try-error", do.interval)){
      fitted.readings$r.lower<-NA
      fitted.readings$r.upper<-NA
    }
    if(all(!class(interval)=="try-error", do.interval)){
      fitted.readings$r.lower<-interval["growth.rate", "2.5%"]
      fitted.readings$r.upper<-interval["growth.rate", "97.5%"]
    }
    
  }
  
  return(fitted.readings)
}

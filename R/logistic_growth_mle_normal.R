#' logistic.growth.mle.norm
#'
#' Maximum likelihood estimation of logistic growth parameters assuming normal error
#'
#' @param time Numeric vector of time since start of measurement (currently does not use date format)
#' @param abundance Abundance (eg. Optical density, counts etc).
#' @keywords growth
#' @examples
#' data(navicula_growth)
#' fitted.readings <- with(navicula_growth,
#'                        logistic_growth_mle_norm(time=Time,
#'                                            abundance=ABS))
#' head(fitted.readings)
#' @export

logistic_growth_mle_norm<-function(abundance, time, upper=2*max(abundance)){
  
  fitted.readings<-data.frame(time,abundance)
  
  r_start <- try(as.numeric(coef(lm(log(abundance)~time))[2]), silent = T)
  if(class(r_start)=="try-error")r_start <- 0.5
  
  start.values = c(K = max(abundance, na.rm = T),
                   r = r_start, 
                   N0 = max(c(min(abundance, na.rm = T),10^-3)),
                   st.dev = 0.5*sd(abundance))
  
  #Log likelyhood function to be minimized
  like.growth<-function(parameters=start.values, time, abundance){
    
    #Parameter extraction
    K<-parameters[[1]]
    r<-parameters[[2]]
    N0<-parameters[[3]]
    st.dev<-parameters[[4]]
    
    
    #Logistic growth model
    Nt<-(K*N0*exp(r*time) ) / (K + N0 * (exp(r*time)-1))
    #Nt<-(N0*K) / (N0 + (K-N0)*exp(-r*t))  #Synonymous model
    
    #log likelihood estimate
    #Nomral distribution
    likelihood<- -sum(dnorm(abundance, Nt, sd=st.dev, log=T))
    
    # Sanity bounds (remove if using "L-BFGS-B" or constrOptim
    if(any(c(Nt<0,
             Nt>upper, 
             K>upper,
             N0<0,
             r<0,
             st.dev<0,
             st.dev>upper))){likelihood<-NA}
    
    
    return(likelihood)
    
  }
  
  try.test<-try({
    
    fit<-optim(par=start.values,
               fn=like.growth,
               time=time,
               abundance=abundance)
    
    
    #extract fit values
    K<-fit$par[1]
    r<-fit$par[2]
    N0<-fit$par[3]
    
    predicted<-(K*N0*exp(r*time) ) / (K + N0 * (exp(r*time)-1))
    
    fitted.readings$N0<-N0
    fitted.readings$K<-K
    fitted.readings$r<-r
    fitted.readings$predicted<-predicted
  })
  
  #Pad with NAs for failed fits                                                             
  if(class(try.test)=="try-error"){
    fitted.readings$N0<-NA
    fitted.readings$K<-NA
    fitted.readings$r<-NA
    fitted.readings$predicted<-NA
  }
  
  return(fitted.readings)                                                     
  
}


# fn_prior_sig1<-function(a_sig, b_sig, sigma1){
#   #define a function,return the minus of the log likelihood
#   sigma1_sq <- sigma1^2
#   return((a_sig+1)*log(sigma1_sq)+b_sig/sigma1_sq-log(sigma1))
#   
# }
# fn_prior_sig2<-function(a_sig, b_sig, sigma2){
#   #define a function,return the minus of the log likelihood
#   sigma2_sq <- sigma2^2
#   return((a_sig+1)*log(sigma2_sq)+b_sig/sigma2_sq-log(sigma2))
#   
# }
# 
# fn_prior_beta01<-function(beta01,sig_beta01){
#   #define a function,return the minus of the log likelihood
#   return(beta01^2/sig_beta01^2)
# }
# fn_prior_beta11<-function(beta11,sig_beta11){
#   #define a function,return the minus of the log likelihood
#   return(beta11^2/sig_beta11^2)
# }
# fn_prior_beta02<-function(beta02,sig_beta02){
#   #define a function,return the minus of the log likelihood
#   return(beta02^2/sig_beta02^2)
# }
# fn_prior_beta12<-function(beta12,sig_beta12){
#   #define a function,return the minus of the log likelihood
#   return(beta12^2/sig_beta12^2)
# }

fn_logprior <- function(paraid, sampleid, paralist, sig_beta, a_sig, b_sig) {
  paraname <- paste(paraid, sampleid, sep = "")
  paraval <- paralist[[paraname]]
  if(paraid == "b0" | paraid == "b1") {
    s <- sig_beta[[paraname]]
    return(-paraval^2 / s^2)
  } else if(paraid == "sigma") {
    return(-( (2*a_sig + 1)*log(paraval) + b_sig/(paraval^2) ))
  } else if(paraid == "xi") {
    return(0)
  } else stop("Error in fn_logprior: Invalid paraid entered")
}
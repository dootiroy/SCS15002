library(msm)
### Metropolis Hastings update functions
## paraid: string identifying which parameter to update
#         valid values: "sigma","xi","b0","b1"
## sampleid: integer, valid values 1 or 2
## data.list: list(X, y1, y2, phi, delta1, delta2, n) 
## current: list of current values of parameters
#           list(sigma1, xi1, sigma2, xi2, b01, b11, b02, b12)

## Since all proposal densities are truncated normal, we need to get upper
#  and lower bounds, which are different for each parameter.
get_bounds <- function(paraid, sampleid, current, data.list) {
  
  # subset the data and current parameter by sampleid
  if(sampleid == 1) {
    y <- data.list$y1
    sigma <- current$sigma1
    xi <- current$xi1
    mu <- current$mu1
    b0 <- current$b01
    b1 <- current$b11
  } else if(sampleid == 2) {
    y <- data.list$y2
    sigma <- current$sigma2
    xi <- current$xi2
    mu <- current$mu2
    b0 <- current$b02
    b1 <- current$b12
  } else stop("Error in get_UB: Invalid sampleid entered")
  X <- data.list$X
  
  # Get the bound based on the paraid
  if(paraid == "sigma") {
    UB <- Inf
    nu <- xi * ( log(y) - mu )
    LB <- max(0, max(-nu))
  } else if(paraid == "xi") {
    eta <- log(y) - mu
    LB <- max( -1, max(-sigma / eta[eta>0]) )
    UB <- min(  1, min(-sigma / eta[eta<0]) )
  } else if(paraid == "b0") {
    a <- log(y) + (sigma/xi) - (b1*X)
    if(xi > 0) {
      UB <- min(a)
      LB <- -Inf
    } else {
      LB <- max(a)
      UB <- Inf
    }
  } else if(paraid == "b1") {
    b <- (log(y) + (sigma/xi) - b0) / X
    if(xi > 0) {
      UB <- min(b)
      LB <- -Inf
    } else {
      LB <- max(b)
      UB <- Inf
    }
  } else stop("Error in get_UB: Invalid paraid entered")
  l <- list(LB=LB, UB=UB)
}

# update a single parameter
# propsd: list of proposed sd's for rtnorm for all parameters
# sig_beta: list of length 4 of sd's for priors of beta.
# a_sig and b_sig are hyperparams for prior for sigma1 $ sigma2, common to both
update_para <- function(paraid, sampleid, current, data.list,
                        propsd, sig_beta, a_sig, b_sig, acc) {

  # Get upper and lower bounds for rtnorm proposal density (constraints)
  bounds <- get_bounds(paraid, sampleid, current, data.list)
  UB <- bounds$UB; LB <- bounds$LB
  
  # Generate proposed value based on proposed density and bounds and sd
  paraname <- paste(paraid, sampleid, sep="")
  curr.para <- current[[paraname]]
  prop.sd <- propsd[[paraname]]
  prop.para <- rtnorm(1, curr.para, prop.sd, lower = LB, upper = UB)
  # Make a parameter list proposed, which differs from current only
  # in the parameter of interest
  proposed <- current
  proposed[[paraname]] <- prop.para
  
  # Calc log(a)
  loga <- L_wrap3(proposed, data.list$X, data.list$y1, data.list$y2, data.list$phi, data.list$delta1,
                  data.list$delta2, data.list$n) - 
          L_wrap3(current, data.list$X, data.list$y1, data.list$y2, data.list$phi, data.list$delta1,
                  data.list$delta2, data.list$n) +
          fn_logprior(paraid, sampleid, proposed, sig_beta, a_sig, b_sig) -
          fn_logprior(paraid, sampleid, current, sig_beta, a_sig, b_sig) +
          dtnorm(curr.para, prop.para, prop.sd, lower = LB, upper = UB) - 
          dtnorm(prop.para, curr.para, prop.sd, lower = LB, upper = UB)
#   ### DEBUGGING #####################################################
#   print(paraname)
#   print(loga)
#   print(current)
#   print(proposed)
#   print(c(LB,UB))
#   print(L_wrap3(proposed, data.list$X, data.list$y1, data.list$y2, data.list$phi, data.list$delta1,
#                 data.list$delta2, data.list$n))
#     if(is.nan(loga)) {
#       save(list = c("current", "proposed", "LB", "UB", "paraname", "data.list"), 
#            file = "DEBUG/check.RData")
#       stop("DEBUG")
#     }
#   ####################################################################
  u <- runif(1)  
  u2 <- log(u) 
  
  # Update
  if(u2 < loga & is.finite(loga)) { 
    current <- proposed 
    acc[[paraname]] <- acc[[paraname]] + 1 
  }
  l <- list(current=current, acc=acc)
}

do_MCMC <- function(m, para0, data.list, propsd, sig_beta, a_sig, b_sig) {
  
  # Initialize counter values to calc acceptance prob
  acc <- list()
  acc$sigma1 <- 0
  acc$sigma2 <- 0
  acc$xi1 <- 0
  acc$xi2 <- 0
  acc$b01 <- 0
  acc$b11 <- 0
  acc$b02 <- 0
  acc$b12 <- 0
  
  # Initialize current parameter
  current <- para0
  current$mu1 <-  para0$b01 + para0$b11*data.list$X
  current$mu2 <-  para0$b02 + para0$b12*data.list$X
  
  # para is the matrix containing MCMC samples for all parameters
  # using updatelist with current ensures mu's dont get wiped out in updation step
  para <- matrix(nrow = m, ncol = length(para0)) 
  updatelist <- paste(c("sigma", "xi", "b0", "b1"), rep(1:2, each=4), sep="")
  paraids <- c("sigma", "xi", "b0", "b1")
  sampleids <- 1:2
  
  # MCMC
  for(i in 1:m) {
    for(j in seq_along(paraids)) {
      for(k in seq_along(sampleids)) {
        ### DEBUG ####################################
        #print(c(i,j,k))
        #######################################
        l <- update_para(paraids[j], sampleids[k], current, data.list, 
                         propsd, sig_beta, a_sig, b_sig, acc)
        current[updatelist] <- l$current[updatelist]
        acc <- l$acc
        # update mu's if necessary
        if(paraids[j] == "b0" | paraids[j] == "b1") {
          if(sampleids[k] == 1) {
            current$mu1 <- current$b01 + current$b11*data.list$X
          } else {
            current$mu2 <- current$b02 + current$b12*data.list$X
          }
        }
      }
    }
    para[i, ] <- unlist(current[updatelist])
  }
  list(para=para, acc=acc)
}
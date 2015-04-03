B <-function(phi,phi1,PARA,y1,y2,delta1,delta2,n)
{
  s <- 0
  N = nrow(PARA)
  for(i in 1:N)
  {
    num = exp(L(PARA[i,3] + PARA[i,4]*X,PARA[i,1],PARA[i,2],PARA[i,7] + 
                  PARA[i,8]*X,PARA[i,5],PARA[i,6],y1,y2,phi,delta1,delta2,n) - 
                L(PARA[i,3] + PARA[i,4]*X,PARA[i,1],PARA[i,2],PARA[i,7] + 
                    PARA[i,8]*X,PARA[i,5],PARA[i,6],y1,y2,phi1,delta1,delta2,n))
    s = s + num
  }
  ratio = s/N
  return(ratio)
}

# Function to maximize to estimate eta
# H is list of matrices
etafn <- function(eta_reduced, H) {
  eta <- c(1, eta_reduced) #eta_1 = 1
  mysum <- 0
  for(j in 1:length(H)) {
    H[[j]] <- t( t(H[[j]]) + eta ) # Add eta_k to kth col of H[[j]] for all j
    H[[j]] <- H[[j]] - H[[j]][,j] # Subtract jth col of H[[j]] from each col of H[[j]] for all j
    H[[j]] <- exp(H[[j]]) # Take exponent of entire matrix
    mysum <- mysum + sum( log( rowSums(H[[j]]) ) ) # take rowSums (sum since its a vector) and sum over j
  }
  -mysum
}

# Function to estimate final phi
B_fin <- function(phi, H, eta, dlist, PARA, skel.points) {
  # H is list of matrices 
  # eta is estimated eta with eta1 = 1
  # dlist is data.list
  # PARA is list of PARA1-4
  mysum <- 0
  N <- nrow(PARA[[1]])
  h_phi <- numeric(N)
  for(j in 1:length(H)) {
    # create the vector of likelihoods using phi
    for(i in 1:N) {
      h_phi[i] <- L(PARA[[j]][i,3] + PARA[[j]][i,4]*X,PARA[[j]][i,1],
                    PARA[[j]][i,2],PARA[[j]][i,7] + PARA[[j]][i,8]*X,
                    PARA[[j]][i,5],PARA[[j]][i,6],y1,y2,phi,
                    delta1,delta2,n)
    }
    H[[j]] <- H[[j]] - h_phi 
    H[[j]] <- t( t( exp(H[[j]]) )/eta ) # Take exponent of entire matrix, divide kth column by eta[k]
    H[[j]] <- H[[j]]*N
    
    mysum <- mysum + sum( 1/rowSums(H[[j]])  ) # take rowSums (sum since its a vector) and sum over j
  }
  mysum
}
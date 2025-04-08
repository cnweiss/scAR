#Fitting ordinal soft-clipping models models to ordinal fear states data

cReLU <- function(x) min(1, max(0,x))
cReLU <- Vectorize(cReLU)

sc <- function(x, pc, tol=0){
  min(1-tol, max(0+tol,x)) + pc*log1p(exp(-abs(x/pc))) - pc*log1p(exp(-abs((1-x)/pc)))
}
sc <- Vectorize(sc)

sci <- function(y, pc, tol=0){
  y + pc*log1p(-exp(-y/pc)) - pc*log1p(-exp(-(1-y)/pc))
}
sci <- Vectorize(sci)




#Ordinal binarization
ordbin <- function(j,d) 1 * (j <= (0:(d-1)))





#Conditional log-likelihood function for soft-clipping AR(p) process (par=c(parmarg,pardep) with parmarg=c(eta0,...,eta_{d-1})).
#Model version "IdD-scAR" with scalar parameters:

ll_scarp_scalar <- function(par,datanum,d, pc){ #states 1:(d+1)
  T <- length(datanum)
  p <- length(par)-d
  eta <- par[1:d]
  pardep <- par[(d+1):(d+p)] #alpha1,...,alphap
  
  value <- 0
  for (t in c((p+1):T)) {
    ft01 <- eta
    for(k in 1:p) ft01 <- ft01 + pardep[k]*obincodes[datanum[t-k],]
    ft01 <- c(0, sc(ft01, pc), 1)
    
    value <- value - log(ft01[datanum[t]+1]-ft01[datanum[t]])
  }
  value
}


#Constraints for optimization A %*% par >= b:
#etas in [0,1], non-decreasing
Amat_scalar <- function(d,p){
  value <- array(0,c(d+p+1, d+p))
  
  value[1,1] <- 1 #first CDF >=0
  for(i in 2:d){ #differences >=0
    value[i,i] <- 1
    value[i,i-1] <- (-1)
  }
  value[d+1,d] <- (-1) #last CDF <=1

  for(i in 1:p){ #AR >=0
    value[d+1+i,d+i] <- 1
  }
  
  value
}

bvec_scalar <- function(d,p){
  value <- rep(0, d+1+p)
  value[d+1] <- (-1)
  value
}






#Conditional log-likelihood function for soft-clipping AR(1,1) process (par=c(parmarg,pardep) with parmarg=c(eta0,...,eta_{d-1})).
#Model version "IdD-scAR" with scalar parameters:

ll_scar11_scalar <- function(par,datanum,d, pc){ #states 1:(d+1)
  T <- length(datanum)
  eta <- par[1:d]
  pardep <- par[(d+1):(d+2)] #alpha1,beta1
  
  #Initialize feedback term:
  feedb <- sc(eta, pc)
  
  value <- 0
  for (t in c(2:T)) {
    feedb <- sc(eta + pardep[1]*obincodes[datanum[t-1],] + pardep[2]*feedb, pc)
    ft01 <- c(0, feedb, 1)
    
    value <- value - log(ft01[datanum[t]+1]-ft01[datanum[t]])
  }
  value
}


#Constraints for optimization A %*% par >= b:
#use Amat_scalar and bvec_scalar with p=2






#Conditional log-likelihood function for soft-clipping AR(p,1) process (par=c(parmarg,pardep) with parmarg=c(eta0,...,eta_{d-1})).
#Model version "IdD-scAR" with scalar parameters:

ll_scarp1_scalar <- function(par,datanum,d, pc){ #states 1:(d+1)
  T <- length(datanum)
  p <- length(par)-d-1
  eta <- par[1:d]
  pardep <- par[-(1:d)] #alpha1,...,alphap,beta1
  
  #Initialize feedback term:
  feedb <- sc(eta, pc)
  # feedb <- hatf
  
  value <- 0
  for (t in c((p+1):T)) {
    ft01 <- eta + pardep[p+1]*feedb
	for(k in 1:p) ft01 <- ft01 + pardep[k]*obincodes[datanum[t-k],]
    feedb <- sc(ft01, pc)
    ft01 <- c(0, feedb, 1)
    
    value <- value - log(ft01[datanum[t]+1]-ft01[datanum[t]])
  }
  value
}


#Constraints for optimization A %*% par >= b:
#use Amat_scalar and bvec_scalar with p+1






#Conditional log-likelihood function for soft-clipping AR(p,2) process (par=c(parmarg,pardep) with parmarg=c(eta0,...,eta_{d-1})).
#Model version "IdD-scAR" with scalar parameters:

ll_scarp2_scalar <- function(par,datanum,d, pc){ #states 1:(d+1)
  T <- length(datanum)
  p <- length(par)-d-2
  eta <- par[1:d]
  pardep <- par[-(1:d)] #alpha1,...,alphap,beta1,beta2
  
  #Initialize feedback term:
  feedb1 <- feedb2 <- sc(eta, pc)
  
  value <- 0
  for (t in c((p+1):T)) {
    ft01 <- eta + pardep[p+1]*feedb1 + pardep[p+2]*feedb2
	for(k in 1:p) ft01 <- ft01 + pardep[k]*obincodes[datanum[t-k],]
	feedb2 <- feedb1
    feedb1 <- sc(ft01, pc)
    ft01 <- c(0, feedb1, 1)
    
    value <- value - log(ft01[datanum[t]+1]-ft01[datanum[t]])
  }
  value
}


#Constraints for optimization A %*% par >= b:
#use Amat_scalar and bvec_scalar with p+2






#Conditional log-likelihood function for soft-clipping AR(p,3) process (par=c(parmarg,pardep) with parmarg=c(eta0,...,eta_{d-1})).
#Model version "IdD-scAR" with scalar parameters:

ll_scarp3_scalar <- function(par,datanum,d, pc){ #states 1:(d+1)
  T <- length(datanum)
  p <- length(par)-d-3
  eta <- par[1:d]
  pardep <- par[-(1:d)] #alpha1,...,alphap,beta1,beta2,beta3
  
  #Initialize feedback term:
  feedb1 <- feedb2 <- feedb3 <- sc(eta, pc)
  
  value <- 0
  for (t in c((p+1):T)) {
    ft01 <- eta + pardep[p+1]*feedb1 + pardep[p+2]*feedb2 + pardep[p+3]*feedb3
	for(k in 1:p) ft01 <- ft01 + pardep[k]*obincodes[datanum[t-k],]
	feedb3 <- feedb2
	feedb2 <- feedb1
    feedb1 <- sc(ft01, pc)
    ft01 <- c(0, feedb1, 1)
    
    value <- value - log(ft01[datanum[t]+1]-ft01[datanum[t]])
  }
  value
}


#Constraints for optimization A %*% par >= b:
#use Amat_scalar and bvec_scalar with p+3






#Conditional log-likelihood function for soft-clipping AR(p,4) process (par=c(parmarg,pardep) with parmarg=c(eta0,...,eta_{d-1})).
#Model version "IdD-scAR" with scalar parameters:

ll_scarp4_scalar <- function(par,datanum,d, pc){ #states 1:(d+1)
  T <- length(datanum)
  p <- length(par)-d-4
  eta <- par[1:d]
  pardep <- par[-(1:d)] #alpha1,...,alphap,beta1,...,beta4
  
  #Initialize feedback term:
  feedb1 <- feedb2 <- feedb3 <- feedb4 <- sc(eta, pc)
  
  value <- 0
  for (t in c((p+1):T)) {
    ft01 <- eta + pardep[p+1]*feedb1 + pardep[p+2]*feedb2 + pardep[p+3]*feedb3 + pardep[p+4]*feedb4
	for(k in 1:p) ft01 <- ft01 + pardep[k]*obincodes[datanum[t-k],]
	feedb4 <- feedb3
	feedb3 <- feedb2
	feedb2 <- feedb1
    feedb1 <- sc(ft01, pc)
    ft01 <- c(0, feedb1, 1)
    
    value <- value - log(ft01[datanum[t]+1]-ft01[datanum[t]])
  }
  value
}


#Constraints for optimization A %*% par >= b:
#use Amat_scalar and bvec_scalar with p+4






#Fear states based on closing of VIX index for period 1990-2006

#End of years:
years <- c(1990:2007)
ends <- c(253,506,760,1013,1265,1517,1771,2024,2276,2528,2780,3028,3280,3532,3784,4036,4287)

#Fear states:
#"VLA" = "at most very low anxiety"
#"LA" = "low anxiety"
#"MA" = "moderate anxiety"
#"MHA" = "at least moderately high anxiety"

states <- c("VLA","LA","MA","MHA") #position in states corresponds to numeric code
nostates <- length(states)
nostates #4

d <- nostates-1

#Binarization:
nbincodes <- diag(1,nostates) #nominal binarization
obincodes <- array(0, c(nostates,d))
for(s in 0:d) obincodes[s+1,] <- ordbin(s,d)


datanum <- c(2,2,2,3,3,3,3,3,3,4,3,3,3,3,4,3,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,2,2,3,3,3,2,2,2,2,2,3,3,3,2,3,2,2,2,3,3,3,3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,3,4,4,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,3,2,2,1,1,1,2,2,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,3,2,3,3,2,3,3,2,2,3,3,2,2,3,2,2,2,3,2,2,2,2,2,2,
2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,1,1,2,2,2,2,2,1,1,1,1,1,1,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,2,3,2,3,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,1,2,1,1,1,1,1,1,1,1,2,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,2,2,2,2,1,2,1,2,1,2,1,1,1,1,1,1,1,1,1,1,2,2,3,3,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,1,1,1,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,2,2,2,2,2,2,3,3,2,2,2,3,3,3,2,2,3,3,2,2,2,2,2,2,2,2,2,2,2,2,1,2,1,1,1,1,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,1,1,1,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,2,2,2,2,2,2,2,2,3,3,2,2,2,3,3,2,2,2,2,2,2,2,2,2,2,2,3,3,3,2,2,2,2,3,3,2,3,2,2,2,2,2,3,3,3,2,2,3,3,3,3,3,2,3,2,2,2,2,2,2,3,3,3,3,2,3,2,2,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,3,2,3,3,3,3,2,3,2,2,3,3,3,3,3,3,3,3,2,3,3,2,2,2,2,2,2,2,2,3,3,3,3,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,3,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,3,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,3,2,3,3,3,3,2,2,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
4,4,4,4,4,4,4,4,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,3,2,2,2,2,2,2,2,2,2,2,2,3,2,3,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,3,3,4,3,3,3,2,3,3,3,3,3,3,3,3,2,3,3,3,2,2,2,3,3,3,3,3,3,3,3,2,3,3,3,3,3,4,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,4,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,4,3,4,4,4,3,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,4,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,3,3,3,4,4,4,3,3,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,3,3,3,3,3,3,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,3,3,3,3,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,2,2,3,3,3,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
3,3,3,2,2,2,3,3,3,3,3,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,3,3,3,3,3,2,3,3,3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,3,3,4,3,3,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,3,3,2,3,3,3,3,2,2,2,2,2,2,2,2,2,2,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,3,2,2,3,3,3,3,3,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,3,3,3,4,4,4,4,3,3,3,3,3,3,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,3,4,4,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,4,4,4,3,4,3,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,3,2,2,3,3,3,3,3,3,3,3,2,3,3,3,3,3,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,3,3,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,3,3,3,4,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
4,4,4,4,4,3,3,3,3,3,3,4,4,3,4,3,3,3,3,4,4,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,3,3,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,3,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,3,2,2,2,3,3,3,2,3,3,3,3,3,3,3,3,4,4,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,4,4,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,3,2,2,3,3,3,2,2,2,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,3,3,3,2,2,2,2,2,2,3,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,2,2,2,3,2,3,3,2,2,2,3,3,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,1,2,2,2,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,1,2,1,1,1,2,2,2,1,2,2,2,2,2,2,1,1,1,1,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,2,2,2,1,2,1,1,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,1,2,1,2,1,1,1,2,2,2,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,1,1,2,1,2,2,1,1,1,2,
1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,1,2,2,1,1,2,2,2,2,2,3,3,3,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,1,1,1,1,1,1,2,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

Tlen <- length(datanum)
Tlen #4287

databin <- nbincodes[datanum,] #assign row vectors
dataobin <- obincodes[datanum,] #assign row vectors





#####################
#Descriptive analysis
#####################


#relative frequencies
hatpi <- colMeans(databin)
hatpi
#0.3272685 0.2857476 0.2192676 0.1677164

hatf <- cumsum(hatpi)[1:d]
hatf
#0.3272685 0.6130161 0.8322836

#Median:
min((0:(d-1))[hatf>=0.5]) #1


#IOV:
IOV <- 4/d * sum(hatf*(1-hatf))
IOV #0.7959717

#CPE:
CPE <- -1/d/log(2) * sum(hatf*log(hatf) + (1-hatf)*log(1-hatf))
CPE #0.8424636

#skew:
skew <- 2/d * sum(hatf) - 1
skew #0.1817122


#Time series plot of data
plot(datanum, type="b", pch=19, cex=0.5, ylim=c(0.5,nostates+.5), xlab="Trading days", ylab="Fear state", yaxt="n", xaxt="n")
axis(side=1, at=c(1, ends+1), labels=years, cex.axis=0.85, las=2)
axis(side=2, at=c(1:nostates), labels=states, las=2)







maxlag <- 15
hatbivprob <- array(0,c(nostates,nostates,maxlag))

for(k in c(1:maxlag)){
  for(i in c(1:nostates)){
    for(j in c(1:nostates)){
      hatbivprob[i,j,k] <- mean(databin[(k+1):Tlen,i]*databin[1:(Tlen-k),j])
    }
  }
}

#... which is to be compared with
indprob <- hatpi %*% t(hatpi)

#nominal kappa:
cohennom <- rep(0,maxlag)
for(k in c(1:maxlag)){
  cohennom[k] <- (sum(diag(hatbivprob[,,k]))-sum(hatpi^2))/(1-sum(hatpi^2))
}
cohennom
#0.8374792 0.7767993 0.7322872 0.7093542 0.6883168 0.6701297 0.6528878 0.6448581 0.6339625 0.6252885 0.6115197 0.6053823 0.5976504 0.5870494 0.5828127


#diagonal of bivariate cdf
hatbivcdf <- array(0, c(nostates-1,maxlag))
for(k in c(1:maxlag)){
  for(i in c(1:d)){
    hatbivcdf[i,k] <- sum(hatbivprob[1:i,1:i,k])
  }
}

#Required for computing asymptotics:
hatfmat <- array(0, c(nostates-1,nostates-1))
for(i in c(1:d)){
  for(j in c(1:d)){
    hatfmat[i,j] <- hatf[min(i,j)]-hatf[i]*hatf[j]
  }} #for i,j



#ordinal kappa:
cohenord <- rep(0,maxlag)
for(k in c(1:maxlag)){
  cohenord[k] <- (sum(hatbivcdf[,k])-sum(hatf^2))/sum(hatf*(1-hatf))
}
cohenord
#0.8994676 0.8618267 0.8327705 0.8169983 0.8047395 0.7920837 0.7809875 0.7749752 0.7677855 0.7621590 0.7526124 0.7477633 0.7413443 0.7337463 0.7304579


#Critical value
alpha <- 0.05
crit2 <- qnorm(1-alpha/2, sd=sqrt(16/d^2/IOV^2 * sum(hatfmat^2)/Tlen))
#0.0214786

#Plot of ordinal Cohen's kappa:
plot(cohenord, type="h", xlab = "h", ylab = expression(paste("Cohen's   ",kappa[o](h))), lwd=4, ylim=c(-0.5,1), lend=3)
abline(h=c(-1/Tlen+crit2,-1/Tlen-crit2), col=gray(0.5), lwd=2, lty=2)
abline(h=-1/Tlen)



#Compute "partial ordinal Cohen's kappa" recursively as in Appendix B.3:
cohenordpart <- rep(0,maxlag)
cohenordpart[1] <- cohenord[1]

a <- array(0, c(2,maxlag))
a[1,1] <- cohenord[1]

for(k in c(1:(maxlag-1))){
  a[2,k+1] <- (cohenord[k+1] - a[1,1:k] %*% cohenord[k:1]) / (1 - a[1,1:k] %*% cohenord[1:k])
  for(j in c(1:k)){
    a[2,j] <- a[1,j]-a[2,k+1]*a[1,k-j+1]
  }
  cohenordpart[k+1] <- a[2,k+1]
  a[1,] <- a[2,]
}
plot(cohenordpart, type="h", xlab = "h", ylab = expression(paste("Cohen's   ",kappa["o; part"](h))), lwd=4, ylim=c(-0.5,1), lend=3)
abline(h=c(crit2,-crit2), col=gray(0.5), lwd=2, lty=2)
abline(h=0)

cohenordpart
#0.89946756 0.27642069 0.13172014 0.12983104 0.09972284 0.06367958 0.05450205 0.07069206 0.04784418 0.04716901 0.01745251 0.04156909 0.02556315 0.01246609 0.03823501






#Fit IdD-scAR(1) model with scalar dependence parameters:
delta <- 0.01 #soft-clipping parameter

#Initialize scAR(1) estimates:
start <- c(hatf*(1-cohenord[1]), cohenord[1])

scar1.s <- suppressWarnings(constrOptim(start, ll_scarp_scalar, NULL, ui=Amat_scalar(d,1), ci=bvec_scalar(d,1), datanum=datanum, d=d, pc=delta))
scar1.s$par
#0.03026081 0.09081260 0.14217978 0.84136280
scar1.s$value
#2082.806

#AIC and BIC:
AIC.scar1.s <- 2*Tlen/(Tlen-1)*scar1.s$value+2*length(scar1.s$par)
BIC.scar1.s <- 2*Tlen/(Tlen-1)*scar1.s$value+log(Tlen)*length(scar1.s$par)
c(AIC.scar1.s, BIC.scar1.s) #4174.583 4200.037









#Fit IdD-scAR(2) model with scalar dependence parameters:
delta <- 0.01 #soft-clipping parameter

#Initialize scAR(2) estimates:
#"moment estimates"

al1MM <- cohenord[1]*(1-cohenord[2])/(1-cohenord[1]^2) 
al2MM <- (cohenord[2]-cohenord[1]^2)/(1-cohenord[1]^2) 
c(al1MM,al2MM) #0.6508361 0.2764207
start <- c(hatf*(1-al1MM-al2MM), al1MM,al2MM)

scar2.s <- suppressWarnings(constrOptim(start, ll_scarp_scalar, NULL, ui=Amat_scalar(d,2), ci=bvec_scalar(d,2), datanum=datanum, d=d, pc=delta))
scar2.s$par
#0.01958477 0.06308512 0.10171479 0.62178283 0.26464989
scar2.s$value
#1875.293


#AIC and BIC:
AIC.scar2.s <- 2*Tlen/(Tlen-2)*scar2.s$value+2*length(scar2.s$par)
BIC.scar2.s <- 2*Tlen/(Tlen-2)*scar2.s$value+log(Tlen)*length(scar2.s$par)
c(AIC.scar2.s, BIC.scar2.s) #3762.337 3794.154





#Fit IdD-scAR(1,1) model with scalar dependence parameters:
delta <- 0.01 #soft-clipping parameter

#Initialize scAR(1,1) estimates:
start <- c(hatf*(1-cohenord[1]), cohenord[1]/2, cohenord[1]/2)

scar11.s <- suppressWarnings(constrOptim(start, ll_scar11_scalar, NULL, ui=Amat_scalar(d,2), ci=bvec_scalar(d,2), datanum=datanum, d=d, pc=delta))
scar11.s$par
#0.001116733 0.011708363 0.023134872 0.374006226 0.606489547
scar11.s$value
#1685.245


#AIC and BIC:
AIC.scar11.s <- 2*Tlen/(Tlen-1)*scar11.s$value+2*length(scar11.s$par)
BIC.scar11.s <- 2*Tlen/(Tlen-1)*scar11.s$value+log(Tlen)*length(scar11.s$par)
c(AIC.scar11.s, BIC.scar11.s) #3381.277 3413.094



# parameters:
etaML11 <- scar11.s$par[1:d]
alphaML11<-scar11.s$par[d+1]
betaML11<-scar11.s$par[d+2]


#PIT:
nobins <- 10
PIT <- array(0, c(2,nobins+1))

for(j in c(1:nobins)){
  u <- j/nobins
  pitval <- 0
  
  feedb <- sc(etaML11,pc=delta)
  for(t in c(2:Tlen)){
    feedb <-sc(etaML11 + alphaML11*obincodes[datanum[t-1],]+ betaML11*feedb, delta)
    ft <- c(0, feedb, 1)
    if(ft[datanum[t]]<u){
      if(ft[(datanum[t]+1)]<u){
        pitval <- pitval+1
      }else{
        pitval <- pitval+ (u-ft[datanum[t]])/(ft[(datanum[t]+1)]-ft[datanum[t]])
      }
    }
  }
  PIT[1,j+1] <- pitval/(Tlen-1)
  PIT[2,j+1] <- PIT[1,j+1]-PIT[1,j]
}

PIT.freq <- as.vector(rep(((1:nobins)-0.5)/nobins, PIT[2,2:(nobins+1)]*1000))
PIT.hist.tpar <- hist(PIT.freq, plot=FALSE, breaks=nobins)
PIT.hist.tpar$density <- PIT.hist.tpar$counts/1000

plot(PIT.hist.tpar, freq=FALSE, main="fitted IdD-scAR(1,1) model", ylab="PIT histogram", xlab="u", col="gray", ylim=c(0,0.15), yaxt="n")
axis(side=2)
abline(h=0.1, lwd=2)



#Marginal calibration:
tabf <- array(0, c(Tlen-1, d))

feedb <- sc(etaML11,pc=delta)
for(t in c(2:Tlen)){
    feedb <-sc(etaML11 + alphaML11*obincodes[datanum[t-1],]+ betaML11*feedb, delta)
	
	tabf[t-1, ] <- feedb
}

#Marginal calibration:
hatf
#0.3272685 0.6130161 0.8322836
colMeans(tabf)
#0.3188074 0.6108630 0.8367334





#Fit IdD-scAR(1,2) model with scalar dependence parameters:
delta <- 0.01 #soft-clipping parameter

#Initialize with scAR(1,1) estimates:
start <- c(0.001, 0.012, 0.023, 0.374, 0.606, 0.01)

scar12.s <- suppressWarnings(constrOptim(start, ll_scarp2_scalar, NULL, ui=Amat_scalar(d,3), ci=bvec_scalar(d,3), datanum=datanum, d=d, pc=delta))
scar12.s$par
#0.00105265 0.01284187 0.02544572 0.44828579 0.26007047 0.27065205
scar12.s$value
#1656.729

#AIC and BIC:
AIC.scar12.s <- 2*Tlen/(Tlen-1)*scar12.s$value+2*length(scar12.s$par)
BIC.scar12.s <- 2*Tlen/(Tlen-1)*scar12.s$value+log(Tlen)*length(scar12.s$par)
c(AIC.scar12.s, BIC.scar12.s) #3326.232 3364.412





#Fit IdD-scAR(1,3) model with scalar dependence parameters:
delta <- 0.01 #soft-clipping parameter

#Initialize with scAR(1,2) estimates:
start <- c(0.001, 0.013, 0.025, 0.448, 0.260, 0.271, 0.01)

scar13.s <- suppressWarnings(constrOptim(start, ll_scarp3_scalar, NULL, ui=Amat_scalar(d,4), ci=bvec_scalar(d,4), datanum=datanum, d=d, pc=delta))
scar13.s$par
#0.000782904 0.012996868 0.025807186 0.482906312 0.266804353 0.023907697 0.205409577
scar13.s$value
#1636.611

#AIC and BIC:
AIC.scar13.s <- 2*Tlen/(Tlen-1)*scar13.s$value+2*length(scar13.s$par)
BIC.scar13.s <- 2*Tlen/(Tlen-1)*scar13.s$value+log(Tlen)*length(scar13.s$par)
c(AIC.scar13.s, BIC.scar13.s) #3287.985 3332.528





#Fit IdD-scAR(1,4) model with scalar dependence parameters:
delta <- 0.01 #soft-clipping parameter

#Initialize with scAR(1,3) estimates:
start <- c(0.001, 0.013, 0.026, 0.48, 0.260, 0.024, 0.205, 0.01)

scar14.s <- suppressWarnings(constrOptim(start, ll_scarp4_scalar, NULL, ui=Amat_scalar(d,5), ci=bvec_scalar(d,5), datanum=datanum, d=d, pc=delta))
scar14.s$par
#0.0007310307 0.0130068441 0.0258548665 0.5020720033 0.2436145520 0.0408705320 0.0732826738 0.1193761928
scar14.s$value
#1628.88

#AIC and BIC:
AIC.scar14.s <- 2*Tlen/(Tlen-1)*scar14.s$value+2*length(scar14.s$par)
BIC.scar14.s <- 2*Tlen/(Tlen-1)*scar14.s$value+log(Tlen)*length(scar14.s$par)
c(AIC.scar14.s, BIC.scar14.s) #3274.520 3325.427


# parameters:
etaML14 <- scar14.s$par[1:d]
alpha1ML14<-scar14.s$par[d+1]
beta1ML14<-scar14.s$par[d+2]
beta2ML14<-scar14.s$par[d+3]
beta3ML14<-scar14.s$par[d+4]
beta4ML14<-scar14.s$par[d+5]


#PIT:
nobins <- 10
PIT <- array(0, c(2,nobins+1))

for(j in c(1:nobins)){
  u <- j/nobins
  pitval <- 0
  
  feedb1 <- feedb2 <- feedb3 <- feedb4 <- sc(etaML14,pc=delta)
  for(t in c(2:Tlen)){
    ft01 <- etaML14 + alpha1ML14*obincodes[datanum[t-1],] + beta1ML14*feedb1 + beta2ML14*feedb2 + beta3ML14*feedb3 + beta4ML14*feedb4
	feedb4 <- feedb3
	feedb3 <- feedb2
	feedb2 <- feedb1
    feedb1 <- sc(ft01, pc=delta)
    ft <- c(0, feedb1, 1)
    if(ft[datanum[t]]<u){
      if(ft[(datanum[t]+1)]<u){
        pitval <- pitval+1
      }else{
        pitval <- pitval+ (u-ft[datanum[t]])/(ft[(datanum[t]+1)]-ft[datanum[t]])
      }
    }
  }
  PIT[1,j+1] <- pitval/(Tlen-1)
  PIT[2,j+1] <- PIT[1,j+1]-PIT[1,j]
}

PIT.freq <- as.vector(rep(((1:nobins)-0.5)/nobins, PIT[2,2:(nobins+1)]*1000))
PIT.hist.tpar <- hist(PIT.freq, plot=FALSE, breaks=nobins)
PIT.hist.tpar$density <- PIT.hist.tpar$counts/1000

plot(PIT.hist.tpar, freq=FALSE, main="fitted IdD-scAR(1,4) model", ylab="PIT histogram", xlab="u", col="gray", ylim=c(0,0.15), yaxt="n")
axis(side=2)
abline(h=0.1, lwd=2)



#Marginal calibration:
tabf <- array(0, c(Tlen-1, d))
  
feedb1 <- feedb2 <- feedb3 <- feedb4 <- sc(etaML14,pc=delta)
for(t in c(2:Tlen)){
    ft01 <- etaML14 + alpha1ML14*obincodes[datanum[t-1],] + beta1ML14*feedb1 + beta2ML14*feedb2 + beta3ML14*feedb3 + beta4ML14*feedb4
	feedb4 <- feedb3
	feedb3 <- feedb2
	feedb2 <- feedb1
    feedb1 <- sc(ft01, pc=delta)
	
	tabf[t-1, ] <- feedb1
}

#Marginal calibration:
hatf
#0.3272685 0.6130161 0.8322836
colMeans(tabf)
#0.3200150 0.6117295 0.8365475









#########################################
#Competitors from different model classes
#########################################






#Conditional log-likelihood function for logit AR(1,1) process (par=c(parmarg,pardep) with parmarg=c(eta0,...,eta_{d-1})).
#Model version with scalar parameters:

ll_logitar11_scalar <- function(par,datanum,d, pc){ #states 1:(d+1)
  T <- length(datanum)
  eta <- par[1:d]
  pardep <- par[(d+1):(d+2)] #alpha1,beta1
  
  #Initialize feedback term:
  feedb <- plogis(eta)
  
  value <- 0
  for (t in c(2:T)) {
    feedb <- plogis(eta + pardep[1]*obincodes[datanum[t-1],] + pardep[2]*feedb)
    ft01 <- c(0, feedb, 1)
    
    value <- value - log(ft01[datanum[t]+1]-ft01[datanum[t]])
  }
  value
}


#Constraints for optimization A %*% par >= b:
#etas non-decreasing
Amat_logit <- function(d,p){
  value <- array(0,c(d+p-1, d+p))
  
  for(i in 1:(d-1)){ #differences >=0
    value[i,i+1] <- 1
    value[i,i] <- (-1)
  }

  for(i in 1:p){ #AR >=0
    value[d-1+i,d+i] <- 1
  }
  
  value
}

bvec_logit <- function(d,p){
  value <- rep(0, d-1+p)
  value
}





#Fit logitAR(1,1) model with scalar dependence parameters:

#Initialize logitAR(1,1) estimates:
start <- c(qlogis(hatf*(1-cohenord[1])), 3,3)

logitar11.s <- suppressWarnings(constrOptim(start, ll_logitar11_scalar, NULL, ui=Amat_logit(d,2), ci=bvec_logit(d,2), datanum=datanum, d=d))
logitar11.s$par
#-3.627461 -3.081534 -2.521130  3.251302  3.140576
logitar11.s$value
#1837.129

#AIC and BIC:
AIC.logitar11.s <- 2*Tlen/(Tlen-1)*logitar11.s$value+2*length(logitar11.s$par)
BIC.logitar11.s <- 2*Tlen/(Tlen-1)*logitar11.s$value+log(Tlen)*length(logitar11.s$par)
c(AIC.logitar11.s, BIC.logitar11.s) #3685.115 3716.932





# parameters:
etaML11 <- logitar11.s$par[1:d]
alphaML11<-logitar11.s$par[d+1]
betaML11<-logitar11.s$par[d+2]


#PIT:
nobins <- 10
PIT <- array(0, c(2,nobins+1))

for(j in c(1:nobins)){
  u <- j/nobins
  pitval <- 0
  
  feedb <- plogis(etaML11)
  for(t in c(2:Tlen)){
    feedb <-plogis(etaML11 + alphaML11*obincodes[datanum[t-1],]+ betaML11*feedb)
    ft <- c(0, feedb, 1)
    if(ft[datanum[t]]<u){
      if(ft[(datanum[t]+1)]<u){
        pitval <- pitval+1
      }else{
        pitval <- pitval+ (u-ft[datanum[t]])/(ft[(datanum[t]+1)]-ft[datanum[t]])
      }
    }
  }
  PIT[1,j+1] <- pitval/(Tlen-1)
  PIT[2,j+1] <- PIT[1,j+1]-PIT[1,j]
}

PIT.freq <- as.vector(rep(((1:nobins)-0.5)/nobins, PIT[2,2:(nobins+1)]*1000))
PIT.hist.tpar <- hist(PIT.freq, plot=FALSE, breaks=nobins)
PIT.hist.tpar$density <- PIT.hist.tpar$counts/1000

plot(PIT.hist.tpar, freq=FALSE, main="fitted IdD-logitAR(1,1) model", ylab="PIT histogram", xlab="u", col="gray", ylim=c(0,0.15), yaxt="n")
axis(side=2)
abline(h=0.1, lwd=2)



#Marginal calibration:
tabf <- array(0, c(Tlen-1, d))

feedb <- plogis(etaML11)
for(t in c(2:Tlen)){
    feedb <- plogis(etaML11 + alphaML11*obincodes[datanum[t-1],]+ betaML11*feedb)
	
	tabf[t-1, ] <- feedb
}

#Marginal calibration:
hatf
#0.3272685 0.6130161 0.8322836
colMeans(tabf)
#0.3148705 0.6106575 0.8397323






#Conditional log-likelihood function for logitAR(p,4) process (par=c(parmarg,pardep) with parmarg=c(eta0,...,eta_{d-1})).
#Model version with scalar parameters:

ll_logitarp4_scalar <- function(par,datanum,d, pc){ #states 1:(d+1)
  T <- length(datanum)
  p <- length(par)-d-4
  eta <- par[1:d]
  pardep <- par[-(1:d)] #alpha1,...,alphap,beta1,...,beta4
  
  #Initialize feedback term:
  feedb1 <- feedb2 <- feedb3 <- feedb4 <- plogis(eta)
  
  value <- 0
  for (t in c((p+1):T)) {
    ft01 <- eta + pardep[p+1]*feedb1 + pardep[p+2]*feedb2 + pardep[p+3]*feedb3 + pardep[p+4]*feedb4
	for(k in 1:p) ft01 <- ft01 + pardep[k]*obincodes[datanum[t-k],]
	feedb4 <- feedb3
	feedb3 <- feedb2
	feedb2 <- feedb1
    feedb1 <- plogis(ft01)
    ft01 <- c(0, feedb1, 1)
    
    value <- value - log(ft01[datanum[t]+1]-ft01[datanum[t]])
  }
  value
}


#Constraints for optimization A %*% par >= b:
#use Amat_logit and bvec_logit with p+4







#Fit logitAR(1,4) model with scalar dependence parameters:

start <- c(-3.7, -3.3, -2.8, 3.5,1.2,0.2,0.7,1.2)

logitar14.s <- suppressWarnings(constrOptim(start, ll_logitarp4_scalar, NULL, ui=Amat_logit(d,5), ci=bvec_logit(d,5), datanum=datanum, d=d))
logitar14.s$par
#-3.7677816 -3.2981825 -2.8154967  3.4910148  1.1687164  0.1562705  0.7075559  1.2351682
logitar14.s$value
#1754.148

#AIC and BIC:
AIC.logitar14.s <- 2*Tlen/(Tlen-1)*logitar14.s$value+2*length(logitar14.s$par)
BIC.logitar14.s <- 2*Tlen/(Tlen-1)*logitar14.s$value+log(Tlen)*length(logitar14.s$par)
c(AIC.logitar14.s, BIC.logitar14.s) #3525.115 3576.022



# parameters:
etaML14 <- logitar14.s$par[1:d]
alpha1ML14 <- logitar14.s$par[d+1]
beta1ML14 <- logitar14.s$par[d+2]
beta2ML14 <- logitar14.s$par[d+3]
beta3ML14 <- logitar14.s$par[d+4]
beta4ML14 <- logitar14.s$par[d+5]


#PIT:
nobins <- 10
PIT <- array(0, c(2,nobins+1))

for(j in c(1:nobins)){
  u <- j/nobins
  pitval <- 0
  
  feedb1 <- feedb2 <- feedb3 <- feedb4 <- plogis(etaML14)
  for(t in c(2:Tlen)){
    ft01 <- etaML14 + alpha1ML14*obincodes[datanum[t-1],] + beta1ML14*feedb1 + beta2ML14*feedb2 + beta3ML14*feedb3 + beta4ML14*feedb4
	feedb4 <- feedb3
	feedb3 <- feedb2
	feedb2 <- feedb1
    feedb1 <- plogis(ft01)
    ft <- c(0, feedb1, 1)
    if(ft[datanum[t]]<u){
      if(ft[(datanum[t]+1)]<u){
        pitval <- pitval+1
      }else{
        pitval <- pitval+ (u-ft[datanum[t]])/(ft[(datanum[t]+1)]-ft[datanum[t]])
      }
    }
  }
  PIT[1,j+1] <- pitval/(Tlen-1)
  PIT[2,j+1] <- PIT[1,j+1]-PIT[1,j]
}

PIT.freq <- as.vector(rep(((1:nobins)-0.5)/nobins, PIT[2,2:(nobins+1)]*1000))
PIT.hist.tpar <- hist(PIT.freq, plot=FALSE, breaks=nobins)
PIT.hist.tpar$density <- PIT.hist.tpar$counts/1000

plot(PIT.hist.tpar, freq=FALSE, main="fitted IdD-logitAR(1,4) model", ylab="PIT histogram", xlab="u", col="gray", ylim=c(0,0.15), yaxt="n")
axis(side=2)
abline(h=0.1, lwd=2)



#Marginal calibration:
tabf <- array(0, c(Tlen-1, d))
  
feedb1 <- feedb2 <- feedb3 <- feedb4 <- plogis(etaML14)
for(t in c(2:Tlen)){
    ft01 <- etaML14 + alpha1ML14*obincodes[datanum[t-1],] + beta1ML14*feedb1 + beta2ML14*feedb2 + beta3ML14*feedb3 + beta4ML14*feedb4
	feedb4 <- feedb3
	feedb3 <- feedb2
	feedb2 <- feedb1
    feedb1 <- plogis(ft01)
	
	tabf[t-1, ] <- feedb1
}

#Marginal calibration:
hatf
#0.3272685 0.6130161 0.8322836
colMeans(tabf)
#0.3163657 0.6096673 0.8380936



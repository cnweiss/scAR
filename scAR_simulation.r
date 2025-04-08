#Simulating and fitting ordinal soft-clipping model

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
#ordbin(2,4) #0 0 1 1





##########################
#Log-likelihood functions:
##########################




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
  # feedb <- hatf
  
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




#Conditional log-likelihood function for soft-clipping AR(p) process (par=c(parmarg,pardep) with parmarg=c(eta0,...,eta_{d-1})).
#Model version "MID-scAR" with non-decreasing diagonal parameters:

ll_scarp_diag <- function(par,datanum,d, pc){ #states 1:(d+1)
  T <- length(datanum)
  p <- (length(par)-d)/d
  eta <- par[1:d]
  pardep <- par[(d+1):(d+p*d)] #alpha1.diag,...,alphap.diag
  
  value <- 0
  for (t in c((p+1):T)) {
    ft01 <- eta
    for(k in 1:p) ft01 <- ft01 + pardep[((k-1)*d+1):((k-1)*d+d)]*obincodes[datanum[t-k],]
    ft01 <- c(0, sc(ft01, pc), 1)
    
    value <- value - log(ft01[datanum[t]+1]-ft01[datanum[t]])
  }
  value
}

#Constraints for optimization A %*% par >= b:
#etas in [0,1], non-decreasing, also non-decreasing diagonals
Amat_diag <- function(d,p){
  value <- array(0,c(d+1+p*d, d+p*d))
  
  value[1,1] <- 1 #first CDF >=0
  for(i in 2:d){ #differences >=0
    value[i,i] <- 1
    value[i,i-1] <- (-1)
  }
  value[d+1,d] <- (-1) #last CDF <=1
  
  #Non-negative and non-decreasing diagonals:
  for(k in 1:p){
    value[d+1+(k-1)*d+1, d+(k-1)*d+1] <- 1
	
    for(i in 2:d){
      value[d+1+(k-1)*d+i, d+(k-1)*d+i] <- 1
      value[d+1+(k-1)*d+i, d+(k-1)*d+i-1] <- (-1)
    }
  }
  
  value
}

bvec_diag <- function(d,p){
  value <- rep(0, d+1+p*d)
  value[d+1] <- (-1)
  value
}





#Conditional log-likelihood function for soft-clipping AR(p) process (par=c(parmarg,pardep) with parmarg=c(eta0,...,eta_{d-1})).
#Model version "CoM-scAR" with identical rows constant:

ll_scarp_constant <- function(par,datanum,d, pc){ #states 1:(d+1)
  T <- length(datanum)
  p <- (length(par)-d)/d
  eta <- par[1:d]
  pardep <- par[(d+1):(d+p)] #alpha1.row,...,alphap.row
  value <- 0
  for (t in c((p+1):T)) {
    ft01 <- eta
    for(k in 1:p) ft01 <- ft01 + sum(pardep[k]*obincodes[datanum[t-k],])
    ft01 <- c(0, sc(ft01, pc), 1)
    
    value <- value - log(ft01[datanum[t]+1]-ft01[datanum[t]])
  }
  value
}




#Conditional log-likelihood function for soft-clipping AR(1,1) process (par=c(parmarg,pardep) with parmarg=c(eta0,...,eta_{d-1})).
#Model version "CoM-scAR" with identical rows constant:

ll_scar11_constant <- function(par,datanum,d, pc){ #states 1:(d+1)
  T <- length(datanum)
  eta <- par[1:d]
  pardep <- par[(d+1):(d+2)] #alpha1.row,...,alphap.row
  
  #Initialize feedback term:
  feedb <- sc(eta, pc)
  # feedb <- hatf
  
  value <- 0
  for (t in c(2:T)) {
    feedb <- sc(eta + sum(pardep[1]*obincodes[datanum[t-1],]) + sum(pardep[2]*feedb), pc)
    ft01 <- c(0, feedb, 1)
    
    value <- value - log(ft01[datanum[t]+1]-ft01[datanum[t]])
  }
  value
}

#Constraints for optimization A %*% par >= b:
#etas in [0,1], non-decreasing
Amat_constant <- function(d,p){
  value <- array(0,c(d+1, d+p))
  
  value[1,1] <- 1 #first CDF >=0
  for(i in 2:d){ #differences >=0
    value[i,i] <- 1
    value[i,i-1] <- (-1)
  }
  #value[d+1,d] <- (-1) #last CDF <=1
  
  value
}

bvec_constant <- function(d,p){
  value <- rep(0, d+1)
  value[d+1] <- (-1)
  value
}





####################################
#Simulating and fitting time series:
####################################





#Model specification example:

#IdD-scAR(1,1) model

#Parameter for soft-clipping function:
delta <- 0.01

#Model specification
eta <- c(0.07, 0.18, 0.3)
d <- length(eta)
Id <- diag(1, d)
alpha <- 0.4
beta <-0.2
A <- alpha*Id
B <- beta*Id

Tlen <- 500
prerun <- 250


#Binarization:
nbincodes <- diag(1,d+1) #nominal binarization
obincodes <- array(0, c(d+1,d))
for(s in 0:d) obincodes[s+1,] <- ordbin(s,d)



#Prerun:
pvec <- diff(c(0,eta,1)) #Initialization
Xnum <- sample(1:(d+1), 1, replace=TRUE, prob=pvec)
feedb <- sc(eta, delta)

for(t in 1:prerun){
	C <- obincodes[Xnum,] #vector C_{t-1}
	feedb <-sc(eta + A %*% C + B %*% feedb, delta)
	ft <- c(0, feedb, 1)
	pvec <- diff(ft) #Initialization

	Xnum <- sample(1:(d+1), 1, replace=TRUE, prob=pvec)	#Xt
}

#Main data
datanum <- rep(NA, Tlen)

for(t in 1:Tlen){
	C <- obincodes[Xnum,] #vector C_{t-1}
	feedb <-sc(eta + A %*% C + B %*% feedb, delta)
	ft <- c(0, feedb, 1)
	pvec <- diff(ft) #Initialization

	Xnum <- sample(1:(d+1), 1, replace=TRUE, prob=pvec)	#Xt
	datanum[t] <- Xnum
}

#Initialize scAR(1,1) estimates:
start <- c(eta, alpha, beta)


scar11.s <- suppressWarnings(constrOptim(start, ll_scar11_scalar, NULL, ui=Amat_scalar(d,2), ci=bvec_scalar(d,2), datanum=datanum, d=d, pc=delta))




#Model specification example:

#MID-scAR(1) model

#Parameter for soft-clipping function:
delta <- 0.01

#Model specification
eta <- c(0.07, 0.18, 0.3)
d <- length(eta)
Id <- diag(1, d)
a.diag <- c(0.35,0.4,0.45)
A <- diag(a.diag)


Tlen <- 100
prerun <- 250


#Binarization:
nbincodes <- diag(1,d+1) #nominal binarization
obincodes <- array(0, c(d+1,d))
for(s in 0:d) obincodes[s+1,] <- ordbin(s,d)


#Prerun:
pvec <- diff(c(0,eta,1)) #Initialization
Xnum <- sample(1:(d+1), 1, replace=TRUE, prob=pvec)

for(t in 1:prerun){
	C <- obincodes[Xnum,] #vector C_{t-1}
	ft <- c(0, sc(c(eta + A %*% C), delta), 1)
	pvec <- diff(ft) #Initialization

	Xnum <- sample(1:(d+1), 1, replace=TRUE, prob=pvec)	#Xt
}

#Main data
datanum <- rep(NA, Tlen)

for(t in 1:Tlen){
	C <- obincodes[Xnum,] #vector C_{t-1}
	ft <- c(0, sc(c(eta + A %*% C), delta), 1)
	pvec <- diff(ft) #Initialization

	Xnum <- sample(1:(d+1), 1, replace=TRUE, prob=pvec)	#Xt
	datanum[t] <- Xnum
}


#Initialize scAR(1) estimates:
start <- c(eta, a.diag)

scar1.d <- suppressWarnings(constrOptim(start, ll_scarp_diag, NULL, ui=Amat_diag(d,1), ci=bvec_diag(d,1), datanum=datanum, d=d, pc=delta))




#Model specification example:

#CoM-scAR(1,1) model

#Parameter for soft-clipping function:
delta <- 0.01

#Model specification
eta <- c(0.25,0.4,0.6)
d <- length(eta)
Id <- diag(1, d)
E <- array(1, c(d,d))
alpha <- -0.1
beta <- -0.1
A <- alpha*E
B <- beta*E

Tlen <- 500
prerun <- 250


#Binarization:
nbincodes <- diag(1,d+1) #nominal binarization
obincodes <- array(0, c(d+1,d))
for(s in 0:d) obincodes[s+1,] <- ordbin(s,d)



#Prerun:
pvec<- diff(c(0,solve(Id-A-B) %*% eta,1)) #Initialization
Xnum <- sample(1:(d+1), 1, replace=TRUE, prob=pvec)

feedb <- sc(eta, delta)
for(t in 1:prerun){
	C <- obincodes[Xnum,] #vector C_{t-1}
	feedb <-sc(eta + A %*% C + B %*% feedb, delta)
	ft <- c(0, feedb, 1)
	pvec <- diff(ft) #Initialization

	Xnum <- sample(1:(d+1), 1, replace=TRUE, prob=pvec)	#Xt
}

#Main data
datanum <- rep(NA, Tlen)

for(t in 1:Tlen){
	C <- obincodes[Xnum,] #vector C_{t-1}
	feedb <-sc(eta + A %*% C + B %*% feedb, delta)
	ft <- c(0, feedb, 1)
	pvec <- diff(ft) #Initialization

	Xnum <- sample(1:(d+1), 1, replace=TRUE, prob=pvec)	#Xt
	datanum[t] <- Xnum
}

#Initialize scARMA(1,1) estimates:
start <- c(eta, alpha , beta)


scar11.c <- suppressWarnings(constrOptim(start, ll_scar11_constant, NULL, ui=Amat_constant(d,2), ci=bvec_constant(d,2), datanum=datanum, d=d, pc=delta))


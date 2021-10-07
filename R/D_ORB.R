# library(statip) # Needed for obtaining kernel density estimates (KDEs)
# library(dad)    # Needed for calculating Hellinger distance between KDEs

####################################
## D_ORB MEASURE OF DISSIMILARITY ##
####################################
# Obtains a measure of the dissimilarity in posterior distributions
# due to ORB

# INPUTS:
# samples.selectionmodel = MCMC samples under ABSORB model
# samples.nobiasmodel = MCMC samples under no-bias model (rho1=rho2=0)

# OUTPUT:
# D measure = measure of dissimilarity in estimates from the two models

D_ORB <- function(samples.ABSORBmodel, samples.nobiasmodel){

  if(dim(as.matrix(samples.ABSORBmodel))[2] != dim(as.matrix(samples.nobiasmodel))[2])
    stop("Error: Please make sure that the dimensions of the two MCMC samples are equal.")

  if(dim(as.matrix(absorb.mod$mu1.samples))[2]==1){

    fx = statip::densityfun(samples.ABSORBmodel)
    fy = statip::densityfun(samples.nobiasmodel)

    g <- function(z) (fx(z)^0.5 - fy(z)^0.5)^2

    # Compute Hellinger distance
    H.sq = stats::integrate(g, lower=-Inf, upper=Inf, subdivisions=500, stop.on.error=FALSE)$value/2
    H.hat = sqrt(H.sq)

    if(H.hat > 1){ # sometimes returns values slightly larger than 1
    H.hat = 1
    }
  } else {
    samples.ABSORBmodel = as.data.frame(samples.ABSORBmodel)
    samples.nobiasmodel = as.data.frame(samples.nobiasmodel)

    H.hat = dad::hellinger(samples.ABSORBmodel, samples.nobiasmodel)
  }

    # Return the Hellinger distance
    return(H.hat)

}

# library(rjags) # implements MCMC using JAGS

##################
## ABSORB model ##
##################
# Implements the ABSORB model

# INPUTS:
# study_size = study sizes for the individual studies: n_1, ..., n_n. Needed for
#              estimating the missing studies' standard errors
# y1 = given vector of observed treatment effects for the first endpoint.
#      If treatment is missing for the ith study, should put an NA in the ith entry of the vector.
# s1 = given vector of within-study standard errors for first endpoint.
#      If treatment is missing for the ith study, should put an NA in the ith entry of the vector.
# y2 = given vector of observed treatment effects for the second endpoint.
#      If treatment is missing for the ith study, should put an NA in the ith entry of the vector.
# s2 = given vector of within-study standard errors for the second endpoint.
#      If treatment is missing for the ith study, should put an NA for that entry of the vector.
# seed = seed. Set this if you want to reproduce exact results
# burn = number of burn-in samples. Default is 10,000
# nmc = number of posterior draws to be saved in each of 3 MCMC chains.
#       Default is 50,000 which returns 150,000 posterior samples total

# OUTPUTS:
# mu1.hat = point estimate for mu1. Returns posterior mean of mu1
# mu2.hat = point estimate for mu2. Returns posterior mean of mu2
# tau1.hat = point estimate for tau1. Returns posterior mean of tau1
# tau2.hat = point estimate for tau1. Returns posterior mean of tau2
# gamma01.hat = point estimate for gamma01. Returns posterior median for gamma01
# gamma02.hat = point estimate for gamma02. Returns posterior median for gamma02
# gamma11.hat = point estimate for gamma11. Returns posterior median for gamma11
# gamma12.hat = point estimate for gamma12. Returns posterior median for gamma12
# rho1.hat = point estimate for rho1. Returns posterior median for rho1
# rho2.hat = point estimate for rho2. Returns posterior median for rho2
# rho3.hat = point estimate for rho3. Returns posterior median for rho3
# rho4.hat = point estimate for rho4. Returns posterior median for rho4

# mu1.samples = mu samples. Used for plotting the posterior for mu1
#               and performing inference for mu1
# mu2.samples = mu2 samples. Used for plotting the posterior for mu2
#               and performing inference for mu2
# tau1.samples = tau1 samples. Used for plotting the posterior for tau1
#                and performing inference for tau1
# tau2.samples = tau2 samples. Used for plotting the posterior for tau2
#                and performing inference for tau2
# gamma01.samples = gamma01 samples. Used for plotting the posterior for gamma01
#                   and performing inference for gamma01
# gamma02.samples = gamma02 samples. Used for plotting the posterior for gamma02
#                   and performing inference for gamma02
# gamma11.samples = gamma11 samples. Used for plotting the posterior for gamma11
#                   and performing inference for gamma11
# gamma12.samples = gamma12 samples. Used for plotting the posterior for gamma12
#                   and performing inference for gamma12
# rho1.samples = rho1 samples. Used for plotting the posterior for rho1
#                and performing inference for rho1
# rho2.samples = rho2 samples. Used for plotting the posterior for rho2
#                and performing inference for rho2
# rho3.samples = rho3 samples. Used for plotting the posterior for rho3
#                and performing inference for rho3
# rho4.samples = rho4 samples. Used for plotting the posterior for rho4
#                and performing inference for rho4

# mu1.chains = matrix of MCMC chains for 3 separate chains of mu1.
#              # of rows = nmc, # of columns = 3
# mu2.chains = matrix of MCMC chains for the separate chains of mu2.
#              # of rows = nmc, # of columns = 3

ABSORB <- function(study_sizes, y1, s1, y2, s2, seed=NULL,
                   burn=10000, nmc=50000){

  ##################
  # INITIAL CHECKS #
  ##################

  if(any(is.nan(y1)) || any(is.nan(y2)))
    stop("ERROR: Effect sizes cannot be NaN.")
  if(any(is.infinite(y1)) || any(is.infinite(y2)))
    stop("ERROR: Effect sizes cannot be Inf or -Inf.")
  if(length(study_sizes) != length(y1) || length(study_sizes) != length(y2))
    stop("ERROR: The vector of study sizes is not the same length as y1 or y2.")
  if(length(which(is.na(study_sizes)) != 0))
    stop("ERROR: All study sizes need to be provided.")
  if(length(y1) != length(y2))
    stop("ERROR: Please ensure that 'y1' and 'y2' are of the same length. If there are missing values in y1 or y2, please put an NA for that entry.")
  # Check that y1 and s1 are of the same length
  if(length(y1) != length(s1))
    stop("ERROR: Please ensure that 'y1' and 's1' are of the same length.")
  if(length(y2) != length(s2))
    stop("ERROR: Please ensure that 'y2' and 's2' are of the same length.")
  if(!identical(which(is.na(y1)),which(is.na(s1))))
    stop("ERROR: There is a mismatch between missing outcomes in y1 and missing standard errors in s1. If an entry in y1 is NA, the corresponding entry in s1 should be NA.")
  if(!identical(which(is.na(y2)),which(is.na(s2))))
    stop("ERROR: There is a mismatch between missing outcomes in y2 and missing standard errors in s2. If an entry in y2 is NA, the correpsonding entry in s2 should be NA.")
  if(any(s1<=0, na.rm=TRUE))
    stop("ERROR: Please ensure that all non-missing standard errors are strictly positive.")
  if(any(s2<=0, na.rm=TRUE))
    stop("ERROR: Please ensure that all non-missing standard errors are strictly positive.")
  if(any(study_sizes<=0))
    stop("ERROR: Please ensure that all study sizes are strictly positive.")
  if(length(unique(y1[!is.na(y1)]))==1)
    stop("ERROR: All non-missing y1's are identical.")
  if(length(unique(y2[!is.na(y2)]))==1)
    stop("ERROR: All non-missing y2's are identical.")

  ########################
  # Pre-process the data #
  ########################

  # Construct the 2xm1 matrix of responses that report both outcomes
  n = length(y1)
  m1 = length(which(!is.na(y1) & !is.na(y2)))
  y.both = matrix(0, nrow=m1, ncol=2)
  s.both = matrix(0, nrow=m1, ncol=2)
  ind.both = 1

  if(n==m1){
    y.both = cbind(y1,y2)
    s.both = cbind(s1,s2)
  } else {
    for(i in 1:n){
      if(!is.na(y1[i]) & !is.na(y2[i])){
        y.both[ind.both,1] = y1[i]
        y.both[ind.both,2] = y2[i]
        s.both[ind.both,1] = s1[i]
        s.both[ind.both,2] = s2[i]

        if(ind.both < m1){
          ind.both = ind.both + 1
        }
      }
    }
  }

  # Check if both columns of y.both are equal
  if(all(y.both[,1]==y.both[,2])){
    stop("After removing missing values, both endpoints are exactly identical.
         Algorithm will not converge if non-missing y1 and y2 are exactly equal.
         Please verify that y1 and y2 are not the the same outcome.")
  }
  # Check if both columns of y.both are nearly identical
  num.identical = length(which(y.both[,1]==y.both[,2]))
  if(num.identical >= m1-4 || cor(y.both[,1],y.both[,2]) >= 0.99){
    stop("After removing missing values, both endpoints are nearly identical.
         Algorithm might not converge. Please verify that y1 and y2 are not the
         same outcome.")
  }

  # Construct the m2x1 vector of studies that only report y1
  m2 = length(which(!is.na(y1) & is.na(y2)))
  if(m2!=0){
    y1.only = rep(0, m2)
    s1.only = rep(0, m2)
    studysizes.y1.only = rep(0, m2)
    ind.firstonly = 1

    for(i in 1:n){
      if(!is.na(y1[i]) & is.na(y2[i])){
        y1.only[ind.firstonly] = y1[i]
        s1.only[ind.firstonly] = s1[i]
        studysizes.y1.only[ind.firstonly] = study_sizes[i]

        if(ind.firstonly < m2){
          ind.firstonly = ind.firstonly + 1
        }
      }
    }

    # Construct the m2x1 vector of the imputed standard errors for the missing s_2's
    s2.imputed = rep(0, m2)
    s2.nonmissing.indices = setdiff(seq(1:n), which(!is.na(y1) & is.na(y2)))
    s2.nonmissing = s2[s2.nonmissing.indices]
    s2.nonmissing.studysizes = study_sizes[s2.nonmissing.indices]
    # Estimate k2.hat
    k2.hat = sum(1/s2.nonmissing^2)/sum(s2.nonmissing.studysizes)

    # Impute the missing standard errors
    ind.firstonly = 1
    for(i in 1:n){
      if(!is.na(y1[i]) & is.na(y2[i])){
        s2.imputed[ind.firstonly] = sqrt(1/((k2.hat)*studysizes.y1.only[ind.firstonly]))

        if(ind.firstonly < m2){
          ind.firstonly = ind.firstonly + 1
        }
      }
    }
  }

  # Construct the m3x1 vector of studies that only report y2
  m3 = length(which(is.na(y1) & !is.na(y2)))
  if(m3!=0){
    y2.only = rep(0, m3)
    s2.only = rep(0, m3)
    studysizes.y2.only = rep(0, m3)
    ind.secondonly = 1

    for(i in 1:n){
      if(is.na(y1[i]) & !is.na(y2[i])){
        y2.only[ind.secondonly] = y2[i]
        s2.only[ind.secondonly] = s2[i]
        studysizes.y2.only[ind.secondonly] = study_sizes[i]

        if(ind.secondonly < m3){
          ind.secondonly = ind.secondonly + 1
        }
      }
    }

    # Construct the m3x1 vector of the imputed standard errors for the missing s_1's
    s1.imputed = rep(0, m3)
    s1.nonmissing.indices = setdiff(seq(1:n), which(is.na(y1) & !is.na(y2)))
    s1.nonmissing = s1[s1.nonmissing.indices]
    s1.nonmissing.studysizes = study_sizes[s1.nonmissing.indices]
    # Estimate k2.hat
    k1.hat = sum(1/s1.nonmissing^2)/sum(s1.nonmissing.studysizes)

    # Impute the missing standard errors
    ind.secondonly = 1
    for(i in 1:n){
      if(is.na(y1[i]) & !is.na(y2[i])){
        s1.imputed[ind.secondonly] = sqrt(1/((k1.hat)*studysizes.y2.only[ind.secondonly]))

        if(ind.secondonly < m3){
          ind.secondonly = ind.secondonly + 1
        }
      }
    }
  }


  #########################
  # Initialization values #
  # for the 3 chains      #
  #########################

  data1 = data.frame(cbind(y1,s1))
  max.s1 = max(s1, na.rm=TRUE)
  data2 = data.frame(cbind(y2,s2))
  max.s2 = max(s2, na.rm=TRUE)


  # obtain initial estimates for mu1
  mu1.init1=mean(data1[,1], na.rm=TRUE) # estimate of mu1 based on all data
  if(m2>=1){
    mu1.init2=mean(y1.only) # estimate of mu1 based only on studies that reported only y1
  } else {
    mu1.init2=mu1.init1
  }
  mu1.init3=mean(y.both[,1]) # estimate of mu1 based only on studies that reported both outcomes

  # moment estimators for tau1
  tau1.init1=abs(sd(data1[,1], na.rm=TRUE)^2-(sd(data1[,2], na.rm=TRUE)^2+mean(data1[,2], na.rm=TRUE)^2))
  if(m2>1){
    tau1.init2=abs(sd(y1.only)^2-(sd(s1.only)^2+mean(s1.only)^2))
  } else if(m2<=1){
    tau1.init2=tau1.init1
  }
  tau1.init3=abs(sd(y.both[,1])^2-(sd(s.both[,1])^2+mean(s.both[,1])^2))

  # obtain initial estimates for mu2
  mu2.init1=mean(data2[,1], na.rm=TRUE)
  if(m3>=1){
    mu2.init2=mean(y2.only)
  } else {
    mu2.init2=mu2.init1
  }
  mu2.init3=mean(y.both[,2])

  # moment estimators for tau2
  tau2.init1=abs(sd(data2[,1], na.rm=TRUE)^2-(sd(data2[,2], na.rm=TRUE)^2+mean(data2[,2], na.rm=TRUE)^2))
  if(m3>1){
    tau2.init2=abs(sd(y2.only)^2-(sd(s2.only)^2+mean(s2.only)^2))
  } else if(m3<=1){
    tau2.init2=tau2.init1
  }
  tau2.init3=abs(sd(y.both[,2])^2-(sd(s.both[,2])^2+mean(s.both[,2])^2))

  # initial values for other parameters
  gamma01.init=-1
  gamma02.init=-1
  gamma11.init=0
  gamma12.init=0
  rho1.init=0
  rho2.init=0
  rho3.init=0
  rho4.init=0


  # List for initial values

  if(is.null(seed)){

    # if no seed is specified
    # initialization for chain 1
    init.vals1 = list(mu1 = mu1.init1,
                     mu2 = mu2.init1,
                     tau1 = tau1.init1,
                     tau2 = tau2.init1,
                     gamma01 = gamma01.init,
                     gamma02 = gamma02.init,
                     gamma11 = gamma11.init,
                     gamma12 = gamma12.init,
                     rho1 = rho1.init,
                     rho2 = rho2.init,
                     rho3 = rho3.init,
                     rho4 = rho4.init)
    # initialization for chain 2
    init.vals2 = list(mu1 = mu1.init2,
                      mu2 = mu2.init2,
                      tau1 = tau1.init2,
                      tau2 = tau2.init2,
                      gamma01 = gamma01.init,
                      gamma02 = gamma02.init,
                      gamma11 = gamma11.init,
                      gamma12 = gamma12.init,
                      rho1 = rho1.init,
                      rho2 = rho2.init,
                      rho3 = rho3.init,
                      rho4 = rho4.init)
    # initialization for chain 3
    init.vals3 = list(mu1 = mu1.init3,
                      mu2 = mu2.init3,
                      tau1 = tau1.init3,
                      tau2 = tau2.init3,
                      gamma01 = gamma01.init,
                      gamma02 = gamma02.init,
                      gamma11 = gamma11.init,
                      gamma12 = gamma12.init,
                      rho1 = rho1.init,
                      rho2 = rho2.init,
                      rho3 = rho3.init,
                      rho4 = rho4.init)

  } else if(!is.null(seed)){

    # if a seed is specified
    seed <- as.integer(seed)
    # initialization for chain 1
    init.vals1 = list(.RNG.name = "base::Marsaglia-Multicarry",
                      .RNG.seed = seed,
                      mu1 = mu1.init1,
                      mu2 = mu2.init1,
                      tau1 = tau1.init1,
                      tau2 = tau2.init1,
                      gamma01 = gamma01.init,
                      gamma02 = gamma02.init,
                      gamma11 = gamma11.init,
                      gamma12 = gamma12.init,
                      rho1 = rho1.init,
                      rho2 = rho2.init,
                      rho3 = rho3.init,
                      rho4 = rho4.init)
    # initialization for chain 2
    init.vals2 = list(.RNG.name = "base::Super-Duper",
                      .RNG.seed = seed,
                      mu1 = mu1.init2,
                      mu2 = mu2.init2,
                      tau1 = tau1.init2,
                      tau2 = tau2.init2,
                      gamma01 = gamma01.init,
                      gamma02 = gamma02.init,
                      gamma11 = gamma11.init,
                      gamma12 = gamma12.init,
                      rho1 = rho1.init,
                      rho2 = rho2.init,
                      rho3 = rho3.init,
                      rho4 = rho4.init)
    # initialization for chain 3
    init.vals3 = list(.RNG.name = "base::Mersenne-Twister",
                      .RNG.seed = seed,
                      mu1 = mu1.init3,
                      mu2 = mu2.init3,
                      tau1 = tau1.init3,
                      tau2 = tau2.init3,
                      gamma01 = gamma01.init,
                      gamma02 = gamma02.init,
                      gamma11 = gamma11.init,
                      gamma12 = gamma12.init,
                      rho1 = rho1.init,
                      rho2 = rho2.init,
                      rho3 = rho3.init,
                      rho4 = rho4.init)
  }


  ##################
  # Fit JAGS model #
  ##################

  if(m2!=0 & m3!=0){

    # Fit model if both endpoints contain missing values

    model_string <- "model{

      # Likelihood

      # For the studies that report BOTH endpoints
      for(i in 1:m1){
        # Sample 2x1 vector z
        z.both[i,1] ~ dnorm(gamma01+gamma11/s.both[i,1],1) T(0,)
        z.both[i,2] ~ dnorm(gamma02+gamma12/s.both[i,2],1) T(0,)

        # Conditional mean for y given z
        y.both.mean[i,1] <- theta.both[i,1] + rho1*s.both[i,1]*(z.both[i,1]-gamma01-gamma11/s.both[i,1])
        y.both.mean[i,2] <- theta.both[i,2] + rho2*s.both[i,2]*(z.both[i,2]-gamma02-gamma12/s.both[i,2])

        # Conditional covariance matrix for y given z
        Sigma.both[i,1,1] <- pow(s.both[i,1],2)*(1-pow(rho1,2))
        Sigma.both[i,1,2] <- rho3*s.both[i,1]*s.both[i,2]
        Sigma.both[i,2,1] <- rho3*s.both[i,1]*s.both[i,2]
        Sigma.both[i,2,2] <- pow(s.both[i,2],2)*(1-pow(rho2,2))

        # Sample y conditionally on z
        y.both[i,1:2] ~ dmnorm.vcov(y.both.mean[i,1:2], Sigma.both[i,1:2,1:2])
      }

      # For the studies that report first endpoint ONLY
      for(i in 1:m2){
        # Sample z1 for the studies that reported first endpoint ONLY
        z1.only[i] ~ dnorm(gamma01+gamma11/s1.only[i],1) T(0,)

        # Conditional mean for y1 given z1
        y1.only.mean[i] <- theta1.only[i] + rho1*s1.only[i]*(z1.only[i]-gamma01-gamma11/s1.only[i])

        # Conditional variance for y1 given z1
        y1.only.var[i] <- pow(s1.only[i],2)*(1-pow(rho1,2))

        # Sample y1 conditionally on z1
        y1.only[i] ~ dnorm(y1.only.mean[i], 1/y1.only.var[i])

        # Sample z2 for the studies that did NOT report first endpoint
        z2.missing[i] ~ dnorm(gamma02+gamma12/s2.imputed[i],1) T(,0)
      }

      # For the studies that report second endpoint ONLY
      for(i in 1:m3){
        # Sample z2 for the studies that reported second endpoint ONLY
        z2.only[i] ~ dnorm(gamma02+gamma12/s2.only[i],1) T(0,)

        # Conditional mean for y2 given z2
        y2.only.mean[i] <- theta2.only[i] + rho2*s2.only[i]*(z2.only[i]-gamma02-gamma12/s2.only[i])

        # Conditional variance for y2 given z2
        y2.only.var[i] <- pow(s2.only[i],2)*(1-pow(rho2,2))

        # Sample y2 conditionally on z2
        y2.only[i] ~ dnorm(y2.only.mean[i], 1/y2.only.var[i])

        # Sample z1 for the studies that did NOT report first endpoint
        z1.missing[i] ~ dnorm(gamma01+gamma11/s1.imputed[i],1) T(,0)
      }

      # Prior for the theta's where both endpoints are reported
      theta.mean[1] <- mu1
      theta.mean[2] <- mu2
      thetacov.both[1,1] <- pow(tau1,2)
      thetacov.both[1,2] <- rho4*tau1*tau2
      thetacov.both[2,1] <- rho4*tau1*tau2
      thetacov.both[2,2] <- pow(tau2,2)

      # Prior for theta's for studies that reported BOTH endpoints
      for(i in 1:m1){
        theta.both[i,1:2] ~ dmnorm.vcov(theta.mean, thetacov.both)
      }

      # Prior for theta1 for studies that reported first endpoint ONLY
      for(i in 1:m2){
        theta1.only[i] ~ dnorm(mu1, 1/pow(tau1,2))
      }

      # Prior for theta2 for studies that reported second endpoint ONLY
      for(i in 1:m3){
        theta2.only[i] ~ dnorm(mu2, 1/pow(tau2,2))
      }

      # Prior for mu1
      mu1 ~ dnorm(0,0.0001)

      # Prior for mu2
      mu2 ~ dnorm(0,0.0001)

      # Prior for tau1. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
      tau1 ~ dt(0,1,1) T(0,)

      # Prior for tau2. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
      tau2 ~ dt(0,1,1) T(0,)

      # Prior for gamma01
      gamma01 ~ dunif(-2,2)

      # Prior for gamma02
      gamma02 ~ dunif(-2,2)

      # Prior for gamma11
      gamma11 ~ dunif(0, max.s1)

      # Prior for gamma12
      gamma12 ~ dunif(0, max.s2)

      # Prior for rho1
      rho1 ~ dunif(-1,1)

      # Prior for rho2
      rho2 ~ dunif(-1,1)

      # Prior for rho3
      rho3 ~ dunif(-1,1)

      # Prior for rho4
      rho4 ~ dunif(-1,1)

    }"

    # Compile the model in JAGS
    model <- rjags::jags.model(textConnection(model_string),
                               data = list(y.both=y.both,s.both=s.both,
                                           y1.only=y1.only,s1.only=s1.only,
                                           y2.only=y2.only,s2.only=s2.only,
                                           s1.imputed=s1.imputed,
                                           s2.imputed=s2.imputed,
                                           m1=m1,m2=m2,m3=m3,
                                           max.s1=max.s1,max.s2=max.s2),
                                           inits = list(init.vals1,init.vals2,init.vals3),
                                           n.chains=3,
                                           quiet=TRUE)

    # Draw samples
    # First use function 'update' to draw warm-up (burnin) samples
    stats::update(model, n.iter=burn, progress.bar="none") # Burnin for 10,000 samples

    # Next use the function 'coda.samples' to produce the next nmc
    # samples which will ultimately used to approximate the posterior.
    full.samples <- rjags::coda.samples(model,
                                        variable.names=c("mu1","mu2","tau1","tau2",
                                                         "gamma01","gamma02","gamma11","gamma12",
                                                         "rho1","rho2","rho3","rho4"),
                                        n.iter=nmc,
                                        progress.bar="none")

  } else if(m2==0 & m3==0) {

    # Fit model when NO endpoints are missing

    model_string <- "model{

      # Likelihood

      # For the studies that report BOTH endpoints
      for(i in 1:m1){
        # Sample 2x1 vector z
        z.both[i,1] ~ dnorm(gamma01+gamma11/s.both[i,1],1) T(0,)
        z.both[i,2] ~ dnorm(gamma02+gamma12/s.both[i,2],1) T(0,)

        # Conditional mean for y given z
        y.both.mean[i,1] <- theta.both[i,1] + rho1*s.both[i,1]*(z.both[i,1]-gamma01-gamma11/s.both[i,1])
        y.both.mean[i,2] <- theta.both[i,2] + rho2*s.both[i,2]*(z.both[i,2]-gamma02-gamma12/s.both[i,2])

        # Conditional covariance matrix for y given z
        Sigma.both[i,1,1] <- pow(s.both[i,1],2)*(1-pow(rho1,2))
        Sigma.both[i,1,2] <- rho3*s.both[i,1]*s.both[i,2]
        Sigma.both[i,2,1] <- rho3*s.both[i,1]*s.both[i,2]
        Sigma.both[i,2,2] <- pow(s.both[i,2],2)*(1-pow(rho2,2))

        # Sample y conditionally on z
        y.both[i,1:2] ~ dmnorm.vcov(y.both.mean[i,1:2], Sigma.both[i,1:2,1:2])
      }

      # Prior for the theta's where both endpoints are reported
      theta.mean[1] <- mu1
      theta.mean[2] <- mu2
      thetacov.both[1,1] <- pow(tau1,2)
      thetacov.both[1,2] <- rho4*tau1*tau2
      thetacov.both[2,1] <- rho4*tau1*tau2
      thetacov.both[2,2] <- pow(tau2,2)

      # Prior for theta's for studies that reported BOTH endpoints
      for(i in 1:m1){
        theta.both[i,1:2] ~ dmnorm.vcov(theta.mean, thetacov.both)
      }

      # Prior for mu1
      mu1 ~ dnorm(0,0.0001)

      # Prior for mu2
      mu2 ~ dnorm(0,0.0001)

      # Prior for tau1. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
      tau1 ~ dt(0,1,1) T(0,)

      # Prior for tau2. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
      tau2 ~ dt(0,1,1) T(0,)

      # Prior for gamma01
      gamma01 ~ dunif(-2,2)

      # Prior for gamma02
      gamma02 ~ dunif(-2,2)

      # Prior for gamma11
      gamma11 ~ dunif(0, max.s1)

      # Prior for gamma12
      gamma12 ~ dunif(0, max.s2)

      # Prior for rho1
      rho1 ~ dunif(-1,1)

      # Prior for rho2
      rho2 ~ dunif(-1,1)

      # Prior for rho3
      rho3 ~ dunif(-1,1)

      # Prior for rho4
      rho4 ~ dunif(-1,1)

    }"

    # Compile the model in JAGS
    model <- rjags::jags.model(textConnection(model_string),
                               data = list(y.both=y.both,s.both=s.both,
                                           m1=m1,
                                           max.s1=max.s1,max.s2=max.s2),
                               inits = list(init.vals1,init.vals2,init.vals3),
                               n.chains=3,
                               quiet=TRUE)

    # Draw samples
    # First use function 'update' to draw warm-up (burnin) samples
    stats::update(model, n.iter=burn, progress.bar="none") # Burnin for 10,000 samples

    # Next use the function 'coda.samples' to produce the next nmc
    # samples which will ultimately used to approximate the posterior.
    full.samples <- rjags::coda.samples(model,
                                        variable.names=c("mu1","mu2","tau1","tau2",
                                                         "gamma01","gamma02","gamma11","gamma12",
                                                         "rho1","rho2","rho3","rho4"),
                                        n.iter=nmc,
                                        progress.bar="none")

  } else if(m2==0 & m3!=0){

    # Fit model where only y1's are missing but no y2's are missing

    model_string <- "model{

      # Likelihood

      # For the studies that report BOTH endpoints
      for(i in 1:m1){
        # Sample 2x1 vector z
        z.both[i,1] ~ dnorm(gamma01+gamma11/s.both[i,1],1) T(0,)
        z.both[i,2] ~ dnorm(gamma02+gamma12/s.both[i,2],1) T(0,)

        # Conditional mean for y given z
        y.both.mean[i,1] <- theta.both[i,1] + rho1*s.both[i,1]*(z.both[i,1]-gamma01-gamma11/s.both[i,1])
        y.both.mean[i,2] <- theta.both[i,2] + rho2*s.both[i,2]*(z.both[i,2]-gamma02-gamma12/s.both[i,2])

        # Conditional covariance matrix for y given z
        Sigma.both[i,1,1] <- pow(s.both[i,1],2)*(1-pow(rho1,2))
        Sigma.both[i,1,2] <- rho3*s.both[i,1]*s.both[i,2]
        Sigma.both[i,2,1] <- rho3*s.both[i,1]*s.both[i,2]
        Sigma.both[i,2,2] <- pow(s.both[i,2],2)*(1-pow(rho2,2))

        # Sample y conditionally on z
        y.both[i,1:2] ~ dmnorm.vcov(y.both.mean[i,1:2], Sigma.both[i,1:2,1:2])
      }

      # For the studies that report second endpoint ONLY
      for(i in 1:m3){
        # Sample z2 for the studies that reported second endpoint ONLY
        z2.only[i] ~ dnorm(gamma02+gamma12/s2.only[i],1) T(0,)

        # Conditional mean for y2 given z2
        y2.only.mean[i] <- theta2.only[i] + rho2*s2.only[i]*(z2.only[i]-gamma02-gamma12/s2.only[i])

        # Conditional variance for y2 given z2
        y2.only.var[i] <- pow(s2.only[i],2)*(1-pow(rho2,2))

        # Sample y2 conditionally on z2
        y2.only[i] ~ dnorm(y2.only.mean[i], 1/y2.only.var[i])

        # Sample z1 for the studies that did NOT report first endpoint
        z1.missing[i] ~ dnorm(gamma01+gamma11/s1.imputed[i],1) T(,0)
      }

      # Prior for the theta's where both endpoints are reported
      theta.mean[1] <- mu1
      theta.mean[2] <- mu2
      thetacov.both[1,1] <- pow(tau1,2)
      thetacov.both[1,2] <- rho4*tau1*tau2
      thetacov.both[2,1] <- rho4*tau1*tau2
      thetacov.both[2,2] <- pow(tau2,2)

      # Prior for theta's for studies that reported BOTH endpoints
      for(i in 1:m1){
        theta.both[i,1:2] ~ dmnorm.vcov(theta.mean, thetacov.both)
      }

      # Prior for theta2 for studies that reported second endpoint ONLY
      for(i in 1:m3){
        theta2.only[i] ~ dnorm(mu2, 1/pow(tau2,2))
      }

      # Prior for mu1
      mu1 ~ dnorm(0,0.0001)

      # Prior for mu2
      mu2 ~ dnorm(0,0.0001)

      # Prior for tau1. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
      tau1 ~ dt(0,1,1) T(0,)

      # Prior for tau2. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
      tau2 ~ dt(0,1,1) T(0,)

      # Prior for gamma01
      gamma01 ~ dunif(-2,2)

      # Prior for gamma02
      gamma02 ~ dunif(-2,2)

      # Prior for gamma11
      gamma11 ~ dunif(0, max.s1)

      # Prior for gamma12
      gamma12 ~ dunif(0, max.s2)

      # Prior for rho1
      rho1 ~ dunif(-1,1)

      # Prior for rho2
      rho2 ~ dunif(-1,1)

      # Prior for rho3
      rho3 ~ dunif(-1,1)

      # Prior for rho4
      rho4 ~ dunif(-1,1)

    }"

    # Compile the model in JAGS
    model <- rjags::jags.model(textConnection(model_string),
                               data = list(y.both=y.both,s.both=s.both,
                                           y2.only=y2.only,s2.only=s2.only,
                                           s1.imputed=s1.imputed,
                                           m1=m1,m3=m3,
                                           max.s1=max.s1,max.s2=max.s2),
                               inits = list(init.vals1,init.vals2,init.vals3),
                               n.chains=3,
                               quiet=TRUE)

    # Draw samples
    # First use function 'update' to draw warm-up (burnin) samples
    stats::update(model, n.iter=burn, progress.bar="none") # Burnin for 10,000 samples

    # Next use the function 'coda.samples' to produce the next nmc
    # samples which will ultimately used to approximate the posterior.
    full.samples <- rjags::coda.samples(model,
                                        variable.names=c("mu1","mu2","tau1","tau2",
                                                         "gamma01","gamma02","gamma11","gamma12",
                                                         "rho1","rho2","rho3","rho4"),
                                        n.iter=nmc,
                                        progress.bar="none")

  } else if(m2!=0 & m3==0){

    # Fit model where only y2's are missing but no y1's are missing

    model_string <- "model{

      # Likelihood

      # For the studies that report BOTH endpoints
      for(i in 1:m1){
        # Sample 2x1 vector z
        z.both[i,1] ~ dnorm(gamma01+gamma11/s.both[i,1],1) T(0,)
        z.both[i,2] ~ dnorm(gamma02+gamma12/s.both[i,2],1) T(0,)

        # Conditional mean for y given z
        y.both.mean[i,1] <- theta.both[i,1] + rho1*s.both[i,1]*(z.both[i,1]-gamma01-gamma11/s.both[i,1])
        y.both.mean[i,2] <- theta.both[i,2] + rho2*s.both[i,2]*(z.both[i,2]-gamma02-gamma12/s.both[i,2])

        # Conditional covariance matrix for y given z
        Sigma.both[i,1,1] <- pow(s.both[i,1],2)*(1-pow(rho1,2))
        Sigma.both[i,1,2] <- rho3*s.both[i,1]*s.both[i,2]
        Sigma.both[i,2,1] <- rho3*s.both[i,1]*s.both[i,2]
        Sigma.both[i,2,2] <- pow(s.both[i,2],2)*(1-pow(rho2,2))

        # Sample y conditionally on z
        y.both[i,1:2] ~ dmnorm.vcov(y.both.mean[i,1:2], Sigma.both[i,1:2,1:2])
      }

      # For the studies that report first endpoint ONLY
      for(i in 1:m2){
        # Sample z1 for the studies that reported first endpoint ONLY
        z1.only[i] ~ dnorm(gamma01+gamma11/s1.only[i],1) T(0,)

        # Conditional mean for y1 given z1
        y1.only.mean[i] <- theta1.only[i] + rho1*s1.only[i]*(z1.only[i]-gamma01-gamma11/s1.only[i])

        # Conditional variance for y1 given z1
        y1.only.var[i] <- pow(s1.only[i],2)*(1-pow(rho1,2))

        # Sample y1 conditionally on z1
        y1.only[i] ~ dnorm(y1.only.mean[i], 1/y1.only.var[i])

        # Sample z2 for the studies that did NOT report first endpoint
        z2.missing[i] ~ dnorm(gamma02+gamma12/s2.imputed[i],1) T(,0)
      }

      # Prior for the theta's where both endpoints are reported
      theta.mean[1] <- mu1
      theta.mean[2] <- mu2
      thetacov.both[1,1] <- pow(tau1,2)
      thetacov.both[1,2] <- rho4*tau1*tau2
      thetacov.both[2,1] <- rho4*tau1*tau2
      thetacov.both[2,2] <- pow(tau2,2)

      # Prior for theta's for studies that reported BOTH endpoints
      for(i in 1:m1){
        theta.both[i,1:2] ~ dmnorm.vcov(theta.mean, thetacov.both)
      }

      # Prior for theta1 for studies that reported first endpoint ONLY
      for(i in 1:m2){
        theta1.only[i] ~ dnorm(mu1, 1/pow(tau1,2))
      }

      # Prior for mu1
      mu1 ~ dnorm(0,0.0001)

      # Prior for mu2
      mu2 ~ dnorm(0,0.0001)

      # Prior for tau1. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
      tau1 ~ dt(0,1,1) T(0,)

      # Prior for tau2. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
      tau2 ~ dt(0,1,1) T(0,)

      # Prior for gamma01
      gamma01 ~ dunif(-2,2)

      # Prior for gamma02
      gamma02 ~ dunif(-2,2)

      # Prior for gamma11
      gamma11 ~ dunif(0, max.s1)

      # Prior for gamma12
      gamma12 ~ dunif(0, max.s2)

      # Prior for rho1
      rho1 ~ dunif(-1,1)

      # Prior for rho2
      rho2 ~ dunif(-1,1)

      # Prior for rho3
      rho3 ~ dunif(-1,1)

      # Prior for rho4
      rho4 ~ dunif(-1,1)

    }"

    # Compile the model in JAGS
    model <- rjags::jags.model(textConnection(model_string),
                               data = list(y.both=y.both,s.both=s.both,
                                           y1.only=y1.only,s1.only=s1.only,
                                           s2.imputed=s2.imputed,
                                           m1=m1,m2=m2,
                                           max.s1=max.s1,max.s2=max.s2),
                               inits = list(init.vals1,init.vals2,init.vals3),
                               n.chains=3,
                               quiet=TRUE)

    # Draw samples
    # First use function 'update' to draw warm-up (burnin) samples
    stats::update(model, n.iter=burn, progress.bar="none") # Burnin for 10,000 samples

    # Next use the function 'coda.samples' to produce the next nmc
    # samples which will ultimately used to approximate the posterior.
    full.samples <- rjags::coda.samples(model,
                                        variable.names=c("mu1","mu2","tau1","tau2",
                                                         "gamma01","gamma02","gamma11","gamma12",
                                                         "rho1","rho2","rho3","rho4"),
                                        n.iter=nmc,
                                        progress.bar="none")

  }

  #############################################
  # Extract samples. Can be used for plotting #
  # and obtaining summary statistics.         #
  #############################################
  full.samples.mat <- as.matrix(do.call(rbind,full.samples))
  mu1.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="mu1")]
  mu2.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="mu2")]
  tau1.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="tau1")]
  tau2.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="tau2")]
  gamma01.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="gamma01")]
  gamma02.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="gamma02")]
  gamma11.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="gamma11")]
  gamma12.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="gamma12")]
  rho1.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="rho1")]
  rho2.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="rho2")]
  rho3.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="rho3")]
  rho4.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="rho4")]

  # Chains for mu1. To check convergence
  mu1.chains = matrix(mu1.samples, ncol=3)
  colnames(mu1.chains) = c("Chain1","Chain2","Chain3")

  # Chains for mu2. To check convergence
  mu2.chains = matrix(mu2.samples, ncol=3)
  colnames(mu2.chains) = c("Chain1","Chain2","Chain3")


  ###########################
  # Extract point estimates #
  ###########################
  mu1.hat <- mean(mu1.samples)
  mu2.hat <- mean(mu2.samples)
  tau1.hat <- mean(tau1.samples)
  tau2.hat <- mean(tau2.samples)
  gamma01.hat <- median(gamma01.samples)
  gamma02.hat <- median(gamma02.samples)
  gamma11.hat <- median(gamma11.samples)
  gamma12.hat <- median(gamma12.samples)
  rho1.hat <- median(rho1.samples)
  rho2.hat <- median(rho1.samples)
  rho3.hat <- median(rho3.samples)
  rho4.hat <- median(rho4.samples)


  ###############
  # Return list #
  ###############
  return(list(mu1.hat = mu1.hat,
              mu2.hat = mu2.hat,
              tau1.hat = tau1.hat,
              tau2.hat = tau2.hat,
              gamma01.hat = gamma01.hat,
              gamma02.hat = gamma02.hat,
              gamma11.hat = gamma11.hat,
              gamma12.hat = gamma12.hat,
              rho1.hat = rho1.hat,
              rho2.hat = rho2.hat,
              rho3.hat = rho3.hat,
              rho4.hat = rho4.hat,
              mu1.samples = mu1.samples,
              mu2.samples = mu2.samples,
              tau1.samples = tau1.samples,
              tau2.samples = tau2.samples,
              gamma01.samples = gamma01.samples,
              gamma02.samples = gamma02.samples,
              gamma11.samples = gamma11.samples,
              gamma12.samples = gamma12.samples,
              rho1.samples = rho1.samples,
              rho2.samples = rho2.samples,
              rho3.samples = rho3.samples,
              rho4.samples = rho4.samples,
              mu1.chains = mu1.chains,
              mu2.chains = mu2.chains
              ))
}

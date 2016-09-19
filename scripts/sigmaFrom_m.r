# Estimate sigma from effective density and isolation by distance slope
#
# De: point estimate of effective density
# Des: samples from the probability distribution of De values
# m: the estimate of the isolation by distance slope
# mse: the standard error on m
# dims: the number of dimensions (1 or 2)


sigmaFrom_m <- function(De, Des=rep(De,1000), m, mse, dims=1){
	niter <- length(Des)

	# generate range of slope estimates
	ms = rnorm(niter, mean = m, sd = mse)

	# calculate sigma
	if(dims==1){
		sigma_point <- sqrt(1/(4*De*m)) # the point estimate for 1D
		sigmas = sqrt(1/(4*Des*ms)) # the range of estimates for 1D
	}
	if(dims==2){
		sigma_point = sqrt(1/(4*pi*De*m)) # the point estimate for 2D
		sigmas = sqrt(1/(4*pi*Des*ms)) # the range of estimates for 2D
	}
	if(!(dims %in% c(1,2))){
		stop('Dimensions need to be 1 or 2')
	}

	return(list(sigma_point=sigma_point, sigmas=sigmas))
}
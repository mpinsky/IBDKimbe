# Estimate sigma from effective density and Migraine estimate of Wright's neighborhood size
#
# De: point estimate of effective density
# Des: samples from the probability distribution of De values


sigmaFrom_NbWright <- function(De, Des=rep(De,1000), NbWright, NbWright_l95, NbWright_u95, dims=1){
	niter <- length(Des)

	a <- c(log(NbWright/NbWright_l95), log(NbWright_u95/NbWright)) # should be about equal if this is a log-normal distribution
	NbWright_sdlog <- mean(a)/1.96
	NbWrights <- rlnorm(niter, mean = log(NbWright), sd = NbWright_sdlog)

	# some checks for debugging
	# print(a) # should be about equal
	# print(quantile(NbWrights, c(0.025, 0.975)))
	# print(c(NbWright_l95, NbWright_u95)) # should be somewhat close to values in previous line

	if(dims == 1){
		sigma_point = sqrt(NbWright/4/De)
		sigmas = sqrt(NbWrights/4/Des)
	} else {
		stop('Script does not handle 2D Migraine yet')
	}

	return(list(sigma_point=sigma_point, sigmas=sigmas))
}
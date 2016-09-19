# Calculate De (effective density) from census density (D) and effective to census population size ratio (Ne/N)
#
# niter: number of iterations to produce
# D: the estimate of census density
# Dse: the standard error on D
# nen: the estimate of the Ne/N ratio (effective population size to census population size)
# nense: the standard error on nen


DeFromD <- function(niter, D, Dse, nen, nense){
	# generate point values and range of values for De
	De <- D*nen # the point estimate for De
	Ds <- rnorm(niter, mean = D, sd = Dse) # the range of estimates for D
	nens <- rnorm(niter, mean=nen, sd = nense) # and for Ne/N
	Des <- Ds*nens # and for De

	return(list(De=De, Des=Des))

}
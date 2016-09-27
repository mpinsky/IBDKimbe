######################################################################
## Resampling approach to estimate 95% CI on dispersal kernel spread
######################################################################

# Read in scripts
source('scripts/DeFromD.R')
source('scripts/DeFromNb.R')
source('scripts/sigmaFrom_m.r')
source('scripts/sigmaFrom_NbWright.r')
source('scripts/summarizeSigmas.R')

# Set parameters
niter <- 1000000 # number of iterations
alpha <- 1.5 # minimum age of maturity (for Waples Nb to Ne calculations)
AL <- 13.1 # adult lifespan (for Waples Nb to Ne calculations)
A <- 225 # length of Kimbe Bay (km)
m <- 0.00014 # isolation by distance slope (genetic distance/km)
mse <- 2.00E-05 # standard error on m



# Estimate sigma (dispersal kernel spread) in various ways

	# Using Waples Nb (# breeders) and IBD slope
	# Waples_Nb value and CI from NeEstimator results
D1 <- DeFromNb(niter, AL, alpha, Waples_Nb=1426.5, Waples_Nbl95=960.5, Waples_Nbu95=2631, A=225) # estimate De
s1 <- sigmaFrom_m(De=D1$De, Des=D1$Des, m=m, mse=mse, dims=1) # estimate sigma
s1$sigma_point # report point value for sigma
summarizeSigmas(sigmas=s1$sigmas, sg=3) # report distribution for sigma

	# Using Waples Nb and Migraine estimate of Wright's neighborhood size
s2 <- sigmaFrom_NbWright(De=D1$De, Des=D1$Des, NbWright=7738, NbWright_l95=6806, NbWright_u95=10527, dims=1)
s2$sigma_point
summarizeSigmas(sigmas=s2$sigmas, sg=3)
			
	# Erroneous calculation using census density
D3 <- DeFromD(niter, D=188.1, Dse=45, nen=1, nense=0)
s3 <- sigmaFrom_m(De=D3$De, Des=D3$Des, m=m, mse=mse, dims=1)
s3$sigma_point
summarizeSigmas(sigmas=s3$sigmas, sg=3)
	
	# Erroneous calculation using 0.1% of census density
D4 <- DeFromD(niter, D=188.1, Dse=45, nen=0.001, nense=0)
s4 <- sigmaFrom_m(De=D4$De, Des=D4$Des, m=m, mse=mse, dims=1)
s4$sigma_point
summarizeSigmas(sigmas=s4$sigmas, sg=3)


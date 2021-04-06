# Calculate sigma (dispersal kernel spread) from number of breeders (Nb)
# Uses Waples et al. 2013 PRSB/Waples et al. 2014 Genetics correction for Nb to Ne (doi: 10.1534/genetics.114.164822)
#
# niter: number of iterations to produce
# AL: adult lifespan (years)
# alpha: minimum age of maturity (years)
# Waples_Nb: linkage disequilibrium estimate of the number of breeders (e.g., output from NeEstimator based on a single cohort)
# Waples_Nbl95: lower 95% confidence interval for Waples_Nb
# Waples_Nbu95: upper 95% confidence interval for Waples_Nb
# A: length (1D) or area (2D) of the region to which the Nb estimate applies


DeFromNb <- function(niter, AL, alpha, Waples_Nb, Waples_Nbl95, Waples_Nbu95, A){
	source('scripts/findChiSqCIs.R') # find paramaters for a chi-sq distribution, given the 95% CIs

	# Fit the Waples et al. model to relate log10(AL/alpha) to the Nb/Ne ratio
	dat = read.csv('data/AllSpecies4_trim.csv') # data from Waples et al. 2013 PRSB
	mod = lm(NbNe ~ logRLAM, data=dat) # fit model

	# Get CIs and samples for our Nb/Ne ratio (using AL and alpha)
	NbNe = predict(mod, new=data.frame(logRLAM = log10(AL/alpha)), interval='prediction', level=0.95) # get the 95%CI prediction intervals. This uses log10 to match the example calculations in the first full paragraph in the left column of p. 775 in Waples et al. 2014 Genetics
	a = c(NbNe[,1] - NbNe[,2], NbNe[,3] - NbNe[,1]) # the CI intervals
	NbNes1 = rnorm(niter, mean = NbNe[,1], sd = mean(a)/1.96) # generate first set of Nb/Ne samples
	NbNes2 = rnorm(niter, mean = NbNe[,1], sd = mean(a)/1.96) # generate second set (need two draws since used 2x in equation)

	# generate point values for De
	Waples_Nb_adj = Waples_Nb/(1.26-0.323*NbNe[,1]) # point estimate for Nb after correction. From Waples et al. 2014 Genetics. We use the equation for "True Nb/Ne" from Table 3 because we have estimated Nb/Ne as a probability distribution (the "Using two traits" version in Table 3 is nothing more than inserting Nb/Ne=0.485+0.758log(AL/alpha) into Eq. 8). 
	Waples_Ne_adj = Waples_Nb_adj/NbNe[,1] # point estimate for Ne after correction. From Eq. Waples et al. 2014. We used the equation for "True Nb/Ne" from Table 3 because we have estimated Nb/Ne as a probability distribution.
	De <- Waples_Ne_adj/A # the point estimate for De

	# generate samples from Nb
	a2 <- findChiSqCIs(Waples_Nb, Waples_Nbl95, Waples_Nbu95, verbose=FALSE) # find the degrees of freedom for a chi-squared distribution that fits our upper and lower 95% CI values
	WaplesNbs = a2*Waples_Nb/rchisq(niter, df=a2) # generate Nb values from chisq distribution (unadjusted Nb)

	# optional debugging to check fitting of the chi-squared distribution
	# print(paste('a2=', a2, 'for i=',i))
	# print(rbind(signif(quantile(WaplesNbs, c(0.025, 0.975)),4), signif(Waples_Nbl95, Waples_Nbu95))) # second line values should be close to values in first line

	# convert range of Nbs to range of Nes and Des using Waples et al. 2013 PRSB equations
	Nbs_adj = WaplesNbs/(1.26-0.323*NbNes1)
	Nes_adj = Nbs_adj/NbNes2
	Des <- Nes_adj/A # and for De
#	Del95 = quantile(Des, 0.025)
#	Deu95 = quantile(Des, 0.975)

	return(list(De=De, Des=Des))
	
}
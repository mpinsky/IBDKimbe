########################################################
# Script to calculate the isolation by distance slope
# Jackknife over populations
########################################################
library(vegan)
source('scripts/readArlFstBatch.R') # to read in an Arlequin Fst batch run file

# load data
fsts = readArlFstBatch('data/FSTdist.sum') # batch run over jackknife options
geo1D = read.csv("data/Dist_1DPCA.csv", row.names=1, header=T) # distance matrix
fst = fsts$"Apercula_noperc21,perc14_2015-06-05.arp" # pull out the run with all loci included
fstnames = c('Malu', 'Muli', 'Restorff', 'Shuman', 'Talele', 'Tarobi', 'Tiwongo', 'Tuare', 'WestChaimain', 'Wulai') # site names in Fst file
	names(fst) = fstnames
	rownames(fst) = fstnames

# order the matrices in the same way
ord = c('Tuare', 'WestChaimain', 'Malu', 'Restorff', 'Shuman', 'Wulai', 'Tarobi', 'Tiwongo', 'Muli', 'Talele') # order to use for matrices
fst = fst[ord, ord]
geo1D = geo1D[ord, ord]

# fill matrices across the diagonal where needed
for(i in 1:nrow(fst)){
	for(j in 1:ncol(fst)){
		if(is.na(fst[i,j])) fst[i,j] = fst[j,i]
		if(is.na(geo1D[i,j])) geo1D[i,j] = geo1D[j,i]
	}
}


# data.frame to hold results from jackknifing over populations
n <- rep(NA, nrow(fst))
jack2 <- data.frame(popdrop = n, int = n, intse = n, slope = n, slopese = n, p = n, r2 = n)

# step through each jackknife option
for(i in 1:length(n)){
	cat(i)

	# drop a population
	thisfst = fst[-(i), -(i)]
	thisdist <- geo1D[-(i), -(i)]

	# linearize fst, remove diagonals
	thisfst[upper.tri(thisfst, diag=T)] = NA
	thisfstlin = thisfst/(1-thisfst)

	x = as.dist(thisdist)[1:length(as.dist(thisdist))]
	y = as.dist(thisfstlin)[1:length(as.dist(thisfstlin))]
	mod=lm(y ~ x)
	mant = mantel(thisdist, thisfstlin, permutations=9999)

	jack2$popdrop[i] <- names(fst)[i]
	jack2$r2[i] <- summary(mod)$r.squared
	jack2$p[i] <- mant$signif
	jack2$int[i] <- mod$coef[1]
	jack2$intse[i] <- summary(mod)$coef[1,2]
	jack2$slope[i] <- mod$coef[2]
	jack2$slopese[i] <- summary(mod)$coef[2,2]

}


out <- data.frame(drop = jack2$popdrop, a = jack2$int, ase = jack2$intse, m = jack2$slope, mse = jack2$slopese)
out

min(out$m)
max(out$m)
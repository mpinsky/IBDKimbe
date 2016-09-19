########################################################
# Script to calculate the isolation by distance slope
# all loci, all populations
########################################################

library(vegan)
source('scripts/readArlFstBatch.R') # to read in an Arlequin Fst batch run file

# read in Fst and geographic distance matrices
geo1D = read.csv("data/Dist_1DPCA.csv", row.names=1, header=T) # geographic distance matrix
fsts = readArlFstBatch('data/FSTdist.sum') # genetic distance matrix: Arlequin batch run jackknifing over loci
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

# remove the upper triangle and diagonals from the matrices, then linearize Fst
geo1D[upper.tri(geo1D, diag=T)] = NA
fst[upper.tri(fst, diag=T)] = NA
fstlin = fst/(1-fst)

# Mantel test
mant = mantel(geo1D, fstlin, permutations=9999)
mant$signif # p-value for isolation by distance pattern

# Regression
x = as.dist(geo1D)[1:length(as.dist(geo1D))]
y = as.dist(fstlin)[1:length(as.dist(fstlin))]
mod=lm(y ~ x)

# slope and 95% CI on slope
summ <- summary(mod)
summ$coefficients[2,1] # slope
summ$coefficients[2,1]-summ$coefficients[2,2]*1.96 # lower 95% CI
summ$coefficients[2,1]+summ$coefficients[2,2]*1.96 # upper 95% CI

# calculate CIs for plotting
CIs = predict(mod, newdata=data.frame(x=x), interval='confidence')
CIs
#save(CIs, file='ibdCIs.csv')
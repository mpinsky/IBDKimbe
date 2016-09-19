########################################################
# Script to calculate the isolation by distance slope
# Jackknife over loci
########################################################
library(vegan)
source('readArlFstBatch.R') # to read in an Arlequin Fst batch run file

# load data
geo1D = read.csv("data/Dist_1DPCA.csv", row.names=1, header=T) # distance matrix
fsts = readArlFstBatch('data/FSTdist.sum') # Fsts: Arlequin batch run jackknifing over loci

fstnames = c('Malu', 'Muli', 'Restorff', 'Shuman', 'Talele', 'Tarobi', 'Tiwongo', 'Tuare', 'WestChaimain', 'Wulai') # site names in Fst file
ord = c('Tuare', 'WestChaimain', 'Malu', 'Restorff', 'Shuman', 'Wulai', 'Tarobi', 'Tiwongo', 'Muli', 'Talele') # order to use for matrices


# data.frame to hold results from jackknifing over loci
n <- rep(NA, length(fsts)-1) # -1 since we skip the first fst matrix (holds all loci)
jack <- data.frame(file = n, slope = n, slopese = n, p = n, r2 = n)

# iterate over each jackknife option
# skip the first one, since it has all loci included
for(i in 2:length(fsts)){
	cat(i)
	thisfst = fsts[[i]]
	names(thisfst) = fstnames
	rownames(thisfst) = fstnames
	thisfst = thisfst[ord, ord]

	# fill across diagonal where needed
	for(j in 1:nrow(thisfst)){
		for(k in 1:ncol(thisfst)){
			if(is.na(thisfst[j,k])) thisfst[j,k] = thisfst[k,j]
		}
	}

	# linearize fst, remove diagonals
	thisfst[upper.tri(thisfst, diag=T)] = NA
	thisfstlin = thisfst/(1-thisfst)

	# fit linear model
	x = as.dist(geo1D)[1:length(as.dist(geo1D))]
	y = as.dist(thisfstlin)[1:length(as.dist(thisfstlin))]
	mod=lm(y ~ x)
	
	# mantel test
	mant = mantel(geo1D, thisfstlin, permutations=9999)

	# save results
	jack$file[i-1] <- names(fsts)[i]
	jack$r2[i-1] <- summary(mod)$r.squared
	jack$p[i-1] <- mant$signif
	jack$slope[i-1] <- mod$coef[2]
	jack$slopese[i-1] <- summary(mod)$coef[2,2]

}

jack$droploci <- gsub('Apercula_noperc21,perc14|Apercula_noperc21,perc14,|_2015-06-05.arp', '', jack$file)


out <- data.frame(drop = jack$droploci, m = jack$slope, mse = jack$slopese)
out

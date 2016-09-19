# Script to reformat the genotype data into file formats used by analysis programs

###############################################
## Write Genepop format for all samples and all loci
## Used by genepop for testing HWP and LD
###############################################
source('scripts/writeGenepop.r')

genowide <- read.csv('data/Apercula_genotypes_2015-06-01.csv', header=TRUE)

loccols <- grep(".1", names(genowide), fixed=TRUE) # columns with genetic data
loctowrite <- unlist(strsplit(names(genowide)[loccols], ".1", fixed=TRUE)) # locus names

genpopname <- "output/Apercula.gen"
header <- 'all populations, no perc21'
genos <- genowide[,sort(c(loccols, loccols+1))] # the genotypes to write out

writeGenepop(file=genpopname, header=header, geno = genos, loci=loctowrite, pops=genowide$Island2)


###############################################
## Write Genepop format for New Recruits only
## Used by NeEstimator
## Drop perc21 and perc14 (out of HWP)
###############################################
source('scripts/writeGenepop.r')

genowide <- read.csv('data/Apercula_genotypes_2015-06-01.csv', header=TRUE)
genowideNR <- genowide[genowide$Stage=='NR',] # trim to new recruits
	dim(genowideNR) # 190 individuals

genowideNR$Island3 = "Onepop" # Group all new recruits together in one population	


loccols <- grep(".1", names(genowideNR), fixed=TRUE) # columns with genetic data
locnms <- unlist(strsplit(names(genowideNR)[loccols], ".1", fixed=TRUE)) # locus names

alwaysdrop <- c('perc21', 'perc14') # to drop the loci out of HWP
loctowrite <- setdiff(locnms, alwaysdrop) # names of loci to write out
loccolstowrite <- which(names(genowideNR) %in% paste(loctowrite, '.1', sep='')) # columns of loci to write

genpopname <- paste("Output/Apercula_NRonepop_no", paste(alwaysdrop, collapse = ","), ".gen", sep="") # output file name
header <- 'group all together' # file header
genos <- genowideNR[,sort(c(loccolstowrite, loccolstowrite+1))] # genotypes to write out

writeGenepop(file=genpopname, header=header, geno = genos, loci=loctowrite, pops=genowideNR$Island3) # write out



#################################################
## Arlequin files							    #
## All samples                                  #
## Jackknife over loci                          #
#################################################
source("scripts/matchall.R")

genowide <- read.csv('data/Apercula_genotypes_2015-06-01.csv', header=TRUE)

pops <- sort(unique(genowide$Island2)) # population names
loccols <- grep(".1", names(genowide), fixed=T) # columns with genotype data
locnms <- unlist(strsplit(names(genowide)[loccols], ".1", fixed=T)) # names of loci

# Write out each region to file, dropping particular loci
alwaysdrop <- c('perc21', 'perc14') # loci to drop in every file (out of HWP)
sometimesdrop <- setdiff(locnms, alwaysdrop) # loci to jackknife over
todrop <- vector('list',length(sometimesdrop)+1) # list to hold vectors of loci to write out (jackknifing over each one in turn)
todrop[[1]] <- alwaysdrop # first list item: keep all the loci (except the two out of HWP)
for(i in 1:length(sometimesdrop)){ # drop the loci out of HWP, plus jackknife over the remainder
	todrop[[i+1]] <- c(alwaysdrop, sometimesdrop[i])
} 

outfiles = character(0) # will hold a list of the files written out

# Write out each Arlequin file
for(d in 1:length(todrop)){
	print(d)
	if(length(todrop[[d]])>0){
		theseloci <- locnms[-matchall(todrop[[d]], locnms)] # the loci to write out, if dropping one or more
	} else {
		theseloci <- locnms # loci to write out if dropping none
	}
	
	thesecols = matchall(paste(theseloci, ".1", sep=""), names(genowide)) # columns for the loci to output in this round
		
	# Open file for writing
	if(length(todrop[[d]])>0){
		name = paste("Apercula_no", paste(todrop[[d]], collapse=","), ".arp", sep="") # file name
	} else {
		name = paste("Apercula.arp", sep="") # file name if dropping none
	}
	outfiles = c(outfiles, name) # keep track of the file written out
	name = paste("output/", name, sep="") # add directory to file names
	file = file(name) # create file handle
	open(file, "w") # open the file

	# Write the header
	cat("[Profile]\n", file= file, append=T)
	cat(paste("\tTitle=\"A. percula Kimbe Bay no ", paste(todrop[[d]], collapse=" or "), "\"\n", sep=""), file= file, append=T)
	cat(paste("\tNbSamples=", length(pops), "\n", sep=""), file= file, append=T)
	cat("\tGenotypicData=1\n", file= file, append=T)
	cat("\tGameticPhase=0\n", file= file, append=T)
	cat("\tRecessiveData=0\n", file= file, append=T)
	cat("\tDataType=STANDARD\n", file= file, append=T) # should this be msat?
	cat("\tLocusSeparator=WHITESPACE\n", file= file, append=T)
	cat("\tMissingData='?'\n", file= file, append=T)
	cat("\tCompDistMatrix=1\n", file= file, append=T)
	
	# Write the data header
	cat("[Data]\n", file=file, append=T)

	cat(paste("\t[[Samples]] # Data for ", length(thesecols), " loci: ", paste(theseloci, collapse = " "), "\n", sep=""), file=file, append=T)

	# Write the samples
	for(j in 1:length(pops)){
		k = which(genowide$Island2==pops[j])
		if(length(k)>0){
			cat(paste("\tSampleName=\"", pops[j], "\"\n", sep=""), file=file, append=T)
			cat(paste("\tSampleSize=", length(k), "\n", sep=""), file=file, append=T)
			cat("\tSampleData=  {\n", file=file, append=T)

			# Write each individual in this population
			for(i in 1:length(k)){
				# allele 1 on first line with sample name
				out = paste("\t", genowide$SampleID[k[i]], "\t1 ", paste(genowide[k[i],thesecols], collapse=" "), "\n", sep="") # format the genotypes for allele 1 of all loci with spaces between alleles
				out = gsub("NA", "?  ", out) # replace NAs with ?
				cat(out, file=file, append=T)
				
				# allele 2 on second line	
				out = paste("\t      \t  ", paste(genowide[k[i],thesecols+1], collapse=" "), "\n", sep="") # format the genotypes for allele 2 of all loci with spaces between alleles
				out = gsub("NA", "?  ", out) # replace NAs with ?
				cat(out, file=file, append=T)
			}
			cat("\t}\n", file=file, append=T)
		}
	}
	close(file)	
}


# write out an Arlequin batch file
name = paste("output/Aclarkii_jackknife_", Sys.Date(), ".arb", sep="")
file = file(name)
open(file, "w")
for(i in 1:length(outfiles)){
	cat(paste(outfiles[i], "\n", sep=""), file=file, append=T)
}
close(file)	
	
	
###############################################
## Genepop format for all samples
## Used by Migraine
###############################################
source('scripts/writeGenepop.r')

genowide <- read.csv('data/Apercula_genotypes_2015-06-01.csv', header=TRUE)
locs <- read.csv('data/Site UTM PCA aggreg.csv', header=TRUE, row.names=1) # x,y locations of sites (rotated to PCA axes)

loccols <- grep(".1", names(genowide), fixed=TRUE) # columns with genetic data
locnms <- unlist(strsplit(names(genowide)[loccols], ".1", fixed=TRUE)) # locus names

alwaysdrop <- c('perc21', 'perc14') # to drop the standard set (out of HWP)
loctowrite <- setdiff(locnms, alwaysdrop) # the names of loci to write
loccolstowrite <- which(names(genowide) %in% paste(loctowrite, '.1', sep='')) # the columns of loci to write

genpopname <- paste("Output/Apercula_migraine1D_no", paste(alwaysdrop, collapse = ","), ".gen", sep="") # file name
header <- 'all populations, no perc21 or perc14' # file header

# calc population coordinates in km (for Migraine)
locs$x <- round((locs$PC1_m - min(locs$PC1_m))/1000 + 100) # convert from m to km, re-center to start at 100
locs$y <- 2 # since 1D, set all y-coords to the same for Migraine

# set up individual names. last one in each population is the x and y coords (for Migraine)
nms <- rep('', nrow(genowide))
poplist <- sort(unique(genowide$Island2))
for(i in 1:length(poplist)){
	k <- which(genowide$Island2 == poplist[i])
	j <- which(locs$pop == poplist[i])
	nms[max(k)] <- paste(locs$x[j], locs$y[j])
}

writeGenepop(file=genpopname, header=header, geno = genowide[,sort(c(loccolstowrite, loccolstowrite+1))], loci=loctowrite, pops=genowide$Island2, nms = nms)



# Reads an Arlequin FSTdist.sum file and returns a list of matrices, each lower triangular
# of the pairwise FSTs between populations
# names of the list correspond to the filenames
# from 090311: added code to detect end of fsts

readArlFstBatch <- function (filename) {
	line = readLines(filename, n=-1)
	offset = 1
	while(line[offset] != "\tfile_name"){ offset = offset+1 } # iterate through the header
	offset = offset + 1

	fin = offset
	while(line[fin] != "" & fin < length(line)){ fin = fin+1 } # find the end of the fsts
	if(length(grep(".arp", line[fin], fixed=TRUE))==0){ # if we ended on a blank line, back up one
		fin = fin-1
	}

	fsts = vector("list", (fin-offset)) # set up list to hold fst matrices
	names = vector("list", (fin-offset))
	for(i in offset:fin){ # read in all files
		lin = strsplit(line[i], "\t")
		name = strsplit(lin[[1]][1], "/") # get the name of this file
		name = name[[1]][length(name[[1]])] # get the name of this file
		name = gsub("\"", "", name) # remove the trailing quotation mark
		names[[i+1-offset]] = name

		len = length(lin[[1]])
		f = as.numeric(lin[[1]][2:len]) # get Fsts, skipping the filename
		f = data.frame(toTri(f))
		fsts[[i+1-offset]] = f
	}
	names(fsts) = names
	return(fsts)
}


# convert a vector to a lower triangular matrix (very memory inefficient in current implementation!)
toTri <- function(x) {
	fulllen = length(x)/2 # max height and width of the final matrix
	out =  matrix(0, nrow = fulllen, ncol = fulllen)
	linelen = 1
	i = 1
	while(i<length(x)){
		out[linelen,1:linelen] = x[i:(i+linelen-1)]
		i = i+linelen
		linelen = linelen + 1	
	}
	out = out[1:(linelen-1), 1:(linelen-1)] # trim matrix to final dims
	out[upper.tri(out)] = NA # set uppper triangular to NA
	return(out)
}
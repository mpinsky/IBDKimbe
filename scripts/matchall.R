# returns indices into table for aall matches of the members of x
# not very efficient right now

matchall = function(x, table){
	ind = numeric(0)
	x= unique(x)
	for(i in 1:length(x)){
		this = match(table, x[i])
		thisind = which(!is.na(this))
		ind = c(ind, thisind)
	}
	ind = sort(ind)
	return(ind)
}
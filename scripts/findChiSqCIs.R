findChiSqCIs <- function(De, Del95, Deu95, verbose = FALSE){
	a2 = 105
	a1 = 100
	while (abs(mean(as.numeric(a1))-a2) > 1){ # inefficient search for a df that fits my CIs
		if(mean(as.numeric(a1)) > a2) a2 = a2 - 1
		if(mean(as.numeric(a1)) < a2) a2 = a2 + 1
		a1 = c(Del95, Deu95)*qchisq(c(0.975, 0.025), df=a2)/De
		if(verbose) print(paste(paste(round(mean(as.numeric(a1)),2), collapse=','), a2))
	}
	return(a2)

}
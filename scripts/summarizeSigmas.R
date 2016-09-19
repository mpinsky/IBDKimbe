# Print some useful statistics for a distribution of sigma values
#
# sigmas: the samples from the sigma distribution
# sg: the number of significant digits to report

summarizeSigmas <- function(sigmas, sg=3){

	sigma_mean = signif(mean(sigmas, na.rm=TRUE), sg)
	sigma_median = signif(median(sigmas, na.rm=TRUE), sg)
	sigma_sd = signif(sd(sigmas, na.rm=TRUE),3)
	sigma_l95 = signif(quantile(sigmas, 0.025, na.rm=TRUE), sg)
	sigma_u95 = signif(quantile(sigmas, 0.975, na.rm=TRUE), sg)

	cat(paste(
		'Mean:         ', sigma_mean, '\n',
		'Median:       ', sigma_median, '\n',
		'SD:           ', sigma_sd, '\n',
		'Lower 95% CI: ', sigma_l95, '\n',
		'Upper 95% CI: ', sigma_u95, sep=''))
		
}
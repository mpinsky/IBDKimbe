# IBDtools
Data and scripts for isolation by distance analyses, Pinsky et al. in review

Basic structure:
* /data: The data from the paper 
* /output: empty, but used by scripts for saving files
* /scripts: utility scripts
* 1_assemblegenotypes.R: output files in the right format for Genepop, Arlequin, NeEstimator, and Migraine
* 2_mantelFst.R: Calculate isolation by distance slope and Mantel test (after running Arlequin)
* 3_mantelFst_jackknifeloci.R: same as 2_, but jackknifing over loci (after running Arlequin)
* 4_mantelFst_jackknifepopulations.R: same as 2_, but jackknifing over populations (after running Arlequin)
* 5_estimate_sigma.r: combine the data together to estimate dispersal kernel spread (sigma)

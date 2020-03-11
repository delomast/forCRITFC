#' remove markers with either all no-call or no variation
#' also generates default error and dropout rate inputs for biallelic SNPs based
#' on what is typically used as input for SNPPIT
#' @param baseline dataframe input for baseline
#' @param mixture dataframe input for mixture
#' @param unsamp dataframe input for unsampledPops
#' @export
filterMarkers <- function(baseline, mixture = NULL, unsamp = NULL){

	mixPresent <- !is.null(mixture)
	unsampPresent <- !is.null(unsamp)

	if(mixPresent && any(colnames(baseline)[3:ncol(baseline)] != colnames(mixture)[2:ncol(mixture)])) stop(
		"marker columns in baseline and mixture do not have the same names")

	# filter markers
	alleleList <- list()
	to_remove <- c()
	for(i in seq(3, ncol(baseline) - 1, 2)){
		if(!any(!is.na(baseline[,i])) ){
			cat("removing marker", colnames(baseline)[i], "because the baseline is missing all genotypes\n")
			to_remove <- c(to_remove, i)
			next
		}
		if(mixPresent && !any(!is.na(mixture[,i-1])) ){
			cat("removing marker", colnames(mixture)[i-1], "because the mixture is missing all genotypes\n")
			to_remove <- c(to_remove, i)
			next
		}
		if(unsampPresent){
			if(!any(!is.na(unsamp[,i])) ){
				cat("removing marker", colnames(unsamp)[i], "because unsamp is missing all genotypes\n")
				to_remove <- c(to_remove, i)
				next
			}
		}
		alleles <- unique(c(baseline[,i], baseline[,i+1]))
		if(mixPresent) alleles <- unique(c(alleles, mixture[,i-1], mixture[,i]))
		if(unsampPresent)	alleles <- unique(c(alleles, unsamp[,i], unsamp[,i+1]))
		alleles <- alleles[!is.na(alleles)]
		if(length(alleles) < 2){
			cat("removing marker", colnames(baseline)[i], "because there is no variation at this locus\n")
			to_remove <- c(to_remove, i)
			next
		} else if (length(alleles) > 2){
			stop("More than 2 alleles at locus", colnames(baseline[i], ". This function is only for ",
																		 "biallelic markers."))
		}
		alleleList[[colnames(baseline)[i]]] <- alleles
	}
	baseline <- baseline[,-c(to_remove, (to_remove + 1))]
	if(mixPresent) mixture <- mixture[,-c((to_remove - 1), to_remove)]
	if(unsampPresent) unsamp <- unsamp[,-c(to_remove, (to_remove + 1))]



	return(list(baseline = baseline,
					mixture = mixture,
					unsampledPops = unsamp
			))
}

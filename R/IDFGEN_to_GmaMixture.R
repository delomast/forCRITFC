#' create the mixture dataframe needed by createGmaInput from IDFGEN objects
#' @param mixturePops either a PopList or a Population
#' @param Markers the marker names (SNPs only) to use, in order of how they are used
#' @export
IDFGEN_to_GmaMixture <- function(mixturePops, Markers){

	if(class(mixturePops) == "Population"){
		mixturePops <- as.PopList(names(mixturePops))
	}
	if(class(mixturePops) != "PopList") stop("mixturePops must be either a Population or a PopList class")

	out <- data.frame()
	for(pop in mixturePops){
		tempInds <- inds(get(pop))
		outTemp <- data.frame(IndividualID = tempInds, stringsAsFactors = FALSE)
		genos <- scores(get(pop))[tempInds,Markers, drop = FALSE]
		for(m in Markers){
			tempGenos <- cbind(substr(genos[,m],1,1), substr(genos[,m],2,2))
			tempGenos[tempGenos[,1] == "0",1] <- NA # "0" as missing data
			tempGenos[tempGenos[,2] == "0",2] <- NA
			colnames(tempGenos) <- c(m, paste0(m, ".A2"))
			outTemp <- cbind(outTemp, tempGenos, stringsAsFactors = FALSE)
		}
		out <- rbind(out, outTemp, stringsAsFactors = FALSE)
	}
	rownames(out) <- NULL
	return(out)
}

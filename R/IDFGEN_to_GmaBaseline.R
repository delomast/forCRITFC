#' create the baseline dataframe needed by createGmaInput from IDFGEN objects
#' @param baselinePops either a PopList or a Population
#' @param Markers the marker names (SNPs only) to use, in order of how they are used
#' @export
IDFGEN_to_GmaBaseline <- function(baselinePops, Markers){

	if(class(baselinePops) == "Population"){
		baselinePops <- as.PopList(names(baselinePops))
	}
	if(class(baselinePops) != "PopList") stop("baselinePops must be either a Population or a PopList class")

	out <- data.frame()
	for(pop in baselinePops){
		tempInds <- inds(get(pop))
		outTemp <- data.frame(Population = pop, IndividualID = tempInds, stringsAsFactors = FALSE)
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

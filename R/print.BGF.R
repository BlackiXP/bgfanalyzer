#' Print a Biogas Fermentation
#'
#' Prints a BGF object to the R console
#'
#' Default print method for a BGF.
#'
#' @param x An object of class BGF
#' @inheritParams base::print
#' @returns print's details of a BGF object to the console
#'
#'
#' @export

print.BGF<-function(x,...){

    l1=sprintf(gettext("'%s' - a BGF with %d fermentation(s)"),x$ExpParam$name,nrow(x$metaData))

    l2=sprintf(gettext("$ExpParam: %d experimental parameters"),length(x$ExpParam))

    l3=sprintf(gettext("$metaData: %d meta variables"),length(colnames(x$metaData)))

    l4=sprintf(gettext("$BioGasData: %d observations of %d fermentation variables"),nrow(x$BioGasData),length(colnames(x$BioGasData)))

  if(isTRUE(all(is.na(x$BioGasData$yield)))){
    l5=gettext("A '$yield' is not calculated yet")
    cat(l1,"\n",l2,l3,l4,"\n",l5,sep="\n")

  }else{
    cat(l1,"\n",l2,l3,l4,"\n",sep="\n")
    print.data.frame(get_yield_summary(x,FALSE,FALSE))
  }
  invisible(x)
  return(0)

}


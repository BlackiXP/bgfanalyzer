#' Split or merge `BGF`'s
#'
#' Functions that allow sub-setting and merging of `BGF``'s.
#'
#' The function `subset_BGF` can be used to subset a `BGF`.
#' The user can either use reactor/ fermentation names (= row names of `metaData`-layer, values in column 'reactor' of `BioGasData`-layer), or reactor layouts (= values in 'Layout' column of `metaData`-layer) as a basis for subsetting.
#' Furthermore, it is possible to give a new name to the resulting subset.
#'
#'
#' @param x,y a `BGF`
#' @param reactor defaults to `NULL`. Can be an `integer` or `character` specifying the reactors/ fermentations to keep based on the reactor name (= row name of the fermentation in `metaData`-layer)
#' @param layout defaults to `NULL`.Can be an `integer` or `character` specifying the reactors/ fermentations to keep based on the reactor layout (= respective value in the 'Layout' column  of the `metaData`-layer)
#' @param name defaults to `NULL`. Can be a `character` providing a new name for the `BGF`
#'
#' @returns a `BGF`
#'
#' @examples
#' # create an example BGF
#' myBGF <- from_AMPTSV2_report(
#'         ReactorLayout = c("2*Blank","Cellulose","2*S1 ctrl","2*S1 7d","2*S1 4d",
#'                           "2*S2 ctrl","2*S2 4d","2*S2 6d"),
#'         BlankLabel = "Blank",
#'         name = "myBGF",
#'         InocToSubRatio=2,
#'         ProcessTemp = 52,
#'         path = base::system.file("extdata","AMPTSV2.csv",package = "bgfanalyzer"))
#' myBGF
#'
#' # subset using 'reactor' argument
#' myBGF_subset1 <- subset_BGF(myBGF,reactor=c(4:9),name="subset 1")
#' myBGF_subset1
#'
#' # subset using 'layout' argument
#' myBGF_subset2 <- subset_BGF(myBGF,layout=c("S2 ctrl","S2 4d","S2 6d"),name="subset 2")
#' myBGF_subset2
#'
#' @export
#'

# subset_BGF() ####
subset_BGF=function(x,reactor=NULL,layout=NULL,name=NULL){

  validate_BGF(x)

  if(isFALSE(is.null(reactor))){
    x$metaData<-x$metaData[reactor,]
    x$metaData$Layout<-factor(x$metaData$Layout,levels = unique(x$metaData$Layout))

    set<-rownames(x$metaData)

    x$BioGasData <- subset(x$BioGasData,x$BioGasData$reactor%in%set)
    x$BioGasData$reactor<-factor(x$BioGasData$reactor,levels = rownames(x$metaData))

  }

  if(isFALSE(is.null(layout))){
    x$metaData<-subset(x$metaData,x$metaData$Layout%in%layout)
    x$metaData$Layout<-factor(x$metaData$Layout,levels = unique(x$metaData$Layout))

    set<-rownames(x$metaData)

    x$BioGasData <- subset(x$BioGasData,x$BioGasData$reactor%in%set)
    x$BioGasData$reactor<-factor(x$BioGasData$reactor,levels = rownames(x$metaData))

  }

  if(isFALSE(is.null(name))) x$ExpParam$name=name

  x <- update_BGF(x)

  return(x)
}

#' @rdname subset_BGF
#'
#' @details
#' The function `merge_BGF` can be used to merge two `BGF`'s.
#' On the way, it will generate new row names for the `metaData`-layer and adjust the values in the 'reactor' column of the `BioGasData`-layer appropriately.
#'
#' In general, it is advised to only merge `BGF`'s that have the same values in their `ExpParam`-layer, as with the default setting only the `ExpParam`-layer of the first `BGF` passed to `merge_BGF` is preserved.
#' It is possible to change the default behavior with the 'mergeParam' argument.
#' In that case, new names for entries in the merged `ExpParam`-layer will be generated if a entry is either present in only one `ÈxpPram`-layer of the `BGF`'s or if the values for the same entry in both `BGF`'s differs.
#' The new names are prolonged by a specific tag (controlled via arguments 'x_tag' and 'y_tag'), so that it is clear where that entry originates from.
#'
#' Furthermore, these tags will be added to the `metaData`- and `BioGasData`-layer.
#'
#' @param mergeParam `logic`; defaults to `FALSE`. Should the `ExpParam`-layer of the `BGF`'s be also merged?
#' @param x_tag,y_tag defaults to `NULL`. Can be a `character` which is used as a tag when merging `ExpParam`-layers
#'
#' @examples
#' # merge the BGF's
#' mergedBGF <- merge_BGF(myBGF_subset1,myBGF_subset2,TRUE)
#' mergedBGF
#'
#' @export
#'

# merge_BGF ####
merge_BGF=function(x,y,mergeParam=FALSE,x_tag=NULL,y_tag=NULL,name=NULL){
  validate_BGF(x)

  validate_BGF(y)
  x_in=x

  if(isTRUE(is.null(x_tag))) x_tag=x$ExpParam$name
  if(isTRUE(is.null(y_tag))) y_tag=y$ExpParam$name

  # adjust rownames in metaData- and 'reactor' column in BioGasData-layer of 'x' and 'y'
  new_reactor=paste0("R",c(1:nrow(x$metaData)))
  names(new_reactor)<-rownames(x$metaData)

  x$BioGasData[,"reactor"]=as.character(x$BioGasData[,"reactor"])

  for(i in c(1:nrow(x$BioGasData))) x$BioGasData[i,"reactor"]=as.character(new_reactor[x$BioGasData[i,"reactor"]])

  rownames(x$metaData)=new_reactor

  new_reactor<-paste0("R",c(1:nrow(x$metaData)+nrow(y$metaData)))
  names(new_reactor) <- rownames(y$metaData)

  y$BioGasData[,"reactor"]=as.character(y$BioGasData[,"reactor"])

  for(i in c(1:nrow(y$BioGasData))) y$BioGasData[i,"reactor"]=as.character(new_reactor[y$BioGasData[i,"reactor"]])

  rownames(y$metaData)=new_reactor

  # adjust colnames in metaData- and 'reactor' column in BioGasData-layer of 'x' and 'y'
  x_meta <- names(x$metaData)
  y_meta <- names(y$metaData)

  x_BG <- names(x$BioGasData)
  y_BG <- names(y$BioGasData)

  missing_in_x_meta <- subset(y_meta,!y_meta%in%x_meta)
  missing_in_y_meta <- subset(x_meta,!x_meta%in%y_meta)

  missing_in_x_BG <- subset(y_BG,!y_BG%in%x_BG)
  missing_in_y_BG <- subset(x_BG,!x_BG%in%y_BG)

  for(i in missing_in_x_meta) x$metaData[,i]=NA
  for(i in missing_in_y_meta) y$metaData[,i]=NA

  for(i in missing_in_x_BG) x$BioGasData[,i]=NA
  for(i in missing_in_y_BG) y$BioGasData[,i]=NA


  # merge metaData- and BioGasData-layer
  y$metaData<-y$metaData[names(x$metaData)]
  y$BioGasData<-y$BioGasData[names(x$BioGasData)]

  x$metaData<-rbind(x$metaData,y$metaData)
  x$BioGasData<-rbind(x$BioGasData,y$BioGasData)

  # merge ExpParam-layers
  if(isTRUE(mergeParam)){

    x_params <- names(x$ExpParam)
    y_params <- names(y$ExpParam)

    match_params <- subset(y_params,y_params%in%x_params)

    match_l <- NULL

    for(i in match_params) match_l[i]=(x$ExpParam[[i]]==y$ExpParam[[i]])

    mismatch<-subset(match_l,match_l==FALSE)
    missing_param_in_x<-subset(y_params,!y_params%in%x_params)
    missing_param_in_y<-subset(x_params,!x_params%in%y_params)

    if(length(mismatch)>0 | length(missing_param_in_x)>0 | length(missing_param_in_y)>0){
      meta_param <- c(rep(x_tag,nrow(x_in$metaData)),rep(y_tag,nrow(y$metaData)))

      x <- add_metaData(x,meta_param,"MergeGroup")

      x$BioGasData[,names(x$metaData)[length(x$metaData)]]=NA

      for(i in c(1:nrow(x$BioGasData))) x$BioGasData$MergeGroup[i]=x$metaData[x$BioGasData$reactor[i],"MergeGroup"]

    }

    if(isFALSE(length(mismatch)==0)){
      for(j in c(1:length(mismatch))){
      nr<-grep(names(mismatch)[j],names(x$ExpParam))
      x$ExpParam[[paste0(names(x$ExpParam)[nr],"_",y_tag)]]=y$ExpParam[[nr]]
      names(x$ExpParam)[nr]<-paste0(names(x$ExpParam)[nr],"_",x_tag)

      }
    }


    if(isFALSE(length(missing_param_in_x)==0)){
     for(j in missing_param_in_x) x$ExpParam[[paste0(j,"_",y_tag)]]=y$ExpParam[[j]]
    }

    if(isFALSE(length(missing_param_in_y)==0)){
      for(j in missing_param_in_y){
        nr<-grep(j,names(x$ExpParam))
        names(x$ExpParam)[nr]<-paste0(names(x$ExpParam)[nr],"_",x_tag)

      }
    }

    if(isTRUE(is.null(name))) name="mergedBGF"

    x$ExpParam[["name"]] <- name
  }

  return(x)
}

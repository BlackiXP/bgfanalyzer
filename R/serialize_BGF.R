#' Transform a continious BGF into a serialized BGF by time
#'
#' A function that allows to split the data of a one fermentation BGF into a serialized form.
#' It splits the BGF's BioGasData after a given hydraulic retention time (HRT) so that after each HRT a new replicate begins.
#' New replicates are added to metaData and it is possible to recalculate the yield
#'
#' @param x a `BGF`
#' @param hrt `numeric` - e.g. the hydraulic retention time or any other time value to be used to split `BioGasData`
#' @param metaData description
#' @param name description
#' @param MeasurmentType description
#' @param calc_yield description
#' @inheritDotParams bgfanalyzer::netGasGC
#'
#' @export


serialize_by_HRT=function(x,hrt,keep_metaData=TRUE,name=NULL,MeasurementType=NULL,calc_yield=NULL,...){
  df<-x$BioGasData
  df$xHRT <- as.numeric(
    substr(
      as.character(
        as.numeric(df$time)/hrt+1
      ),1,1))

  for(i in c(2:max(df$xHRT))) {
    df$time[which(df$xHRT==i)]=df$time[which(df$xHRT==i)]-min(df$time[which(df$xHRT==i)])
    df$reactor<-factor(df$reactor,levels = c(unique(as.character(df$reactor)),paste0("R",i)))
    df$reactor[which(df$xHRT==i)]=paste0("R",i)
    df$product[which(df$xHRT==i)]=df$product[which(df$xHRT==i)]-df$product[which(df$xHRT==i)][1]

  }

  df<- df[,-which(names(df)=="xHRT")]

  x$BioGasData <- df

  x <- bgfanalyzer::update_BGF(x)

  x <- bgfanalyzer::calculate_flow_from_volume(x)

  x <- bgfanalyzer::relative_production(x)

  if(isTRUE(keep_metaData)) {
    x$metaData[c(2:nrow(x$metaData)),]=x$metaData[1,]
    x <- netGasGC(x,...)
    if(is.character(calc_yield)) calc_yield <- which(names(x$metaData)==calc_yield)
    if(is.numeric(calc_yield)) x <- bgfanalyzer::calc_yield(x,calc_yield)
  }else {
    x$metaData[c(1:nrow(x$metaData)),c(1:3)] <- x$metaData[1,c(1:3)]
    x$metaData <- x$metaData[,c(1:3)]
    if(isFALSE(is.null(calc_yield))){
      stopifnot(is.numeric(calc_yield))
      x <- bgfanalyzer::add_metaData(x,calc_yield,lab = "yield_base")
      x <- bgfanalyzer::calc_yield(x,"yield_base")
      x <- bgfanalyzer::summarize_yield(x)
    }
  }

  if(isFALSE(is.null(name))) x$ExpParam$name <- name
  if(isFALSE(is.null(MeasurementType))) x$ExpParam$MeasurementType <- MeasurementType

  return(x)
}

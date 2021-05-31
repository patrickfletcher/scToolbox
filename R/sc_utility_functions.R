# utility functions for Seurat
# Patrick Fletcher 2021

#TODO: if useful, contribute to Seurat?
#TODO: adhere to Seurat's style and coding conventions
#TODO: proper docs so that one can do ?function_name

#threshold computations

#############################
#' calculateGroupedThr

# Ideas: could store the result back into Seurat object's "misc" list
#
# TODO: group.by=NULL -> use all cells if not supplied
# TODO: default kup/kdown as NULL, but force user to supply at least one of them

#' a function to compute thresholds on feature values in a Seurat object.
#' 
#' @details 
#' thr_hi -> center(x)+kup*spread(x)
#' thr_lo -> center(x)-kdown*spread(x)
#' where x is a vector of values, center(x) is a measure of central tendency, and spread(x) is a
#' measure of spread.
#'
#' @param object A Seurat object
#' @param var A valid variable to compute a threshold on (must work with FetchData)
#' @param group.by A valid metadata variable to group cells with (must work with FetchData)
#' @param centerfun = "median" function to use as measure of central tendancy
#' @param spreadfun = "mad" function to use as measure of spread
#' @param kdown = NULL Number of units of spread below center for thr_lo
#' @param kup = NULL  Number of units of spread above center for thr_hi
#' @param clamplo = NULL Force minimum possible value for thr_lo
#' @param clamphi = NULL Force maximum possible value for thr_hi
#' @param omitZeros = FALSE Omit zero values before computing stats?
#' @return A dataframe containing the statistics computed and the thresholds computed from them
#' @export
calculateGroupedThr <- function(
  object, 
  var, 
  group.by, 
  centerfun = "median",
  spreadfun = "mad", 
  kdown = NULL,
  kup = NULL,
  clamplo = NULL,
  clamphi = NULL,
  omitZeros = FALSE
) {
  
  if (is.null(kdown) && is.null(kup)) {
    stop('At least one of kup or kdown must be and argument. Stopping...', call. = TRUE)
  }
  
  #extract the data and grouping variable **no error checking, uses active.assay
  df <- FetchData(object = object, c(var,group.by)) #FetchData returns data.frame
  
  if (omitZeros) {
    df <- df[which(df[,1]>0), ] #make an option to omit zeros?
  }
  
  data <-c(df[,1])
  grp <-list(c(df[,2])) #put grp values into a list of size 1 for aggregate
  
  mu_df <- aggregate(x=data, by=grp, centerfun) #returns data.frame with rownames=group.by names
  sig_df <- aggregate(x=data, by=grp, spreadfun)
  
  # also get the max/min
  max_df <- aggregate(x=data, by=grp, max)
  min_df <- aggregate(x=data, by=grp, min)
  
  max_val <- max_df[,2]
  min_val <- min_df[,2]
  
  levs <- mu_df[,1]
  mu <- mu_df[,2]
  sig <- sig_df[,2]
  
  #return a data.frame with rows=group levels, stats, and threshold.
  thr_df <- data.frame(mu, sig)
  rownames(thr_df)<-levs
  colnames(thr_df)<-c(centerfun, spreadfun)
  
  #now append the threshold column(s)
  if (!is.null(kdown)) {
    thr_lo <- mu - kdown*sig
    #clamp the thresholds to data range
    thr_lo = pmax(thr_lo, min_val)
    #finally apply user specified clamp values if requested
    if (!is.null(clamplo)) {
      thr_lo = pmax(thr_lo, clamplo)
    }
    thr_df$thr_lo <- thr_lo
  }
  
  if (!is.null(kup)) {
    thr_hi <- mu + kup*sig
    thr_hi = pmin(thr_hi, max_val)
    if (!is.null(clamphi)) {
      thr_hi = pmin(thr_hi, clamphi)
    }
    thr_df$thr_hi <- thr_hi
  }
  
  return(thr_df)
}

########################
#' lcalculateGroupedThr
#' 
#' A function that list-applies [calculateGroupedThr()], then combines results into a single data.frame
#' 
#' @param objectList a list of Seurat objects
#' @param var a list of features compatible with [FetchData()]
lcalculateGroupedThr <- function(
  objectList,
  var, 
  group.by, 
  ...
) {
  
  thr.list <- lapply(X=objectList, FUN = calculateGroupedThr, var=var, group.by=group.by, ...)
  
  # TODO: Should make sure row names are unique, in case group names are reused by objects (e.g., cluster ID)
  # extract all the rownames from each result, check if unique, if not prepend the object@project.name
  
  #collapse the list into a single dataframe. 
  thr.prctMT <- do.call(rbind, unname(thr.list))
}


#########################################
#' applyThresholds
#'
#' A function to subset cells of a Seurat object based on thresholds computed with calculateGroupedThr.
#' 
#' object,          Seurat object
#' featureThrList,
#'
#' Can specify either or both of thr_lo and thr_hi. This method will always keep cells with
#' feature>thr_lo and feature <thr_hi.
applyThresholds <- function(
  object, 
  featureThrList, 
  group.by, 
  rejectEquality=TRUE
) {
  
  group_df <- FetchData(object = object, vars = group.by)
  group <- as.factor(group_df[,1])
  groupnames <- levels(group)
  
  fnames <- names(featureThrList)
  
  #extract the relevant data:
  celldata <- FetchData(object = object, vars = names(featureThrList))
  
  rejected <- rep(FALSE, length(group))
  #extract the conditions for each group. Begin with rejected==FALSE, reject if any condition fails (use |)
  for (gname in groupnames){
    is_gname <- which(group == gname) #get this subset
    reject_gname <- rep(FALSE, length(is_gname))
    for (fname in fnames){
      thr_df <- featureThrList[[fname]]
      if ('thr_lo' %in% colnames(thr_df)){
        reject_gname <- reject_gname | celldata[is_gname,fname] < thr_df[gname,'thr_lo']  #< or <= ???
      }
      if ('thr_hi' %in% colnames(thr_df)){
        reject_gname <- reject_gname | celldata[is_gname,fname] > thr_df[gname,'thr_hi'] 
      }
    }
    rejected[is_gname] <- rejected[is_gname] | reject_gname
  }
  
  #subset the object
  object <- object[, which(x = rejected == FALSE)]
  return(object)
}


###################
#plotting functions

library(patchwork)

#########################################
# VlnPlotThr: violinplots with thresholds
# 
# object,         Seurat object
# featureThrList, named list: name=feature, value=dataframe of thresholds for that feature, computed by calculateGroupedThr
# group.by,       Same group.by as used to compute thresholds
# lwdth=2/3,
# hicol='red',
# locol='blue',
# segsize=1,
# ...             extra args for VlnPlot
#
VlnPlotThr <- function(
  object, 
  featureThrList,
  group.by,
  lwdth=2/3,
  hicol='red',
  locol='blue',
  segsize=1,
  ... 
) {
  
  group_df <- FetchData(object = object, vars = group.by)
  group <- as.factor(group_df[,1])
  groupnames <- levels(group)
  
  #check that groupnames are in threshold df
  if (!all(groupnames %in% rownames(featureThrList[[1]]))) {
    stop('group.by levels do not match threshold table - use same group.by as was used in calculateGroupedThr. Stopping...', call. = TRUE)
  }
  
  fnames <- names(featureThrList)
  
  plots <- VlnPlot(object = object, features = fnames, group.by=group.by, ...)
  
  #add the thresholds as horizontal bars on top of the violins
  for (f in 1:length(fnames)) {
    fname <- fnames[f]
    thr_df <- featureThrList[[fname]]
    for (i in 1:length(groupnames)) {
      gname <- groupnames[i]
      if ('thr_lo' %in% colnames(thr_df)){
        yval <- thr_df[gname, 'thr_lo']
        plots[[f]] <- plots[[f]] + geom_segment(
                                      x = i-lwdth/2, y = yval, xend = i+lwdth/2, yend = yval, 
                                      colour = locol, size = segsize)
      }
      if ('thr_hi' %in% colnames(thr_df)){
        yval <- thr_df[gname, 'thr_hi']
        plots[[f]] <- plots[[f]] + geom_segment(
                                      x = i-lwdth/2, y = yval, xend = i+lwdth/2, yend = yval, 
                                      colour = hicol, size = segsize)
      }
    }
  }
  return(plots)
}

#########################################
# FeatureScatterThr: FeatureScatter with thresholds
# 
# object,         Seurat object
# feature1,       feature1, must be a feature included in featureThrList 
# feature2,       feature2, must be a feature included in featureThrList
# featureThrList, named list: name=feature, value=dataframe of thresholds for that feature, computed by calculateGroupedThr
# group.by,       Same group.by as used to compute thresholds
# lwdth=2/3,
# hicol='red',
# locol='blue',
# segsize=1,
# ...             extra args for VlnPlot
#
# need to make separate plots for each group in group.by otherwise it would be cluttered
FeatureScatterThr <- function(
  object,
  feature1,
  feature2,
  featureThrList,
  group.by,
  lwdth=2/3,
  hicol='red',
  locol='blue',
  size=1,
  ...
) {
  
  #check that feature1 is in featureThrList
  if (!(feature1 %in% names(featureThrList))) {
    stop('group.by levels do not match threshold table - use same group.by as was used in calculateGroupedThr. Stopping...', call. = TRUE)
  }
  #check that feature2 is in featureThrList
  if (!(feature2 %in% names(featureThrList))) {
    stop('group.by levels do not match threshold table - use same group.by as was used in calculateGroupedThr. Stopping...', call. = TRUE)
  }
  
  group_df <- FetchData(object = object, vars = group.by)
  group <- as.factor(group_df[,1])
  groupnames <- levels(group)
  
  #check that groupnames are in threshold df
  if (!all(groupnames %in% rownames(featureThrList[[1]]))) {
    stop('group.by levels do not match threshold table - use same group.by as was used in calculateGroupedThr. Stopping...', call. = TRUE)
  }
  
  # make a separate plot for each group
  plots <- lapply(X=groupnames, FUN = function(gname) {
    
    group.cells = colnames(p.cells$p5)[group==gname]
    
    plot <- FeatureScatter(object = object, feature1 = feature1, feature2 = feature2, 
                               cells=group.cells, plot.cor=FALSE,...)
    
    #add the thresholds as horizontal/vertical lines
    #feature 1 (vertical thresholds)
    thr_df <- featureThrList[[feature1]]
    if ('thr_lo' %in% colnames(thr_df)){
      xval <- thr_df[gname, 'thr_lo']
      plot <- plot + geom_vline( xintercept = xval, colour = locol, size = size)
    }
    if ('thr_hi' %in% colnames(thr_df)){
      xval <- thr_df[gname, 'thr_hi']
      plot <- plot + geom_vline(xintercept = xval, colour = hicol, size = size)
    }
    
    #feature 2
    thr_df <- featureThrList[[feature2]]
    if ('thr_lo' %in% colnames(thr_df)){
      yval <- thr_df[gname, 'thr_lo']
      plot <- plot + geom_hline( yintercept = yval, colour = locol, size = size)
    }
    if ('thr_hi' %in% colnames(thr_df)){
      yval <- thr_df[gname, 'thr_hi']
      plot <- plot + geom_hline(yintercept = yval, colour = hicol, size = size)
    }
    # plot <- plot + ggtitle(gname) + theme(legend.position = "none")
    
  } )
  
  plots <- wrap_plots(plots, ncol = length(x = groupnames))
  
  return(plots)
}
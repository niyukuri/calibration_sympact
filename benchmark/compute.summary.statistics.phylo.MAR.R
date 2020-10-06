#' Computing phylogenetic realted summary statistics
#'
#' @param simpact.trans.net Transmission network and record produced by \code{\link{advanced.transmission.network.builder()}}
#' @param datalist.agemix Data list of simpact output produced by \code{\link{readthedata()}}
#' @param work.dir Working directory
#' @param dirfasttree Directory where fastTree soaftware is called from
#' @param sub.dir.rename Sub-directory where simpact output are stored
#' @param limitTransmEvents Consider a transmission network which counts individuals more than the value (numeric)
#' @param seq.cov Sequence coverage
#' @param seq.gender.ratio Proportion of women in the selected population (women/(women + men))
#' @param age.group.15.25 Consider individuals with age greater than 15 and less than 25
#' @param age.group.25.40 Consider individuals with age greater than 25 and less than 40
#' @param age.group.40.50 Consider individuals with age greater than 40 and less than 50
#' @param endpoint Only transmission events that took place before this point in simulation time, are captured in the transmission network
#' @param timewindow Time window in which the experience is carried out
#' @param cut.off Cut off value for constructing pairings based on tMRCA
#' @export





compute.summary.statistics.phylo.MAR <- function(simpact.trans.net = simpact.trans.net,
                                                 datalist.agemix = datalist.agemix,
                                                 work.dir = work.dir,
                                                 sub.dir.rename = sub.dir.rename,
                                                 dirfasttree = work.dir,
                                                 limitTransmEvents = 7,
                                                 seq.cov = 100,
                                                 seq.gender.ratio = 0.7,
                                                 age.group.15.25 = c(15,25),
                                                 age.group.25.40 = c(25,40),
                                                 age.group.40.50 = c(40,50),
                                                 endpoint = 40,
                                                 timewindow = c(30,40),
                                                 cut.off = 7){
  
  
  
  
  dirseqgen <- work.dir
  
  
  # Function for linear mixed effect models in transmission clusters
  ###################################################################
  
  
  mAr.IDs <- IDs.Seq.Age.Groups(simpact.trans.net = simpact.trans.net,
                                limitTransmEvents = limitTransmEvents,
                                timewindow = timewindow,
                                seq.cov = seq.cov,
                                seq.gender.ratio = seq.gender.ratio,
                                age.group.15.25 = age.group.15.25,
                                age.group.25.40 = age.group.25.40,
                                age.group.40.50 = age.group.40.50)
  
  # nrow(agemixing.df.IDs) > length(unique(agemixing.df.IDs$parent)) & length(unique(agemixing.df.IDs$parent)) > 1 
  
  dirfasttree <- dirfasttree
  
  
  
  if(length(mAr.IDs)>10){
    
    
    choose.sequence.ind(pool.seq.file = paste0(sub.dir.rename,"/C.Epidemic.fas"),
                        select.vec = mAr.IDs,
                        name.file = paste0(sub.dir.rename,"/",paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta")))
    
    
    mAr.IDs.tree.calib <- phylogenetic.tree.fasttree.par(dir.tree = dirfasttree,
                                                         sub.dir.rename = sub.dir.rename,
                                                         fasttree.tool = "FastTree", # FastTreeMP
                                                         calendar.dates = "samplingtimes.all.csv",
                                                         simseqfile = paste0("cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta"),
                                                         count.start = 1977,
                                                         endsim = endpoint,
                                                         clust = FALSE) # TRUE
    
    
    
    tree.cal.cov.35.IDs <- read.tree(paste0(sub.dir.rename, paste0("/calibrated.tree.cov.",seq.cov, ".mAr.IDs.C.Epidemic.Fasta.tree")))
    
    
    tree.cal <- tree.cal.cov.35.IDs
    
    # Statistics for calibrated tree ------------------
    
    
    
    # Mean height of internal nodes
    
    nodeHeights_1 <- phytools::nodeHeights(tree.cal) # similar to node.depth.edgelength(tree) - Branch lengths:
    
    
    # It's clear from a casual inspection of the matrix that each parent node height (in the right column) 
    # is represented twice and only twice. Thus, if we exclude the root node (zero height), 
    # we can just take the mean of nodeHeights_1[,1].
    
    
    mean_nodeHeights_1 <- mean(sort(nodeHeights_1[,1])[3:nrow(nodeHeights_1)]) # important
    
    median_nodeHeights_1 <- median(sort(nodeHeights_1[,1])[3:nrow(nodeHeights_1)])
    
    sd_nodeHeights_1 <- sd(sort(nodeHeights_1[,1])[3:nrow(nodeHeights_1)])
    
    
    colless.phylo_1 <- phyloTop::colless.phylo(tree.cal, normalise = TRUE)
    
    sackin.phylo_1 <- phyloTop::sackin.phylo(tree.cal, normalise = TRUE)
    
    
    branch_leng_summary_1 <- summary(tree.cal$edge.length) # summary(tree.cal) # https://rdrr.io/cran/ape/man/summary.phylo.html
    
    
    Depths_1 <- phyloTop::getDepths(tree.cal) # depth of tips and nodes
    
    mean.tipsDepths_1 <- mean(Depths_1$tipDepths)
    
    mean.nodesDepths_1 <- mean(Depths_1$nodeDepths)
    
    
    median.tipsDepths_1 <- median(Depths_1$tipDepths)
    
    median.nodesDepths_1 <- median(Depths_1$nodeDepths)
    
    
    sd.tipsDepths_1 <- sd(Depths_1$tipDepths)
    
    sd.nodesDepths_1 <- sd(Depths_1$nodeDepths)
    
    
    maxHeight_1 <- phyloTop::maxHeight(tree.cal, normalise = TRUE)
    
    
    
    features.calib.tree.numeric <- as.numeric(c(mean_nodeHeights_1, median_nodeHeights_1, sd_nodeHeights_1,
                                                colless.phylo_1, sackin.phylo_1,
                                                branch_leng_summary_1,
                                                mean.tipsDepths_1, median.tipsDepths_1, sd.tipsDepths_1,
                                                mean.nodesDepths_1, median.nodesDepths_1, sd.nodesDepths_1,
                                                maxHeight_1))
    
    
    
    names.all <- c("mean_nodeHeights", "median_nodeHeights", "sd_nodeHeights",
                   "colless", "sackin",
                   "min_BL", "1Q_BL", "median_BL", "mean_BL", "3Q_BL", "max_BL",
                   "mean.tipsDepths", "median.tipsDepths", "sd.tipsDepths",
                   "mean.nodesDepths", "median.nodesDepths", "sd.nodesDepths",
                   "maxHeight")
    
    names(features.calib.tree.numeric) <- names.all
    
  }else{
    
    names.all <- c("mean_nodeHeights", "median_nodeHeights", "sd_nodeHeights",
                   "colless", "sackin",
                   "min_BL", "1Q_BL", "median_BL", "mean_BL", "3Q_BL", "max_BL",
                   "mean.tipsDepths", "median.tipsDepths", "sd.tipsDepths",
                   "mean.nodesDepths", "median.nodesDepths", "sd.nodesDepths",
                   "maxHeight")
    
    features.calib.tree.numeric <- rep(NA, length(names.all))
    
    names(features.calib.tree.numeric) <- names.all
    
  }
  
  return(features.calib.tree.numeric)
  
}

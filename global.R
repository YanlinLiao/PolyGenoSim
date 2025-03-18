###########################Call Environments########################
#install all useful package
if("shiny" %in% rownames(installed.packages())){
  library(shiny)
}else{
  install.packages("shiny")
  library(shiny)
}

if("shinydashboard" %in% rownames(installed.packages())){
  library(shinydashboard)
}else{
  install.packages("shinydashboard")
  library(shinydashboard)
}


if("plyr" %in% rownames(installed.packages())){
  library(plyr)
}else{
  install.packages("plyr")
  library(plyr)
}


# if("xlsx" %in% rownames(installed.packages())){
#   library(xlsx)
# }else{
#   install.packages("xlsx")
#   library(xlsx)
# }


if("MASS" %in% rownames(installed.packages())){
  library(MASS)
}else{
  install.packages("MASS")
  library(MASS)
}
###########################Check the file########################
if(!all(c("PedigreeSim.jar","lib") %in% dir())) stop("Lack depdent pacakges")


###########################Useful Functions########################
#readDatafile: from fitPolyTools
readDatfile <- function (file, header = TRUE, sep = "\t", check.names = FALSE) {
  return(read.table(file = file, header = header, sep = sep, 
                    check.names = check.names))
}

#call pedigreesim to run
run.PedigreeSim <- function(parfile, folder,path.to.PedigreeSim="PedigreeSim.jar") {
  ##To execute (e.g.), type: run.PedigreeSim("myparfile.par")
  path <- getwd()
  setwd(folder)
  ps <- system2(command = "java",
                args = c("-jar",
                         path.to.PedigreeSim,
                         parfile),
                stdout = TRUE,
                stderr = TRUE)
  setwd(path)
  
} #run.PedigreeSim()







## General function to create input files for PedSim for cross-ploidy populations

#' @param markermat Matrix of 3 columns containing columns P1 dosage, P2 dosage and N markers per LG
#' @param LG_number Number of linkage groups to simulate
#' @param filename Stem of the filename to use
#' @param ploidy Ploidy of parent 1
#' @param ploidy Ploidy of parent 2, assumed to be different from parent 1
#' @param prefP Preferential pairing parameter (assumed to be a single value across all LGs)
#' @param quads The rate of quadrivalent formation in parental meiosis
#' @param len Average length of the linkage groups in cM
#' @param centro Position of the centromeres, in cM
#' @param popnum The size of the F1 mapping population
#' @param mapfun The mapping function to use, either Haldane's or Kosambi's
#' @param marker.positions If specified, the position (in cM) of markers on the chromosomes, must equal number of markers per chromosome.
PedSim_input_complex <- function(markermat = matrix(c(1,0,10,
                                                      0,1,10,
                                                      1,1,10,
                                                      1,2,10,
                                                      2,0,5,
                                                      3,1,5,
                                                      4,1,5), ncol = 3,byrow = T),
                                 LG_number = 1,
                                 filename = "PedTry",
                                 folder = 'test',
                                 ploidy = 6,
                                 ploidy2 = 4,
                                 prefP_P1 = 0,
                                 prefP_P2 = 0.5,
                                 quads_P1 = 1,
                                 quads_P2 = 0,
                                 len = 100,
                                 centro = 50,
                                 popnum = 100,
                                 mapfun = 'haldane',
                                 marker.positions = NULL){
  dir.create(folder)
  file.copy(c("lib","PedigreeSim.jar"),  paste0(folder), recursive = TRUE)
  
  if(ncol(markermat) != 3 | !is.matrix(markermat)) stop("marker_type.mat should be a matrix of 3 columns")
  if(!is.null(marker.positions) & length(marker.positions) != sum(markermat[,3]))
    stop("marker.positions should be a vector of same length as the sum of column 3 of markermat!")
  mapfun <- match.arg(mapfun,choices = c("haldane","kosambi"))
  
  ## Test that the markermat is compatible with parental ploidies:
  if(any(markermat[,1] > ploidy) | any(markermat[,2] > ploidy2))
    stop(paste0("Supplied markermat is incompatible with parental ploidies!\nPossible P1 dosages: ",
                paste(seq(ploidy),collapse=","),"\nPossible P2 dosages: ",paste(seq(ploidy2),collapse=", ")))
  
  ## P1 and P2 meioses have to be simulated separately, then combined. Add dummy nulliplex parent to each:
  P1markermat <- P2markermat <- markermat
  P1markermat[,2] <- P2markermat[,1] <- 0
  
  nmark <- sum(markermat[,3])
  
  mapdat <- lapply(1:LG_number, function(lg){
    
    P1genmat <- matrix(0,ncol = 2*ploidy, nrow = nmark)
    P2genmat <- matrix(0,ncol = 2*ploidy2, nrow = nmark)
    
    ## randomise the positions of the markers
    mark.pos <- sample(seq(0,len,0.01),nmark)
    
    counter <- 1
    
    for(r in 1:nrow(markermat)){
      for(j in counter:(counter + markermat[r,3] - 1)){
        # message(paste(r,j))
        if(P1markermat[r,1] > 0) P1genmat[j,(1:ploidy)[sample((1:ploidy),P1markermat[r,1])]] <- 1
        if(P2markermat[r,2] > 0) P2genmat[j,((ploidy2 + 1):(2*ploidy2))[sample((1:ploidy2),P2markermat[r,2])]] <- 1
      }
      counter <- counter + markermat[r,3]
    }
    
    ## Generate marker names:
    mtypes <-unlist(apply(markermat,1,function(x) rep(paste0(x[1],"x",x[2]),x[3])))
    
    # phase <- paste0("_h",sapply(lapply(seq(nrow(genmat)),function(k) which(genmat[k,]!=0)),
    # function(ch) paste0(ch,collapse = "")))
    
    mnames <- paste0("LG",lg,"_",mark.pos,"_",mtypes)#,phase)
    
    P1gen.df <- as.data.frame(cbind("marker" = mnames,P1genmat))
    P2gen.df <- as.data.frame(cbind("marker" = mnames,P2genmat))
    
    colnames(P1gen.df)[2:ncol(P1gen.df)] <- c(paste0("P1_",1:ploidy),
                                              paste0("P2_",1:ploidy))
    
    colnames(P2gen.df)[2:ncol(P2gen.df)] <- c(paste0("P1_",1:ploidy2),
                                              paste0("P2_",1:ploidy2))
    
    map.df <- data.frame(marker = mnames, chromosome = lg, position = mark.pos)
    
    #Re-order:
    P1gen.df <- P1gen.df[order(mark.pos),]
    P2gen.df <- P2gen.df[order(mark.pos),]
    
    map.df <- map.df[order(mark.pos),]
    
    return(list("P1gen" = P1gen.df,
                "P2gen" = P2gen.df,
                "map" = map.df))
  })
  
  # Generate chromfile. For now have same settings per LG (easy to alter this afterwards if needed)
  # chrom.out <- data.frame(chromosome = 1:LG_number,
  #                         length = rep(len,LG_number),
  #                         centromere = rep(centro,LG_number),
  #                         prefPairing = rep(prefP,LG_number),
  #                         quadrivalents = rep(quads,LG_number))
  
  P1chrom.out <- data.frame(chromosome = 1:LG_number,
                            length = rep(len,LG_number),
                            centromere = rep(centro,LG_number),
                            prefPairing = rep(prefP_P1,LG_number),
                            quadrivalents = rep(quads_P1,LG_number))
  P2chrom.out <- data.frame(chromosome = 1:LG_number,
                            length = rep(len,LG_number),
                            centromere = rep(centro,LG_number),
                            prefPairing = rep(prefP_P2,LG_number),
                            quadrivalents = rep(quads_P2,LG_number))
  
  ## Extract and compile map and gen df:
  map.out <- do.call(rbind,lapply(mapdat, function(x) x$map))
  
  P1gen.out <- do.call(rbind,lapply(mapdat, function(x) x$P1gen))
  P2gen.out <- do.call(rbind,lapply(mapdat, function(x) x$P2gen))
  
  ## Write out the files to the working directory
  P1genpath <- paste0(folder,"/",filename,c("P1.gen"))
  P2genpath <- paste0(folder,"/",filename,c("P2.gen"))
  
  mappath <- paste0(folder,"/",filename,c(".map"))
  # chrompath <- paste0(folder,"/",filename,c(".chrom"))
  
  P1chrompath <- paste0(folder,"/",filename,c("P1.chrom"))
  P2chrompath <- paste0(folder,"/",filename,c("P2.chrom"))
  
  P1parpath <- paste0(folder,"/",filename,c("P1.par"))
  P2parpath <- paste0(folder,"/",filename,c("P2.par"))
  
  write.table(P1gen.out, P1genpath, quote = FALSE, sep = " ", row.names = FALSE)
  write.table(P2gen.out, P2genpath, quote = FALSE, sep = " ", row.names = FALSE)
  write.table(map.out, mappath, quote = FALSE, sep = " ", row.names = FALSE)
  # write.table(chrom.out, chrompath, quote = FALSE, sep = " ", row.names = FALSE)
  write.table(P1chrom.out, P1chrompath, quote = FALSE, sep = " ", row.names = FALSE)
  write.table(P2chrom.out, P2chrompath, quote = FALSE, sep = " ", row.names = FALSE)
  
  ## Finally write parfile
  for(parent in 1:2){
    fileConn<- file(c(P1parpath,P2parpath)[parent])
    writeLines(c(paste("PLOIDY =",c(ploidy,ploidy2)[parent]),
                 paste("MAPFUNCTION =",toupper(mapfun)),
                 "MISSING = NA",
                 paste0("CHROMFILE = ",filename,'P',parent,".chrom"),
                 "POPTYPE = F1",
                 paste0("POPSIZE = ", popnum),
                 paste0("MAPFILE = ", filename, ".map"),
                 paste0("FOUNDERFILE = ", filename, "P",parent,".gen"),
                 paste0("OUTPUT = ", filename,"P",parent,"_out"),
                 paste0("NATURALPAIRING = ",0)), fileConn)
    close(fileConn)
  }
  
  Sys.sleep(0.5) #make sure files are written out fully
  
  path <- getwd()
  setwd(folder)
  ps <- system2(command = "java",
                args = c("-jar",
                         "PedigreeSim.jar",
                         strsplit(P1parpath,"/")[[1]][2]),
                stdout = TRUE,
                stderr = TRUE)
  ps <- system2(command = "java",
                args = c("-jar",
                         "PedigreeSim.jar",
                         strsplit(P2parpath,"/")[[1]][2]),
                stdout = TRUE,
                stderr = TRUE)
  
  ## Load the allele dosages from both parents:
  P1dose <- as.matrix(read.table(paste0(filename, "P1_out_alleledose.dat"),
                                 header = TRUE, stringsAsFactors = FALSE,
                                 sep = "\t",row.names=1))
  P2dose <- as.matrix(read.table(paste0(filename, "P2_out_alleledose.dat"),
                                 header = TRUE, stringsAsFactors = FALSE,
                                 sep = "\t",row.names=1))
  
  ## Any 0 x 0 markers were converted by PedSim to ploidy x ploidy. Convert back to 0 x 0 first:
  P1_00 <- which(P1dose[,1]==ploidy & P1dose[,2] == ploidy)
  P2_00 <- which(P2dose[,1]==ploidy2 & P2dose[,2] == ploidy2)
  
  if(length(P1_00) > 0) P1dose[P1_00,] <- 0
  if(length(P2_00) > 0) P2dose[P2_00,] <- 0
  
  ## Add the dosages together:
  F1dose <- P1dose + P2dose
  
  write.table(cbind("marker" = rownames(F1dose),F1dose), paste0(filename,"_out_alleledose.dat"), sep = "\t", row.names=FALSE)
  
  ## Generate founderalleles file:
  P1founder <- as.matrix(read.table(paste0(filename, "P1_out_founderalleles.dat"),
                                    header = TRUE, stringsAsFactors = FALSE,
                                    sep = "\t",row.names=1))
  P2founder <- as.matrix(read.table(paste0(filename, "P2_out_founderalleles.dat"),
                                    header = TRUE, stringsAsFactors = FALSE,
                                    sep = "\t",row.names=1))
  
  P1founder.X <- P1founder[,-seq(2*ploidy)]
  P2founder.X <- P2founder[,-seq(2*ploidy2)]
  
  ## Generate new colnames for P2founder.1
  popNames <- colnames(F1dose)[3:ncol(F1dose)]
  
  P1founder.1 <- P1founder.X[,sort(c(matrix(1:ncol(P1founder.X),ncol=ploidy,byrow = TRUE)[,1:(ploidy/2)]))]
  P2founder.1 <- P2founder.X[,sort(c(matrix(1:ncol(P2founder.X),ncol=ploidy2,byrow = TRUE)[,(ploidy2/2 + 1):(ploidy2)]))]
  
  ## P1founder.1 numbering is ok - starts from 0 as expected. P2founder.1 numbering is not right and must be adjusted.
  d <- ploidy - ploidy2
  P2founder.1 <- P2founder.1 + d
  
  colnames(P2founder.1) <- sort(c(sapply((ploidy/2 + 1):((ploidy + ploidy2)/2),
                                         function(n) paste(popNames,n,sep="_"))))
  
  founder.out <- cbind(P1founder.1,P2founder.1)
  founder.out <- founder.out[,order(colnames(founder.out))]
  
  ## Add parental columns back:
  P2parcols <- do.call(rbind,lapply(1:nrow(P2founder), function(cl) ploidy:(ploidy + ploidy2 - 1)))
  colnames(P2parcols) <- paste0("P2_",1:ploidy2)
  founder.out <- cbind("marker" = rownames(P1founder),P1founder[,1:ploidy],P2parcols,founder.out)
  
  write.table(founder.out, paste0(filename,"_out_founderalleles.dat"), sep = "\t", row.names=FALSE)
  
  ## Generate the .gen file for the cross-ploidy population:
  F1gen.out <- cbind(P1gen.out[,1:(ploidy+1)],P2gen.out[,(ploidy2 + 2):(2*ploidy2 + 1)])
  write.table(F1gen.out, paste0(filename,".gen"), quote = FALSE, sep = " ", row.names = FALSE)
  
  ## Read in the .hsa and .hsb files:
  P1hsa <- read.table(paste0(filename, "P1_out.hsa"),header=FALSE,stringsAsFactors = FALSE)
  P1hsb <- read.table(paste0(filename, "P1_out.hsb"),header=FALSE,stringsAsFactors = FALSE)
  P2hsa <- read.table(paste0(filename, "P2_out.hsa"),header=FALSE,stringsAsFactors = FALSE)
  P2hsb <- read.table(paste0(filename, "P2_out.hsb"),header=FALSE,stringsAsFactors = FALSE)
  #if they are not in the same length, make them the same
  asjust_length <- function(tmp1, tmp2){
    combine_list <- list(tmp1,tmp2)
    colLen <- unlist(lapply(combine_list,ncol))
    if(diff(colLen) != 0){
      chosen <- which.min(colLen);NTchosen <- which.max(colLen)
      NA_matrix <- t(matrix(rep(NA,nrow(combine_list[[chosen]]) * abs(diff(colLen))),abs(diff(colLen))))
      chosen_tmp <- cbind(combine_list[[chosen]],NA_matrix)
      colnames(chosen_tmp) <- colnames(combine_list[[NTchosen]])
      combine_list[[chosen]] <- chosen_tmp
    }
    return(list(combine_list[[1]],combine_list[[2]]))
  } #add NA to the one with less columns
  hsa_list <- asjust_length(tmp1 = P1hsa,
                            tmp2 = P2hsa)
  P1hsa <- hsa_list[[1]]; P2hsa <- hsa_list[[2]]
  hsb_list <- asjust_length(tmp1 = P1hsb,
                            tmp2 = P2hsb)
  P1hsb <- hsb_list[[1]]; P2hsb <- hsb_list[[2]]
  
  
  ## Get rid of P2 from the P1 files, P1 from the P2 files, and combine...
  P1hsa.hd <- P1hsa[1:(ploidy*LG_number),]; P2hsa.hd <- P2hsa[(ploidy2*LG_number+1):(2*ploidy2*LG_number),]
  P1hsb.hd <- P1hsb[1:(ploidy*LG_number),]; P2hsb.hd <- P2hsb[(ploidy2*LG_number+1):(2*ploidy2*LG_number),]
  
  P1hsa <- P1hsa[-(1:(2*ploidy*LG_number)),]; P2hsa <- P2hsa[-(1:(2*ploidy2*LG_number)),]
  P1hsb <- P1hsb[-(1:(2*ploidy*LG_number)),]; P2hsb <- P2hsb[-(1:(2*ploidy2*LG_number)),]
  
  ## Add difference in ploidy (d) to P2hsa
  P2hsa.hd[,ploidy] <- P2hsa.hd[,ploidy] + d
  P2hsa[,(ploidy + 1):ncol(P2hsa)] <- P2hsa[,(ploidy + 1):ncol(P2hsa)] + d
  
  P1set <- which(rep(c(rep(T,ploidy/2),rep(F,ploidy/2)),popnum*LG_number))
  P2set <- which(rep(c(rep(F,ploidy2/2),rep(T,ploidy2/2)),popnum*LG_number))
  
  P12hsa <- rbind(P1hsa[P1set,],P2hsa[P2set,])
  P12hsb <- rbind(P1hsb[P1set,],P2hsb[P2set,])
  reord <- order(P12hsa[,1],P12hsa[,2])
  
  P12hsa <- rbind(P1hsa.hd,P2hsa.hd,P12hsa[reord,])
  P12hsb <- cbind("",rbind(P1hsb.hd,P2hsb.hd,P12hsb[reord,])) #Space in first column is in original .hsb files
  
  write.table(P12hsa, paste0(filename,"_out.hsa"), sep = "\t", quote=FALSE, row.names=FALSE, col.names = FALSE)
  write.table(P12hsb, paste0(filename,"_out.hsb"), sep = "\t", quote = FALSE, row.names=FALSE, col.names = FALSE)
  
  ## Write the .ped file (can use either of the parental .ped files)
  pedfl <- read.table(paste0(filename,"P1_out.ped"), header = TRUE, stringsAsFactors = FALSE)
  write.table(pedfl, paste0(filename,"_out.ped"), quote = FALSE, row.names = FALSE, sep = "\t")
  
  
  ## Generate the genotypes file
  P1genotypes <- as.matrix(read.table(paste0(filename, "P1_out_genotypes.dat"),
                                      header = TRUE, stringsAsFactors = FALSE,
                                      sep = "\t",row.names=1))
  P2genotypes <- as.matrix(read.table(paste0(filename, "P2_out_genotypes.dat"),
                                      header = TRUE, stringsAsFactors = FALSE,
                                      sep = "\t",row.names=1))
  
  P1genotypes.parental <- P1genotypes[,1:ploidy]; P2genotypes.parental <- P2genotypes[,(ploidy2+1):(2*ploidy2)]
  P1genotypes <- P1genotypes[,-seq(2*ploidy)]; P2genotypes <- P2genotypes[,-seq(2*ploidy2)]
  
  if(any(rownames(P1genotypes) != rownames(P2genotypes))) stop("Error in combining genotype.dat files!")
  
  P1set <- which(rep(c(rep(T,ploidy/2),rep(F,ploidy/2)),popnum))
  P2set <- which(rep(c(rep(F,ploidy2/2),rep(T,ploidy2/2)),popnum))
  
  colnames(P2genotypes)[P2set] <- paste0(substr(colnames(P2genotypes[,P2set]),1,nchar(colnames(P2genotypes[,P2set]))-1),(ploidy/2+1):(ploidy/2+ploidy2/2))
  
  P12genotypes <- cbind(P1genotypes[,P1set],P2genotypes[,P2set])
  P12genotypes <- cbind(P1genotypes.parental,P2genotypes.parental,P12genotypes[,order(colnames(P12genotypes))])
  
  write.table(cbind("marker" = rownames(P12genotypes),P12genotypes),
              paste0(filename,"_out_genotypes.dat"), sep = "\t", row.names=FALSE)
  
  ## Tidy-up - delete the intermediate files:
  file.remove(paste0(filename,"P1.gen"));file.remove(paste0(filename,"P2.gen"))
  file.remove(paste0(filename,"P1.par"));file.remove(paste0(filename,"P2.par"))
  file.remove(paste0(filename,"P1_out.hsa"));file.remove(paste0(filename,"P2_out.hsa"))
  file.remove(paste0(filename,"P1_out.hsb"));file.remove(paste0(filename,"P2_out.hsb"))
  file.remove(paste0(filename,"P1_out.ped"));file.remove(paste0(filename,"P2_out.ped"))
  file.remove(paste0(filename,"P1_out_founderalleles.dat"));file.remove(paste0(filename,"P2_out_founderalleles.dat"))
  file.remove(paste0(filename,"P1_out_alleledose.dat"));file.remove(paste0(filename,"P2_out_alleledose.dat"))
  file.remove(paste0(filename,"P1_out_genotypes.dat"));file.remove(paste0(filename,"P2_out_genotypes.dat"))
  
  
  print("Thank you. Please check your R working directory.")
  setwd(path)
  
} #PedSim_input_complex


## General function to create input files for PedSim for any ploidy level (assumed parental ploidies are equal)

#' @param markermat Matrix of 3 columns containing columns P1 dosage, P2 dosage and N markers per LG
#'
#' @param LG_number Number of linkage groups to simulate
#' @param filename Stem of the filename to use
#' @param ploidy Ploidy of the species, assumed to be even and equal between parents
#' @param prefP Preferential pairing parameter (assumed to be a single value across all LGs)
#' @param quads The rate of quadrivalent formation in parental meiosis
#' @param len Average length of the linkage groups in cM
#' @param centro Position of the centromeres, in cM
#' @param popnum The size of the F1 mapping population
#' @param mapfun The mapping function to use, either Haldane's or Kosambi's
#' @param marker.positions If specified, the position (in cM) of markers on the chromosomes, must equal number of markers per chromosome.
#' @param folder If specified, the path to the folder to which files should be written. by default \code{NULL}
#' @param genfilename Option to specify .gen file
#' @param mapfilename Option to specify .map file
#' @param chromfilename Option to specify .chrom file
#' @param poptype The population type that you would like to use, can be 'F1', 'F2', 'BC'
#' @param runPedigreeSim Should PedigreeSim be run directly (using \code{system} call) - by default \code{TRUE}
PedSim_input <- function(markermat = matrix(c(1,0,100,
                                              0,1,100,
                                              1,1,75,
                                              1,2,50,
                                              2,0,40), ncol = 3,byrow = T),
                         LG_number = 1,
                         filename = "PedSim_input",
                         ploidy = 4,
                         prefP = 0,
                         quads = 0,
                         len = 100,
                         centro = 50,
                         popnum = 200,
                         mapfun = c("haldane","kosambi"),
                         marker.positions = NULL,
                         folder = NULL,
                         poptype = c('F1','F2','BC'),
                         genfilename = NULL,
                         mapfilename = NULL,
                         chromfilename = NULL,
                         runPedigreeSim = TRUE){
  dir.create(folder)
  file.copy(c("lib","PedigreeSim.jar"),  paste0(folder), recursive = TRUE)
  if(ncol(markermat) != 3 | !is.matrix(markermat)) stop("marker_type.mat should be a matrix of 3 columns")
  if(!is.null(marker.positions) & length(marker.positions) != sum(markermat[,3]))
    stop("marker.positions should be a vector of same length as the sum of column 3 of markermat!")
  mapfun <- match.arg(mapfun)
  
  mapdat <- lapply(1:LG_number, function(lg){
    nmark <- sum(markermat[,3])
    genmat <- matrix(0,ncol = 2*ploidy, nrow = nmark)
    
    ## randomise the positions of the markers
    mark.pos <- sample(seq(0,len,0.01),nmark)
    
    counter <- 1
    
    for(r in 1:nrow(markermat)){
      for(j in counter:(counter + markermat[r,3] - 1)){
        # message(paste(r,j))
        if(markermat[r,1] > 0) genmat[j,(1:ploidy)[sample((1:ploidy),markermat[r,1])]] <- 1
        if(markermat[r,2] > 0) genmat[j,((ploidy + 1):(2*ploidy))[sample((1:ploidy),markermat[r,2])]] <- 1
      }
      counter <- counter + markermat[r,3]
    }
    
    ## Generate marker names:
    mtypes <-unlist(apply(markermat,1,function(x) rep(paste0(x[1],"x",x[2]),x[3])))
    
    phase <- paste0("_h",sapply(lapply(seq(nrow(genmat)),function(k) which(genmat[k,]!=0)),
                                function(ch) paste0(ch,collapse = "")))
    
    mnames <- paste0("LG",lg,"_",mark.pos,"_",mtypes,phase)
    
    gen.df <- as.data.frame(cbind("marker" = mnames,genmat))
    colnames(gen.df)[2:ncol(gen.df)] <- c(paste0("P1_",1:ploidy),
                                          paste0("P2_",1:ploidy))
    
    map.df <- data.frame(marker = mnames, chromosome = lg, position = mark.pos)
    
    #Re-order:
    gen.df <- gen.df[order(mark.pos),]
    map.df <- map.df[order(mark.pos),]
    
    return(list("gen" = gen.df,
                "map" = map.df))
  })
  
  # Generate chromfile. For now have same settings per LG (easy to alter this afterwards if needed)
  chrom.out <- data.frame(chromosome = 1:LG_number,
                          length = rep(len,LG_number),
                          centromere = rep(centro,LG_number),
                          prefPairing = rep(prefP,LG_number),
                          quadrivalents = rep(quads,LG_number))
  
  ## Extract and compile map and gen df:
  map.out <- do.call(rbind,lapply(mapdat, function(x) x$map))
  gen.out <- do.call(rbind,lapply(mapdat, function(x) x$gen))
  
  ## Write out the files to the working directory
  genpath <- ifelse(is.null(folder),paste0(filename,c(".gen")),file.path(folder,paste0(filename,c(".gen"))))
  mappath <- ifelse(is.null(folder),paste0(filename,c(".map")),file.path(folder,paste0(filename,c(".map"))))
  chrompath <- ifelse(is.null(folder),paste0(filename,c(".chrom")),file.path(folder,paste0(filename,c(".chrom"))))
  parpath <- ifelse(is.null(folder),paste0(filename,c(".par")),file.path(folder,paste0(filename,c(".par"))))
  
  if(is.null(genfilename)) write.table(gen.out, genpath, quote = FALSE, sep = " ", row.names = FALSE)
  if(is.null(mapfilename)) write.table(map.out, mappath, quote = FALSE, sep = " ", row.names = FALSE)
  if(is.null(chromfilename)) write.table(chrom.out, chrompath, quote = FALSE, sep = " ", row.names = FALSE)
  
  if(is.null(genfilename)){
    genfilename <- paste0(filename, ".gen")
  }
  
  if(is.null(mapfilename)){
    mapfilename <- paste0(filename, ".map")
  }
  
  if(is.null(chromfilename)){
    chromfilename <- paste0(filename, ".chrom")
  }
  
  ## Finally write parfile
  fileConn<-file(parpath)
  writeLines(c(paste("PLOIDY =",ploidy),
               paste("MAPFUNCTION =",toupper(mapfun)),
               "MISSING = NA",
               paste0("CHROMFILE = ",chromfilename),
               paste0('POPTYPE = ',poptype),
               paste0("POPSIZE = ", popnum),
               paste0("MAPFILE = ", mapfilename),
               paste0("FOUNDERFILE = ", genfilename),
               paste0("OUTPUT = ", filename, "_out"),
               "NATURALPAIRING = 0"), fileConn)
  close(fileConn)
  
  Sys.sleep(0.5) #make sure files are written out fully
  path <- getwd()
  setwd(folder)
  ps <- system2(command = "java",
                args = c("-jar",
                         "PedigreeSim.jar",
                         strsplit(parpath,"/")[[1]][2]),
                stdout = TRUE,
                stderr = TRUE)
  setwd(path)
  print("Thank you. Please check your R working directory.")
  
} #PedSim_input


## General function to create input files for PedSim for selfing population 

#' @param markermat Matrix of 2 columns containing columns parental dosage and N markers per LG
#'
#' @param LG_number Number of linkage groups to simulate
#' @param filename Stem of the filename to use
#' @param ploidy Ploidy of the species, assumed to be even and equal between parents
#' @param prefP Preferential pairing parameter (assumed to be a single value across all LGs)
#' @param quads The rate of quadrivalent formation in parental meiosis
#' @param len Average length of the linkage groups in cM
#' @param centro Position of the centromeres, in cM
#' @param popnum The size of the F1 mapping population
#' @param mapfun The mapping function to use, either Haldane's or Kosambi's
#' @param marker.positions If specified, the position (in cM) of markers on the chromosomes, must equal number of markers per chromosome.
#' @param folder If specified, the path to the folder to which files should be written. by default \code{NULL}
#' @param genfilename Option to specify .gen file
#' @param mapfilename Option to specify .map file
#' @param chromfilename Option to specify .chrom file
#' @param runPedigreeSim Should PedigreeSim be run directly (using \code{system} call) - by default \code{TRUE}
PedSim_Selfing <- function(markermat = matrix(c(0,100,
                                                1,100,
                                                2,75,
                                                3,50,
                                                4,40), ncol = 2,byrow = T),
                           LG_number = 1,
                           filename = "PedSim_input",
                           ploidy = 4,
                           prefP = 0,
                           quads = 0,
                           len = 100,
                           centro = 50,
                           popnum = 200,
                           mapfun = c("haldane","kosambi"),
                           marker.positions = NULL,
                           folder = NULL,
                           genfilename = NULL,
                           mapfilename = NULL,
                           chromfilename = NULL,
                           runPedigreeSim = TRUE){
  dir.create(folder)
  file.copy(c("lib","PedigreeSim.jar"),  paste0(folder), recursive = TRUE)
  if(ncol(markermat) != 2 | !is.matrix(markermat)) stop("marker_type.mat should be a matrix of 2 columns")
  if(!is.null(marker.positions) & length(marker.positions) != sum(markermat[,2]))
    stop("marker.positions should be a vector of same length as the sum of column 2 of markermat!")
  mapfun <- match.arg(mapfun)
  
  mapdat <- lapply(1:LG_number, function(lg){
    nmark <- sum(markermat[,2])
    genmat <- matrix(0,ncol = ploidy, nrow = nmark)
    
    ## randomise the positions of the markers
    mark.pos <- sample(seq(0,len,0.01),nmark)
    
    counter <- 1
    
    for(r in 1:nrow(markermat)){
      for(j in counter:(counter + markermat[r,2] - 1)){
        # message(paste(r,j))
        if(markermat[r,1] > 0) genmat[j,(1:ploidy)[sample((1:ploidy),markermat[r,1])]] <- 1
      }
      counter <- counter + markermat[r,2]
    }
    
    ## Generate marker names:
    mtypes <-unlist(apply(markermat,1,function(x) rep(x[1],x[2])))
    mnames <- paste0("LG",lg,"_",mark.pos,"_",mtypes)
    
    gen.df <- as.data.frame(cbind("marker" = mnames,genmat))
    colnames(gen.df)[2:ncol(gen.df)] <- c(paste0("P1_",1:ploidy))
    
    map.df <- data.frame(marker = mnames, chromosome = lg, position = mark.pos)
    
    #Re-order:
    gen.df <- gen.df[order(mark.pos),]
    map.df <- map.df[order(mark.pos),]
    
    return(list("gen" = gen.df,
                "map" = map.df))
  })
  
  # Generate chromfile. For now have same settings per LG (easy to alter this afterwards if needed)
  chrom.out <- data.frame(chromosome = 1:LG_number,
                          length = rep(len,LG_number),
                          centromere = rep(centro,LG_number),
                          prefPairing = rep(prefP,LG_number),
                          quadrivalents = rep(quads,LG_number))
  
  ## Extract and compile map and gen df:
  map.out <- do.call(rbind,lapply(mapdat, function(x) x$map))
  gen.out <- do.call(rbind,lapply(mapdat, function(x) x$gen))
  
  ## Write out the files to the working directory
  genpath <- ifelse(is.null(folder),paste0(filename,c(".gen")),file.path(folder,paste0(filename,c(".gen"))))
  mappath <- ifelse(is.null(folder),paste0(filename,c(".map")),file.path(folder,paste0(filename,c(".map"))))
  chrompath <- ifelse(is.null(folder),paste0(filename,c(".chrom")),file.path(folder,paste0(filename,c(".chrom"))))
  parpath <- ifelse(is.null(folder),paste0(filename,c(".par")),file.path(folder,paste0(filename,c(".par"))))
  
  if(is.null(genfilename)) write.table(gen.out, genpath, quote = FALSE, sep = " ", row.names = FALSE)
  if(is.null(mapfilename)) write.table(map.out, mappath, quote = FALSE, sep = " ", row.names = FALSE)
  if(is.null(chromfilename)) write.table(chrom.out, chrompath, quote = FALSE, sep = " ", row.names = FALSE)
  
  if(is.null(genfilename)){
    genfilename <- paste0(filename, ".gen")
  }
  
  if(is.null(mapfilename)){
    mapfilename <- paste0(filename, ".map")
  }
  
  if(is.null(chromfilename)){
    chromfilename <- paste0(filename, ".chrom")
  }
  
  ## Finally write parfile
  fileConn<-file(parpath)
  writeLines(c(paste("PLOIDY =",ploidy),
               paste("MAPFUNCTION =",toupper(mapfun)),
               "MISSING = NA",
               paste0("CHROMFILE = ",chromfilename),
               'POPTYPE = S1',
               paste0("POPSIZE = ", popnum),
               paste0("MAPFILE = ", mapfilename),
               paste0("FOUNDERFILE = ", genfilename),
               paste0("OUTPUT = ", filename, "_out"),
               "NATURALPAIRING = 0"), fileConn)
  close(fileConn)
  
  Sys.sleep(0.5) #make sure files are written out fully
  path <- getwd()
  setwd(folder)
  ps <- system2(command = "java",
                args = c("-jar",
                         "PedigreeSim.jar",
                         strsplit(parpath,"/")[[1]][2]),
                stdout = TRUE,
                stderr = TRUE)
  setwd(path)
  print("Thank you. Please check your R working directory.")
  
} #PedSim_Selfing


## General function to create input files for PedSim for any ploidy level (assumed parental ploidies are equal)

#' @param markermat Matrix of 3 columns containing columns P1 dosage, P2 dosage and N markers per LG
#'
#' @param LG_number Number of linkage groups to simulate
#' @param filename Stem of the filename to use
#' @param ploidy Ploidy of the species, assumed to be even and equal between parents
#' @param prefP Preferential pairing parameter (assumed to be a single value across all LGs)
#' @param quads The rate of quadrivalent formation in parental meiosis
#' @param len Average length of the linkage groups in cM
#' @param centro Position of the centromeres, in cM
#' @param popnum The size of the F1 mapping population
#' @param mapfun The mapping function to use, either Haldane's or Kosambi's
#' @param marker.positions If specified, the position (in cM) of markers on the chromosomes, must equal number of markers per chromosome.
#' @param folder If specified, the path to the folder to which files should be written. by default \code{NULL}
#' @param genfilename Option to specify .gen file
#' @param mapfilename Option to specify .map file
#' @param chromfilename Option to specify .chrom file
#' @param runPedigreeSim Should PedigreeSim be run directly (using \code{system} call) - by default \code{TRUE}
PedSim_pedprovided <- function(markermat = matrix(c(1,0,100,
                                                    0,1,100,
                                                    1,1,75,
                                                    1,2,50,
                                                    2,0,40), ncol = 3,byrow = T),
                               LG_number = 1,
                               filename = "PedSim_test",
                               ploidy = 4,
                               prefP = 0,
                               quads = 0,
                               len = 100,
                               centro = 50,
                               mapfun = c("haldane","kosambi"),
                               marker.positions = NULL,
                               folder = NULL,
                               genfilename = NULL,
                               mapfilename = NULL,
                               chromfilename = NULL,
                               runPedigreeSim = TRUE){
  dir.create(folder)
  file.copy(c("lib","PedigreeSim.jar"),  paste0(folder), recursive = TRUE)
  if(ncol(markermat) != 3 | !is.matrix(markermat)) stop("marker_type.mat should be a matrix of 3 columns")
  if(!is.null(marker.positions) & length(marker.positions) != sum(markermat[,3]))
    stop("marker.positions should be a vector of same length as the sum of column 3 of markermat!")
  mapfun <- match.arg(mapfun)
  
  mapdat <- lapply(1:LG_number, function(lg){
    nmark <- sum(markermat[,3])
    genmat <- matrix(0,ncol = 2*ploidy, nrow = nmark)
    
    ## randomise the positions of the markers
    mark.pos <- sample(seq(0,len,0.01),nmark)
    
    counter <- 1
    
    for(r in 1:nrow(markermat)){
      for(j in counter:(counter + markermat[r,3] - 1)){
        # message(paste(r,j))
        if(markermat[r,1] > 0) genmat[j,(1:ploidy)[sample((1:ploidy),markermat[r,1])]] <- 1
        if(markermat[r,2] > 0) genmat[j,((ploidy + 1):(2*ploidy))[sample((1:ploidy),markermat[r,2])]] <- 1
      }
      counter <- counter + markermat[r,3]
    }
    
    ## Generate marker names:
    mtypes <-unlist(apply(markermat,1,function(x) rep(paste0(x[1],"x",x[2]),x[3])))
    
    phase <- paste0("_h",sapply(lapply(seq(nrow(genmat)),function(k) which(genmat[k,]!=0)),
                                function(ch) paste0(ch,collapse = "")))
    
    mnames <- paste0("LG",lg,"_",mark.pos,"_",mtypes,phase)
    
    gen.df <- as.data.frame(cbind("marker" = mnames,genmat))
    colnames(gen.df)[2:ncol(gen.df)] <- c(paste0("P1_",1:ploidy),
                                          paste0("P2_",1:ploidy))
    
    map.df <- data.frame(marker = mnames, chromosome = lg, position = mark.pos)
    
    #Re-order:
    gen.df <- gen.df[order(mark.pos),]
    map.df <- map.df[order(mark.pos),]
    
    return(list("gen" = gen.df,
                "map" = map.df))
  })
  
  # Generate chromfile. For now have same settings per LG (easy to alter this afterwards if needed)
  chrom.out <- data.frame(chromosome = 1:LG_number,
                          length = rep(len,LG_number),
                          centromere = rep(centro,LG_number),
                          prefPairing = rep(prefP,LG_number),
                          quadrivalents = rep(quads,LG_number))
  
  ## Extract and compile map and gen df:
  map.out <- do.call(rbind,lapply(mapdat, function(x) x$map))
  gen.out <- do.call(rbind,lapply(mapdat, function(x) x$gen))
  
  ## Write out the files to the working directory
  genpath <- ifelse(is.null(folder),paste0(filename,c(".gen")),file.path(folder,paste0(filename,c(".gen"))))
  mappath <- ifelse(is.null(folder),paste0(filename,c(".map")),file.path(folder,paste0(filename,c(".map"))))
  chrompath <- ifelse(is.null(folder),paste0(filename,c(".chrom")),file.path(folder,paste0(filename,c(".chrom"))))
  parpath <- ifelse(is.null(folder),paste0(filename,c(".par")),file.path(folder,paste0(filename,c(".par"))))
  
  if(is.null(genfilename)) write.table(gen.out, genpath, quote = FALSE, sep = " ", row.names = FALSE)
  if(is.null(mapfilename)) write.table(map.out, mappath, quote = FALSE, sep = " ", row.names = FALSE)
  if(is.null(chromfilename)) write.table(chrom.out, chrompath, quote = FALSE, sep = " ", row.names = FALSE)
  
  if(is.null(genfilename)){
    genfilename <- paste0(filename, ".gen")
  }
  
  if(is.null(mapfilename)){
    mapfilename <- paste0(filename, ".map")
  }
  
  if(is.null(chromfilename)){
    chromfilename <- paste0(filename, ".chrom")
  }
  
  ## Finally write parfile
  fileConn<-file(parpath)
  writeLines(c(paste("PLOIDY =",ploidy),
               paste("MAPFUNCTION =",toupper(mapfun)),
               "MISSING = NA",
               paste0("CHROMFILE = ",chromfilename),
               paste0('PEDFILE = ',filename,'.ped'),
               paste0("MAPFILE = ", mapfilename),
               paste0("FOUNDERFILE = ", genfilename),
               paste0("OUTPUT = ", filename, "_out"),
               "NATURALPAIRING = 0"), fileConn)
  close(fileConn)
  
  Sys.sleep(0.5) #make sure files are written out fully
  path <- getwd()
  setwd(folder)
  ps <- system2(command = "java",
                args = c("-jar",
                         "PedigreeSim.jar",
                         strsplit(parpath,"/")[[1]][2]),
                stdout = TRUE,
                stderr = TRUE)
  setwd(path)
  print("Thank you. Please check your R working directory.")
  
} #PedSim_pedprovided

#generate N random number sum up to M
rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}

#read .dat file
readDatfile <- function (file, header = TRUE, sep = "\t", check.names = FALSE) {
  return(read.table(file = file, header = header, sep = sep, 
                    check.names = check.names))
}

#generate random code
random_code <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

####SeqSim
# SeqSim <- function(dsg= dosages,
#                    ploidy = 4,
#                    ploidy2 = NULL,
#                    ideal_depth = 60,  #ideal read depth
#                    mrk_effect = 0.43,
#                    individual_effect = 0.96,
#                    dispersion = 1.8, #for binomial model
#                    seq_e = 0.01,
#                    sys_b, #affect A/a unbalancement for hybridization (yunb)|affect XY unbalancement for hybridization (yunb)
#                    od,
#                    seed_number = 1)# overdispersion parameter used in rbetabinomial distribution
# {
#   #############If there is no previous input, give user chances to make their choice
#   if(is.null(seed_number)){
#     seed_number = round(runif(1,0,1000))
#   }
#   
#   #Read Depth Simulation
#   set.seed(seed_number)
#   V <- rnorm(nrow(dsg),0,mrk_effect) # marker level effect
#   set.seed(seed_number)
#   W <- rnorm(ncol(dsg),0, individual_effect)  # individual effect
#   
#   s_V <- quantile(V,probs = 0.95)
#   s_W <- quantile(W,probs = 0.95)
#   bt <- log(ideal_depth)      # overall mean
#   
#   if(bt < s_V |bt < s_W){
#     stop("Please choose a higher ideal_depth")
#   }
#   
#   l <- sapply(V, function(x){
#     exp(bt + x + W)
#   }) 
#   
#   depvar <- t(sapply(1:length(V), function(i){
#     rnbinom(n = nrow(l), # per marker
#             mu = l[,i], #get the mean for each individual
#             size = dispersion) # dispersion in overall
#   }))
#   
#   
#   ############  Modeling read counts ############
#   #### step 1 - convert dosages into read counts ratios: genotype/ploidy           
#   if(is.null(ploidy2)){
#     ploidy_temp <- replicate(ncol(dsg), ploidy)
#   }else{
#     ploidy_temp <- c(ploidy, ploidy2, replicate(ncol(dsg) - 2, (ploidy + ploidy2)/2))
#   }
#   thr <- sapply(1:ncol(dsg), function(x) {
#     round(dsg[,x]/ploidy_temp[x],3)
#   })
#   
#   #### step 2 - simulate the sequencing error amomng all reads                                            
#   # seqE is a number between 0.5 – 1%     
#   thr1 <- thr * (1 - seq_e) + (1 - thr) * seq_e
#   
#   #### step 3 - simulate the systematic bias from read mapping step                                     
#   ##In Updog, they treat the systematic bias (observed/non-observed) as a fixed term for all markers.
#   ##In SNParraySim, the effect could vary between markers, it is more realistic because observation is independent per marker for all individuals
#   set.seed(seed_number)
#   yunb <- rbeta(nrow(dsg),sys_b[1],sys_b[2])  #the number of markers because vary in marker level, it is in fact d/c+d
#   # yunb <- replicate(nrow(dsg),sys_b)
#   yunbset <- replicate(ncol(dsg),yunb)
#   # thr2 <- matrix(thr1*yunbset,nrow = nrow(dsg)) 
#   thr2 <- thr1/(yunbset* (1 - thr1) + thr1)  #same equation as updog
#   
#   #### step 4-  overdispersion
#   set.seed(seed_number)
#   # overdis <- replicate(ncol(dsg),od)
#   overdis <- rbeta(ncol(dsg),od[1],od[2])
#   overdisset <- t(replicate(nrow(dsg),overdis))
#   para <- matrix(adjust(thr2,overdisset),nrow = nrow(dsg))
#   
#   #### step 5 - simulate read counts from rbeta distribution                                             Updog
#   counts <- matrix(rbinom(n = length(depvar), size = depvar, para),nrow= nrow(dsg))    #binomial model [Yanlin], do not takes account in the overal error because error has been modeled more specifically in each different steps
#   
#   #define colnames & rownames
#   colnames(counts)  <- individuals <- colnames(depvar) <- colnames(dsg)
#   rownames(counts) <- markers <- rownames(depvar) <- rownames(dsg)
#   
#   ##### step 6 - prepare genotype calling template
#   Template <- do.call(rbind, lapply(markers, function(m) data.frame("MarkerName" = m,
#                                                                     "SampleName"  = names(depvar[m,]),
#                                                                     "ReadDepth" = as.numeric(depvar[m,individuals]),
#                                                                     "A" = as.numeric(counts[m,individuals]),
#                                                                     "a" = as.numeric(depvar[m,individuals] - counts[m,individuals]))))
#   Template$ratio <- Template$A/Template$ReadDepth
#   
#   ##### step 7 - prepare output
#   output <- list()
#   output[["ReadDepth"]] <- depvar
#   output[["ReadCounts"]] <- counts
#   output[["GenotypeCallTemplate"]] <- Template
#   Parameter <- list()
#   Parameter[["SequencingError"]] <- seq_e
#   Parameter[["SystematicBias"]] <- yunbset
#   Parameter[["Overdispersion"]] <- overdisset
#   output[["Parameter"]] <- Parameter
#   
#   return(output)
# }

#GenoSim
GenoSim <- function(dsg,
                    choice = "B", # A.SNP array, B.sequence reads
                    ploidy = 4, 
                    ploidy2 = 4,
                    seed_number = 1,
                    od = c(1,80), #overdispersion from beta distribution [0.1-1]
                    ale_b = c(100,1), #allelic bias (ale_b = Y/X) [0.2-0.8]
                    avint = 1000, # for total intensities
                    bmcv = 0.05, # for total intensities
                    b = c(20,100), # for total intensities
                    l = 1, #amount of background (l = sign/backg) [0.5-2]
                    m = 0.5, #balance between marker (SNP array)
                    ideal_depth = 60,  #ideal read depth
                    mrk_effect = 0.43, # for depth
                    individual_effect = 0.96,# for depth
                    dispersion = 1.8,# for depth
                    seq_e = 0.01){
  #functions
  adjust <- function(mu,rho){
    alpha <- mu * (1-rho)/rho
    beta <- (1 - mu) * (1-rho)/rho
    adjusttt <- rbeta(length(mu),shape1 = alpha,shape2 = beta) #random generation of beta distribution
    return(adjusttt)
  }
  
  #make sure which one you wan to simulate: SNP array OR Sequence reads
  writeLines("########Welcome to use GenoSim########")
  if(!choice %in% c("A","B")) writeLines("Do you want to simulate A.SNP array, B.sequence reads:")
  while(!choice %in% c("A","B")){
    if(count > 2) writeLines("Please check your choice!")
    choice <- readline(prompt="Please make your choice (A or B):")
    count <- count + 1
  }
  
  #set seed 
  if(is.null(seed_number)){
    seed_number = round(runif(1,0,1000))
  }
  
  #ploidy set-up
  if(is.null(ploidy2)){
    ploidy_temp <- replicate(ncol(dsg), ploidy)
  }else{
    ploidy_temp <- c(ploidy, ploidy2, replicate(ncol(dsg) - 2, (ploidy + ploidy2)/2))
  }
  
  
  ##Calculate Genetic attributes: convert dosages into read counts ratios: genotype/ploidy           
  p <- sapply(1:ncol(dsg), function(x) {
    round(dsg[,x]/ploidy_temp[x],3)
  })
  
  #Total intensity  vs. Depth
  if(choice == 'A'){ #SNP array
    writeLines(paste0("You would like to simulate SNP array"))
    #######Simulate Total Intensities#######    
    set.seed(seed_number)
    mav <- rnorm(nrow(dsg), avint, avint*bmcv)
    mcv <- rbeta(1, b[1], b[2])
    n <- matrix(rnorm(nrow(dsg)*ncol(dsg), mav, mav*mcv), nrow = nrow(dsg))
    
    #######Simulate background intensities#######
    b <- rowMeans(n)/(l+1) # background intensity
    s <- n - b             # signal intensity
  } else { #Sequence Reads
    writeLines(paste0("You would like to simulate sequence reads"))
    #######Simulate Read depth#######    
    #with some variaion, negative binomial
    set.seed(seed_number)
    V <- rnorm(nrow(dsg),0,mrk_effect) # marker level effect
    set.seed(seed_number)
    W <- rnorm(ncol(dsg),0, individual_effect)  # individual effect
    
    s_V <- quantile(V,probs = 0.95)
    s_W <- quantile(W,probs = 0.95)
    bt <- log(ideal_depth)      # overall mean
    
    if(bt < s_V |bt < s_W){
      stop("Please choose a higher ideal_depth")
    }
    
    l <- sapply(V, function(x){
      exp(bt + x + W)
    }) 
    
    dep <- t(sapply(1:length(V), function(i){
      rnbinom(n = nrow(l), # per marker
              mu = l[,i], #get the mean for each individual
              size = dispersion) # dispersion in overall
    }))
    
    #######Sequencing error#######    
    p <- p * (1 - seq_e) + (1 - p) * seq_e
  }
  
  #######allelic bias#######
  # pu <- p*ale_b/(p*ale_b+(1-p))  # add allelic bias to p A/B
  set.seed(seed_number)
  aleB <- rbeta(ncol(dsg),ale_b[1],ale_b[2])
  aleBdisset <- t(replicate(nrow(dsg),aleB))
  pu <- p/ (aleBdisset * (1.0 - p) + p) #A/A+B
  
  #######overdispersion#######
  set.seed(seed_number)
  overdis <- rbeta(ncol(dsg),od[1],od[2])
  overdisset <- t(replicate(nrow(dsg),overdis/50))
  para <- matrix(adjust(mu = pu,rho = overdisset),nrow = nrow(dsg))
  
  #######Generate output#######
  if(choice == 'A'){
    set.seed(seed_number)
    Ys <- matrix(rbinom(n = length(n), size = round(s), prob = para),nrow = nrow(dsg))
    #add background signal
    Y <- Ys + matrix(rbinom(n = length(Ys), size = round(b), prob = m),
                     nrow = nrow(Ys))
    # obatin X
    X <- n - Y
    
    colnames(X)  <- individuals <- colnames(Y) <- colnames(dsg)
    rownames(X) <- markers <- rownames(Y) <- rownames(dsg)
    Template <- do.call(rbind, lapply(markers, function(m){data.frame("MarkerName" = m,
                                                                      "SampleName"  = names(dsg[m,]),
                                                                      "Total" = as.numeric(X[m,]) + as.numeric(Y[m,]),
                                                                      "X" = as.numeric(X[m,]),
                                                                      "Y" = as.numeric(Y[m,]),
                                                                      "ratio" = as.numeric(Y[m,])/(as.numeric(X[m,]) + as.numeric(Y[m,])))}))
    output <- list()
    output[["X"]] <- X
    output[["Y"]] <- Y
    output[["GenotypeCallTemplate"]] <- Template
    Parameter <- list()
    Parameter[["BackgroundIntensities"]] <- n
    Parameter[["AllelicBias"]] <- ale_b
    Parameter[["Overdispersion"]] <- overdis
    output[["Parameter"]] <- Parameter
  }else{
    set.seed(seed_number)
    counts <- matrix(rbinom(n = length(dep), size = dep, para),nrow= nrow(dsg))    #binomial model [Yanlin], do not takes account in the overal error because error has been modeled more specifically in each different steps
    
    colnames(counts)  <- individuals <- colnames(dep) <- colnames(dsg)
    rownames(counts) <- markers <- rownames(dep) <- rownames(dsg)
    Template <- do.call(rbind, lapply(markers, function(m){data.frame("MarkerName" = m,
                                                                      "SampleName"  = names(dep[m,]),
                                                                      "ReadDepth" = as.numeric(dep[m,individuals]),
                                                                      "A" = as.numeric(counts[m,individuals]),
                                                                      "a" = as.numeric(dep[m,individuals] - counts[m,individuals]),
                                                                      "ratio" = as.numeric(counts[m,individuals])/as.numeric(dep[m,individuals]))}))
    output <- list()
    output[["ReadDepth"]] <- dep
    output[["ReadCounts"]] <- counts
    output[["GenotypeCallTemplate"]] <- Template
    Parameter <- list()
    Parameter[["SequencingError"]] <- seq_e
    Parameter[["AllelicBias"]] <- ale_b
    Parameter[["Overdispersion"]] <- overdis
    output[["Parameter"]] <- Parameter
  }
  return(output)
}

GenoSim_oneMarker <- function(dsg,
                              choice = "B", # A.SNP array, B.sequence reads
                              ploidy = 4, 
                              ploidy2 = 4,
                              seed_number = 1,
                              od =od,
                              ale_b = ale_b, #allelic bias (ale_b = Y/X) [0.2-0.8]
                              avint = avint, # for total intensities
                              mcv = mcv, # for total intensities
                              l = l, #amount of background (l = sign/backg) [0.5-2]
                              m = m, #balance between marker (SNP array)
                              ideal_depth = 60,  #ideal read depth
                              size = size,# for depth
                              seq_e = 0.01){# for depth
  #functions
  library(MASS)
  adjust <- function(mu,rho){
    alpha <- mu * (1-rho)/rho
    beta <- (1 - mu) * (1-rho)/rho
    adjusttt <- rbeta(length(mu),shape1 = alpha,shape2 = beta) #random generation of beta distribution
    return(adjusttt)
  }
  
  #make sure which one you wan to simulate: SNP array OR Sequence reads
  writeLines("########Welcome to use GenoSim########")
  if(!choice %in% c("A","B")) writeLines("Do you want to simulate A.SNP array, B.sequence reads:")
  while(!choice %in% c("A","B")){
    if(count > 2) writeLines("Please check your choice!")
    choice <- readline(prompt="Please make your choice (A or B):")
    count <- count + 1
  }
  
  #set seed 
  if(is.null(seed_number)){
    seed_number = round(runif(1,0,1000))
  }
  
  #ploidy set-up
  if(is.null(ploidy2)){
    ploidy_temp <- replicate(ncol(dsg), ploidy)
  }else{
    ploidy_temp <- c(ploidy, ploidy2, replicate(ncol(dsg) - 2, (ploidy + ploidy2)/2))
  }
  
  
  ##Calculate Genetic attributes: convert dosages into read counts ratios: genotype/ploidy           
  p <- sapply(1:ncol(dsg), function(x) {
    round(dsg[,x]/ploidy_temp[x],3)
  })
  
  #Total intensity  vs. Depth
  if(choice == 'A'){ #SNP array
    writeLines(paste0("You would like to simulate SNP array"))
    #######Simulate Total Intensities#######    
    set.seed(seed_number)
    # mav <- rnorm(nrow(dsg), avint, avint*bmcv)
    # mcv <- rbeta(1, b[1], b[2])
    # n <- matrix(rnorm(nrow(dsg)*ncol(dsg), mav, mav*mcv), nrow = nrow(dsg))
    n <- matrix(rnorm(ncol(dsg),avint,mcv),nrow = 1)
    
    #######Simulate background intensities#######
    b <- rowMeans(n)/(l+1) # background intensity
    # s <- n - as.numeric(b)             # signal intensity
  } else { #Sequence Reads
    writeLines(paste0("You would like to simulate sequence reads"))
    set.seed(seed_number)
    dep <- matrix(rnbinom(n = ncol(dsg), # per marker
                          mu = ideal_depth, #get the mean for each individual
                          size = size),nrow = 1)
    
    #######Sequencing error#######    
    p <- p * (1 - seq_e) + (1 - p) * seq_e
  }
  
  #######allelic bias#######
  # pu <- p*ale_b/(p*ale_b+(1-p))  # add allelic bias to p A/B
  set.seed(seed_number)
  # aleB <- ale_b
  # aleBdisset_A <- t(replicate(nrow(dsg),aleB)) #A/A+B
  # aleBdisset_B <- 1 - aleBdisset_A + 0.0001
  # aleBdisset <- aleBdisset_A/aleBdisset_B
  # pu <- p/ (aleBdisset * (1.0 - p) + p) #A/A+B
  # pu1 <- p*aleB/(p*aleB+(1-p))
  # if(choice == 'A'){
  pu <- p*ale_b/(p*ale_b+(1-p))
  # }else{
  # pu <- p/ (ale_b * (1.0 - p) + p) #A/A+B
  # }
  
  #######overdispersion#######
  set.seed(seed_number)
  # overdis <- rbeta(ncol(dsg),od[1],od[2])
  overdisset <- t(replicate(nrow(dsg),od/10))
  para <- matrix(adjust(mu = pu,rho = overdisset),nrow = nrow(dsg))
  
  #######Generate output#######
  if(choice == 'A'){
    set.seed(seed_number)
    Xs <- matrix(rbinom(n = length(n), size = floor(n), prob = para),nrow = nrow(dsg))
    #add background signal
    X <- Xs + matrix(rbinom(n = length(Xs), size = as.numeric(round(b)), prob = as.numeric(m)),
                     nrow = nrow(Xs))
    # obatin X
    Ys <- n - Xs
    Y <- Ys + matrix(rbinom(n = length(Ys), size = as.numeric(round(b)), prob = 1- as.numeric(m)),
                     nrow = nrow(Ys))
    
    colnames(X)  <- individuals <- colnames(Y) <- colnames(dsg)
    rownames(X) <- markers <- rownames(Y) <- rownames(dsg)
    
    # if(ale_b < 1){
    Template <- do.call(rbind, lapply(markers, function(m){data.frame("MarkerName" = m,
                                                                      "SampleName"  = names(dsg[m,]),
                                                                      "Total" = as.numeric(X[m,]) + as.numeric(Y[m,]),
                                                                      "X" = as.numeric(Y[m,]),
                                                                      "Y" = as.numeric(X[m,]),
                                                                      "ratio" = as.numeric(X[m,])/(as.numeric(X[m,]) + as.numeric(Y[m,])))}))
    # }else{
    #   Template <- do.call(rbind, lapply(markers, function(m){data.frame("MarkerName" = m,
    #                                                                     "SampleName"  = names(dsg[m,]),
    #                                                                     "Total" = as.numeric(X[m,]) + as.numeric(Y[m,]),
    #                                                                     "X" = as.numeric(X[m,]),
    #                                                                     "Y" = as.numeric(Y[m,]),
    #                                                                     "ratio" = as.numeric(Y[m,])/(as.numeric(X[m,]) + as.numeric(Y[m,])))}))
    #   
    # }
    
    output <- list()
    output[["X"]] <- X
    output[["Y"]] <- Y
    output[["GenotypeCallTemplate"]] <- Template
    Parameter <- list()
    Parameter[["BackgroundIntensities"]] <- n
    Parameter[["AllelicBias"]] <- ale_b
    Parameter[["Overdispersion"]] <- od
    output[["Parameter"]] <- Parameter
  }else{
    set.seed(seed_number)
    counts <- matrix(rbinom(n = length(dep), size = dep, para),nrow= nrow(dsg))    #binomial model [Yanlin], do not takes account in the overal error because error has been modeled more specifically in each different steps
    
    colnames(counts)  <- individuals <- colnames(dep) <- colnames(dsg)
    rownames(counts) <- markers <- rownames(dep) <- rownames(dsg)
    # counts <- dep - counts
    Template <- do.call(rbind, lapply(markers, function(m){data.frame("MarkerName" = m,
                                                                      "SampleName"  = names(dep[m,]),
                                                                      "ReadDepth" = as.numeric(dep[m,individuals]),
                                                                      "A" = as.numeric(dep[m,individuals] - counts[m,individuals]),
                                                                      "a" = as.numeric(counts[m,individuals]),
                                                                      "ratio" = as.numeric(counts[m,individuals])/as.numeric(dep[m,individuals]))}))
    output <- list()
    output[["ReadDepth"]] <- dep
    output[["ReadCounts"]] <- counts
    output[["GenotypeCallTemplate"]] <- Template
    Parameter <- list()
    Parameter[["SequencingError"]] <- seq_e
    Parameter[["AllelicBias"]] <- ale_b
    Parameter[["Overdispersion"]] <- od
    output[["Parameter"]] <- Parameter
  }
  return(output)
}


#read description
description <- read.dcf("./DESCRIPTION")

#get profile URL for author
authors <- unlist(strsplit(description[,"Author"],
                 split=","))

#set css style
css <- paste(".tab-content {
                              margin-bottom: 100px;
                              min-height: 700px;
                              }\n",
             sprintf(".container {
                              background-color: %s
                              }", "#F7F7F7"), ".footer {
                              position: fixed;
                              bottom: 0;
                              left: 0;
                              padding 5px 20px;
                              width: 100%;
                              /* Set the fixed height of the footer here */
                              height: 60px;\n",
                      # sprintf("background-color: %s");
             sprintf("background-color: #6f6f6f;
                              background-image: url('WUR_RGB_standard.png');
                              background-size: contain;
                              background-repeat:no-repeat;
                              background-position:right bottom;
                              list-style: none;
                              }", "#F7F7F7")
             )


#set footer
footer <- fluidRow(
  div(
    class="footer",
    column(2, offset=1,
           a("Wageningen UR",
               target="_top",style="color:white"),
           br(),
           # a("https://www.wur.nl/nl/Onderzoek-Resultaten/Onderzoeksinstituten/plant-research/Plant-Breeding.htm",style="color:white")
           a('maintainer: yanlin.liao@wur.nl')
           # br(),
           # paste0("Version ", description[,"Version"])
           )
  )
)


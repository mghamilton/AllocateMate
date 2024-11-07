#mhamilton@cgiar.org
#Feb 2021

#is.wholenumber function
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol #from https://github.com/ProjectMOSAIC/mosaic/blob/master/R/is.wholenumber.R

#check pedigree function
check.ped <- function(ped) {
  if(sum(colnames(ped) %in% c("ID", "DAM", "SIRE")) != 3) {
    stop("colnames of ped must be: ID, DAM and SIRE")
  }
  
  if(!is.character(ped$ID)) {
    stop("Fields in ped ID, DAM and SIRE must be of class character")
  }
  
  if(!is.character(ped$DAM)) {
    stop("Fields in ped ID, DAM and SIRE must be of class character")
  }
  
  if(!is.character(ped$SIRE)) {
    stop("Fields in ped (ID, DAM and SIRE) must be of class character")
  }
  
  if(sum(ped == "") > 0) {
    stop("Fields in ped must not contain a blank spaces.  Replace with 0 or NA.")
  }
  
}

check.ped2 <- function(ped) {
  tmp <- unique(c(ped[,"DAM"], ped[,"SIRE"]))
  tmp <- tmp[tmp != 0]
  tmp <- tmp[!is.na(tmp)]
  if(sum(!tmp %in% ped[,"ID"])) {
    stop("The ID column in \'ped\' does not contain the identifiers of all DAMs and SIREs")
  }
  rm(tmp)
  
  if (sum(ped[, 1] == 0 | ped[, 1] == "0" | is.na(ped[, 1])) > 0) {
    stop("Missing value in the ID column in \'ped\'")
  }
  
  if (sum(ped[(ped[, 2] != 0 & ped[, 2] != "0" & !is.na(ped[, 2])), 2] %in% 
          ped[(ped[, 3] != 0 & ped[, 3] != "0" & !is.na(ped[, 3])), 3]) > 0) {
    warning("Dams appearing as Sires - selfing in ped")
  }
  
  if (sum(duplicated(ped[, 1])) > 0) {
    stop("Some individuals appear more than once in \'ped\'")
  }
  
  if(sum(
    ((ped[, 2] == 0 | ped[, 2] == "0" | is.na(ped[, 2])) & 
     (ped[, 3] != 0 & ped[, 3] != "0" & !is.na(ped[, 3]))) | 
    ((ped[, 3] == 0 | ped[, 3] == "0" | is.na(ped[, 3])) & 
     (ped[, 2] != 0 & ped[, 2] != "0" & !is.na(ped[, 2]))) 
  ) > 0) {
    stop("If an individual is not a founder (i.e. both parents are unknown), both the SIRE and DAM must be specified in \'ped\'.  It may be necessary to define new founders")
  }
}

check.max_F <- function(max_F) {
  if(max_F > 1 | max_F < 0) {
    stop("max_F must be no less than 0 and no greater than 1")
  }
}

check.method <- function(method) {
  if(!method %in% c("min_F", "assortative")) {
    stop("method must be either min_F or assortative")
  }
}

#Check parents function
check.parents <- function(parents) {
  
  if(sum(colnames(parents) %in% c("ID", "SEX", "EBV", "N_AS_PARENT")) != 4) {
    stop("colnames of parents must be: ID, SEX, EBV and N_AS_PARENT")
  }
  
  parents <- parents[,c("ID", "SEX", "EBV", "N_AS_PARENT")]
  
  if(!is.character(parents$ID)) {
    stop("ID must be of type character")
  }
  
  if(sum(!parents$SEX %in% c("M","F")) > 0) {
    stop("SEX must be of type character comprised of M or F, for male and female respectively")
  }
  
  if(!is.numeric(parents$EBV)) {
    stop("EBV must be of type numeric")
  }
  
  if(sum(!is.wholenumber(parents$N_AS_PARENT)) > 0) {
    stop("N_AS_PARENT must be a vector of type integer")
  }
  
  if(sum(parents$N_AS_PARENT < 1) > 0) {
    stop("N_AS_PARENT must contain integers greater than or equal to 1")
  }
  
  if(sum(is.na(parents$ID)) > 0) {
    stop("ID field of 'parents' must not contain missing values")
  }
  
  if(sum(is.na(parents$SEX)) > 0) {
    stop("SEX field of 'parents' must not contain missing values")
  }
  
  if(sum(is.na(parents$N_AS_PARENT)) > 0) {
    stop("N_AS_PARENT field of 'parents'  must not contain missing values")
  }
  
  if(sum(parents[parents$SEX == "F","N_AS_PARENT"]) !=  sum(parents[parents$SEX == "M","N_AS_PARENT"])) {
    stop("Sum of N_AS_PARENT for females must equal the sum of N_AS_PARENT for males")
  }
}

#Check n_fam_crosses function

check.n_fam_crosses <- function(n_fam_crosses) {
  if(sum(!is.wholenumber(n_fam_crosses)) > 0) {
    stop("n_fam_crosses must be of type integer")
  }
}

#Check H function 
check.H <- function(H) {
  
  if(!is.matrix(H)) {
    stop("H is not a matrix")
  }
  
  rownames(H) <- as.character(rownames(H))
  colnames(H) <- as.character(colnames(H))  
  
  if(sum(!rownames(H) == colnames(H)) > 0) {
    stop("rownames of H do not match colnames of H")
  }
  
  if(max(H[upper.tri(H)]) > 3) {
    stop("One or more off-diagonal elements of H are greater than 3")
  }
  
  if(min(H[upper.tri(H)]) < -1) {
    stop("One or more off-diagonal elements of H are less than -1")
  }
}

reduce.ped <- function(ped, parents) {
  
  ancestors <- parents$ID
  ancestors_prev <- 0
  
  while(ancestors_prev != length(ancestors)) {
    tmp <- ped$ID %in% ancestors
    ancestors_prev <- length(ancestors)
    ancestors <- unique(c(ancestors,as.character(ped[tmp, "DAM"]), as.character(ped[tmp, "SIRE"])))
    rm(tmp)
  }
  
  ped <- ped[ped$ID %in% ancestors,]
  
  return(ped)
}


summarise.fam <- function(families, parents) {
  
  #Summary
  summary <- data.frame(SELECTED = c("N", "Y", "All"),
                        COUNT_FAMS = c(nrow(families[families$SELECTED=='N',]),                 nrow(families[families$SELECTED=='Y',]),                 nrow(families)),
                        MEAN_EBV   = c(mean(families[families$SELECTED=='N',"EBV"], na.rm = T), mean(families[families$SELECTED=='Y',"EBV"], na.rm = T), mean(families[,"EBV"],na.rm = T)),
                        SD_EBV     = c(sd(families[families$SELECTED=='N',"EBV"], na.rm = T),   sd(families[families$SELECTED=='Y',"EBV"], na.rm = T),   sd(families[,"EBV"],na.rm = T)),
                        MIN_EBV    = c(min(families[families$SELECTED=='N',"EBV"], na.rm = T),  min(families[families$SELECTED=='Y',"EBV"], na.rm = T),  min(families[,"EBV"],na.rm = T)),
                        MAX_EBV    = c(max(families[families$SELECTED=='N',"EBV"], na.rm = T),  max(families[families$SELECTED=='Y',"EBV"], na.rm = T),  max(families[,"EBV"],na.rm = T)),
                        MEAN_F     = c(mean(families[families$SELECTED=='N',"F"], na.rm = T),   mean(families[families$SELECTED=='Y',"F"], na.rm = T),   mean(families[,"F"], na.rm = T)),
                        SD_F       = c(sd(families[families$SELECTED=='N',"F"], na.rm = T),     sd(families[families$SELECTED=='Y',"F"], na.rm = T),     sd(families[,"F"], na.rm = T)),
                        MIN_F      = c(min(families[families$SELECTED=='N',"F"], na.rm = T),    min(families[families$SELECTED=='Y',"F"], na.rm = T),    min(families[,"F"], na.rm = T)),
                        MAX_F      = c(max(families[families$SELECTED=='N',"F"], na.rm = T),    max(families[families$SELECTED=='Y',"F"], na.rm = T),    max(families[,"F"], na.rm = T))
  )
  
  crosses <- list(summary      = summary,
                  all_families = families,
                  optimal_families = families[families[,"SELECTED"] == "Y",])
  
  if(sum(parents$N_AS_PARENT)/2 != summary[summary$SELECTED == "Y","COUNT_FAMS"]) {
    stop("Optimal solution not found.  You may need to increase max_F.")
  }
  
  return(crosses)
}

generate.fams <- function(H, parents, ped, max_F) {
  
 # if("dplyr" %in% installed.packages()[, "Package"] == F) {install.packages("dplyr")}   
 #  library(dplyr) 
  
  #Data checks
  check.H(H)
  check.parents(parents)
  check.ped(ped)
  check.max_F(max_F)
  
  if(sum(!parents$ID %in% rownames(H) > 0)) {
    stop("H must contain parents in rownames and colnames")
    
  }
  
  #H matrix with dams in rows and sires in columns
  H <- H[rownames(H) %in% parents[parents$SEX == "F","ID"],
         colnames(H) %in% parents[parents$SEX == "M","ID"]]
  H <- H[!duplicated(rownames(H)),
         !duplicated(colnames(H))]
  
  #Loop to create data frame with possible families
  families <- NULL
  for(row in 1:nrow(H)) {
    for(col in 1:ncol(H)) {  
      families <- rbind(families, cbind(rownames(H)[row], colnames(H)[col], H[row,col]/2))
    }
  }
  colnames(families) <- c("DAM","SIRE", "F")
  families   <- as.data.frame(families)
  families$F <- as.numeric(families$F)
  
  families$FAMILY <- 1:nrow(families) #Create FAMILY ID
  families <- families[c("FAMILY","DAM","SIRE", "F")]
  
  #merge EBV data
  ebvs <- dplyr::left_join(families, parents[,c("ID", "EBV")], by = c("DAM"  = "ID"))
  ebvs <- dplyr::left_join(ebvs, parents[,c("ID", "EBV")], by = c("SIRE" = "ID"))
  ebvs$EBV <- rowMeans(ebvs[,c("EBV.x","EBV.y")], na.rm = T)
  ebvs[is.na(ebvs$EBV.x),"EBV"] <- NA
  ebvs[is.na(ebvs$EBV.y),"EBV"] <- NA
  colnames(ebvs)[colnames(ebvs) == "EBV.x"] <- "dam_ebv"
  colnames(ebvs)[colnames(ebvs) == "EBV.y"] <- "sire_ebv"
  families <- dplyr::left_join(families, ebvs[,c("FAMILY", "dam_ebv", "sire_ebv", "EBV")], by = "FAMILY")
  rm(ebvs)
  
  #parental FAMILY combinations
  tmp <- ped
  tmp[is.na(tmp$DAM), "DAM"] <- paste0("_dam_",1:sum(is.na(tmp$DAM)))
  tmp[is.na(tmp$SIRE), "SIRE"] <- paste0("_sire_",1:sum(is.na(tmp$SIRE)))
  tmp$fam <- as.numeric(as.factor(paste0(tmp$DAM, "_", tmp$SIRE)))  
  tmp <- tmp[,c("ID", "fam")]
  tmp$ID <- as.character(tmp$ID)
  
  colnames(tmp) <- c("DAM", "dam_fam")
  families  <- dplyr::left_join(families , tmp[,c("DAM","dam_fam")], by= "DAM")
  
  colnames(tmp) <- c("SIRE", "sire_fam")
  families  <- dplyr::left_join(families , tmp[,c("SIRE","sire_fam")], by= "SIRE")
  
  families <- families[order(as.numeric(families$FAMILY) , decreasing = FALSE), ] 
  
  families$fam_combn <- paste0(pmin(families$dam_fam, families$sire_fam, na.rm = T), "_", pmax(families$dam_fam, families$sire_fam, na.rm = T))
  families$fam_combn <- as.factor(as.numeric(as.factor(families$fam_combn)))                                                                    
  families$dam_fam   <- as.factor(families$dam_fam)
  families$sire_fam  <- as.factor(families$sire_fam)
  
  families_all <- families
  families <- families[families$F <= max_F,] #remove families with excessive F
  
  tmp1 <- aggregate(!is.na(families[,"DAM"]),by=list(families[,"DAM"]), FUN = "sum") 
  colnames(tmp1) <- c("ID", "N_possible_families_after_max_F_constraint_applied")
  tmp2 <- aggregate(!is.na(families[,"SIRE"]),by=list(families[,"SIRE"]), FUN = "sum") 
  colnames(tmp2) <- c("ID", "N_possible_families_after_max_F_constraint_applied")
  
  parents <- dplyr::left_join(parents, rbind(tmp1, tmp2), by = "ID")
  rm(tmp1, tmp2)
  
  parents[is.na(parents$N_possible_families_after_max_F_constraint_applied),"N_possible_families_after_max_F_constraint_applied"] <- 0
  
  tmp <- parents$ID[parents$N_AS_PARENT > parents$N_possible_families_after_max_F_constraint_applied]
  if(length(tmp) > 0 ) {
    print(parents)    
    stop(paste("max_F too small given N_AS_PARENT values. Check: ",paste(tmp, sep=" ")))
  }
  
  return(families_all)
}

solve_lp <- function(families, parents, n_fam_crosses, max_F, min_trait) {
  
#  if("lpSolveAPI" %in% installed.packages()[, "Package"] == F) {install.packages("lpSolveAPI")}   
#  library(lpSolveAPI) 
  
#  if("dplyr" %in% installed.packages()[, "Package"] == F) {install.packages("dplyr")}   
#  library(dplyr) 
  
  #Data checks
  check.parents(parents)
  check.n_fam_crosses(n_fam_crosses)
  
  fam_combns   <- as.matrix(levels(families$fam_combn))
  N_fam_combns <- length(fam_combns)  
  
  EBV_mean <- sum(parents$EBV * parents$N_AS_PARENT) / sum(parents$N_AS_PARENT)
  families$EBV_dev_squared <- -abs(families$EBV - EBV_mean)^2
  
  #Linear Programming to minimise F#########################################
  #http://www.icesi.edu.co/CRAN/web/packages/lpSolveAPI/vignettes/lpSolveAPI.pdf
  
  print("Creating lpSolve linear program model object")
  mate_lp <- lpSolveAPI::make.lp(nrow(parents)+N_fam_combns, nrow(families))  
  
  #creates an lpSolve linear program model object with nrow(parents) + levels of fam_combn constraints and nrow(families) decision variables 
  
  for(fam in 1:nrow(families)) {
    par_count_temp <- NULL
    
    #Count times par is a parent in fam.  Will equal 0 (neither SIRE nor DAM), 1 (SIRE or DAM) or 2 (if self)
    for(par in 1:nrow(parents)){  
      par_count <- as.matrix((families[fam,"DAM"] == parents[par,1]) + (families[fam,"SIRE"] == parents[par,1]))     
      par_count_temp <- as.matrix(cbind(par_count_temp,par_count))
    }
    par_count_temp <- as.vector(par_count_temp)
    
    par_fam_count_temp <- as.vector(1*(fam_combns == as.numeric(families[fam,"fam_combn"])))    
    
    #vector of counts for the number of times par is a parent in fam
    lpSolveAPI::set.column(mate_lp, fam, c(par_count_temp,par_fam_count_temp))
    
    #Constrain FAMILY count to 0 or 1 (i.e. binary)   
    lpSolveAPI::set.type(mate_lp, fam, "binary") 
  }
  
  lpSolveAPI::set.objfn(mate_lp, as.numeric(families[,min_trait]))
  lpSolveAPI::set.constr.type(mate_lp, c(rep("=",nrow(parents)),rep("<=",N_fam_combns)))
  
  tmp <- rep(n_fam_crosses,length(fam_combns))
  tmp[fam_combns %in% families[families$F > max_F,"fam_combn"]] <- 0 #exclude if F > max_F
  lpSolveAPI::set.rhs(mate_lp, c(parents[,"N_AS_PARENT"],tmp))
  rm(tmp)
  
  dimnames(mate_lp) <- list(c(parents[,"ID"],(-1*as.numeric(fam_combns))),families[,"FAMILY"])
  
  #print("Writing linear program")
  #write.lp(mate_lp, "mate_allocation_linear_program.txt", type = c("lp", "mps", "freemps"), use.names = c(TRUE, TRUE))
  
  print("Solving linear program")
  #Solve linear program
  solved <- solve(mate_lp) #0 indicates that the model was successfully solved.
  
  if(solved != 0) {
    stop ("Linear program not solved.  Try relaxing \'max_F\' constraint (or \'n_fam_crosses\' if using allocate.mate.ped) , or altering \'N_AS_PARENT\' values to maximise the number of possible mating pairs with coefficients of coancestry below max_F.")
  }
  
  selected <- data.frame(lpSolveAPI::get.variables(mate_lp))
  selected$FAMILY <- rownames(selected)
  colnames(selected)[1] <- "SELECTED"
  families <- merge(families, selected, by = "FAMILY",all = FALSE)
  families$SELECTED[families$SELECTED == 0] <- 'N'
  families$SELECTED[families$SELECTED == 1] <- 'Y'
  
  #Sort
  families <- families[order(as.numeric(families$FAMILY) , decreasing = FALSE), ] 
  families <- families[order((families$SELECTED) , decreasing = TRUE),  ] 
  
  families <- families[,c("SIRE",	"DAM",	"F",	"EBV", "SELECTED")]  
  
  #Summary
  crosses <- summarise.fam(families = families, parents = parents)
  
  return(crosses)
  
}


ped.order <- function (pedigree) {
  #order pedigree to ensure each individual is listed in ID is before it is listed as a SIRE or DAM
  pedigree <- pedigree[order(pedigree[,"ID"]),]
  pedigree[,"ID_GEN"] <- NA
  pedigree[(pedigree[, 2] == 0 | pedigree[, 2] == "0" | is.na(pedigree[, 2])) & 
             (pedigree[, 3] == 0 | pedigree[, 3] == "0" | is.na(pedigree[, 3])),"ID_GEN"] <- 0

  iteration <- 0
  nrow_ped  <- nrow(pedigree)
  while(sum(is.na(pedigree[,"ID_GEN"])) > 0) {
    gen_known   <- pedigree[!is.na(pedigree["ID_GEN"]),]
    gen_unknown <- pedigree[is.na(pedigree["ID_GEN"]),]
    tmp <- gen_known[,c("ID", "ID_GEN")]
    
    colnames(tmp) <- c("DAM", "DAM_GEN")
    gen_unknown <- merge(gen_unknown, tmp, by = "DAM", all.x = T)
    colnames(tmp) <- c("SIRE", "SIRE_GEN")
    gen_unknown <- merge(gen_unknown, tmp, by = "SIRE", all.x = T)
    gen_unknown[,"ID_GEN"] <-  (gen_unknown[,"DAM_GEN"] + gen_unknown[,"SIRE_GEN"])/2 + 1
    gen_unknown <- gen_unknown[,colnames(gen_known)]
    pedigree <- rbind(gen_known, gen_unknown)
    rm(tmp, gen_known, gen_unknown)
    
    if(iteration > nrow_ped) {
      stop("Unknown issue when ordering pedigree so that ID listed before it is a DAM or SIRE")
    }
    iteration + 1
  }
  pedigree <- pedigree[order(pedigree[,"ID_GEN"]),]
  
  return(pedigree[,c("ID", "DAM", "SIRE")])
}

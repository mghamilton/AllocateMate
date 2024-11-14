#' allocate.mate.ped
#' 
#' @description 
#' This function generates a mating list for a set of parents.  
#' The mating list can be generated i) to minimise the average inbreeding coefficient (F) of families generated or ii) according to assortative mating principles.
#' Inputs include a list of parents and a 3-column pedigree file specifying the ancestry of these candidates.
#' @param ped is a 3-column data frame with the following columns (class in parentheses):
#' \itemize{
#'  \item{'ID' is the individual identifier of parents and their ancestors (character).} 
#'  \item{'DAM' is the identifier of the individual's dam (NA if unknown) (character).} 
#'  \item{'SIRE' is the identifier of the individual's sire (NA if unknown) (character).} 
#' }
#' @param parents data frame with the following columns (class in parentheses):
#' \itemize{
#'  \item{'ID' is the individual identifier (character).} 
#'  \item{'SEX' is the sex of the individual - 'M' or 'F', for male and female respectively (character).} 
#'  \item{'EBV' Estimated breeding value (EBV) of the individual - can not be NA if applying assortative mating as the method (numeric).} 
#'  \item{'N_AS_PARENT' The number of families the indivdiual is to contribute to. Can be zero if backup parents are to be included in the mating list output (integer).} 
#'  \item{'FAM' is the full-sibling family identifier. This column may be omitted (character).} 
#' }
#' @param max_F is the maximum inbreeding coefficient allowed in the offspring of parents (numeric between 0 and 1)
#' @param method either 'min_F' (to minimise the average inbreeding in offspring) or 'assortative' (to apply assortative mating) (character)
#' @param n_fam_crosses is the maximum number of siblings from one family allowed to be crossed with siblings of another family.  This represents a constraint on the generation of 'double first cousins'. (integer)
#' @return 'summary' is a data frame containing a summary of all possible families generated from crosses between parents: 
#' \itemize{
#'  \item{'SELECTED' N represents unselected families, Y represents selected families (i.e. mate allocations to be generated) and All represents all possible families.}           
#'  \item{'COUNT_FAMS' count of families.}
#'  \item{MEAN_EBV' mean of family EBVs.}
#'  \item{'SD_EBV' standard deviation of family EBVs.}
#'  \item{'MIN_EBV' minimum of family EBVs.} 
#'  \item{'MIN_EBV' maximum of family EBVs.}  
#'  \item{MEAN_EBV' mean of family inbreeding coefficients (F).}
#'  \item{'SD_EBV' standard deviation of family inbreeding coefficients (F).}
#'  \item{'MIN_EBV' minimum of family inbreeding coefficients (F).} 
#'  \item{'MIN_EBV' maximum of family inbreeding coefficients (F).}  
#' } 
#' @return 'all_families' is a data frame containing details of all possible families able to be generted from matings between parents: 
#' \itemize{
#'  \item{'SIRE' male parent.}           
#'  \item{'DAM' female parent.} 
#'  \item{'F' inbreeding coefficient of family members (i.e. the 'family F').}       
#'  \item{'EBV' mean of parent EBVs (i.e. the 'family EBV').}
#'  \item{'SELECTED' N represents an unselected family (i.e. family is not to be generated), Y represents a selected family (i.e. mate allocated - family to be generated).}           
#' } 
#' @return 'optimal_families' is a data frame containing details of mate allocations (i.e.families to be generated from matings between parents): 
#' \itemize{
#'  \item{'SIRE' male parent.}           
#'  \item{'DAM' female parent.} 
#'  \item{'F' inbreeding coefficient of family members (i.e. the 'family F').}       
#'  \item{'EBV' mean of parent EBVs (i.e. the 'family EBV').}          
#' } 
#' @examples
#' #Retrieve example data
#' ped <- AllocateMate::ped
#' parents <- AllocateMate::parents
#' 
#' output <- allocate.mate.ped(ped, parents, max_F = 0.02, method = "min_F", n_fam_crosses = 1) 
#' output$summary
#' head(output$all_families)
#' head(output$optimal_families)
#' @import lpSolveAPI
#' @import AGHmatrix
#' @import dplyr
#' @export
# @import nadiv

allocate.mate.ped <- function(ped, parents, max_F = 1, method = "min_F", n_fam_crosses = 1) {
  #mhamilton@cgiar.org
  #Feb 2021
  
  #if("nadiv" %in% installed.packages()[, "Package"] == F) {install.packages("nadiv")}   
  #library(nadiv)
  #library(lpSolveAPI)
  #library(AGHmatrix)
  #library(dplyr)
  
  #split parents
  if(!"N_AS_PARENT" %in% colnames(parents)) {
    stop("Column N_AS_PARENT is not present in \'parents\'.")
  } 
  
  if(!is.integer(parents$N_AS_PARENT)) {
    stop("N_AS_PARENT must be a vector of type integer")
  } 
  
  if(sum(is.na(parents$N_AS_PARENT)) != 0) {
    stop("N_AS_PARENT contains missing values")
  } 
  
  all_candidates <- parents
  parents <- parents[parents$N_AS_PARENT > 0,]
  
  #  if(nrow(parents) == nrow(all_candidates)) { 
  #    all_candidates <-  NULL
  #  }
  
  check.parents(parents)  
  check.ped2(ped)
  ped <- reduce.ped(ped = ped, parents = all_candidates)
  check.n_fam_crosses(n_fam_crosses)
  check.max_F(max_F)
  check.method(method)
  
  # if(!is.null(all_candidates)) {
  all_candidates <- check.all_candidates(ped, parents, all_candidates)
  #  }
  
  #ped <- nadiv::prunePed(ped = ped, phenotyped = parents$ID)
  ped[ped$DAM  == 0 & !is.na(ped$DAM), "DAM"]  <- NA
  ped[ped$SIRE == 0 & !is.na(ped$SIRE),"SIRE"] <- NA
  # ped <- nadiv::prepPed(ped = ped)
  ped <- ped.order(pedigree = ped)
  ped[is.na(ped$DAM), "DAM"]  <- 0
  ped[is.na(ped$SIRE),"SIRE"] <- 0
  
  # H <- makeA(ped)  #unstable with large pedigrees, use Amatrix instead
  H <- AGHmatrix::Amatrix(ped)
  H <- H[rownames(H) %in% parents$ID, colnames(H) %in% parents$ID]
  H <- as.matrix(H)
  
  families <- generate.fams(H = H, ped = ped, parents = parents, max_F = max_F) 
  
  if(method == "assortative") {
    output <- solve_lp(families = families, parents = parents, n_fam_crosses = n_fam_crosses, max_F = max_F, min_trait = "EBV_dev_squared")
  }
  
  if(method == "min_F") {
    output <- solve_lp(families = families, parents = parents, n_fam_crosses = n_fam_crosses, max_F = max_F, min_trait = "F")
  }
  
  output$mating_list <- get_optimal_all_candidates(optimal_families = output$optimal_families, 
                                                   all_candidates = all_candidates)
  output$mating_list <- left_join(output$mating_list, 
                                  output$optimal_families[,c("SIRE", "DAM", "F", "EBV")], 
                                  by = c("SIRE", "DAM"))
  
  crosses <- data.frame(CROSS = unique(output$mating_list$CROSS),
                        CROSSES_WITH_SIRE_FAM = NA,
                        CROSSES_WITH_DAM_FAM = NA)
  
  for(cross in crosses$CROSS) {
    tmp_fam     <-  unique(output$mating_list[output$mating_list$CROSS == cross,"SIRE_FAM"])
    tmp_crosses <- output$mating_list[output$mating_list$SIRE_FAM == tmp_fam |
                                        output$mating_list$DAM_FAM == tmp_fam, "CROSS"]
    crosses[crosses$CROSS == cross, "CROSSES_WITH_SIRE_FAM"]  <-  paste(tmp_crosses[order(tmp_crosses)],  collapse = " ")
    tmp_fam     <-  unique(output$mating_list[output$mating_list$CROSS == cross,"DAM_FAM"])
    tmp_crosses <- output$mating_list[output$mating_list$SIRE_FAM == tmp_fam |
                                        output$mating_list$DAM_FAM == tmp_fam, "CROSS"]
    crosses[crosses$CROSS == cross, "CROSSES_WITH_DAM_FAM"]  <-  paste(tmp_crosses[order(tmp_crosses)],  collapse = " ")
  rm(tmp_fam, tmp_crosses)
  }
  output$mating_list <- left_join(output$mating_list, crosses)
  rm (crosses, cross)
  
  output$A_matrix <- H[rownames(H) %in% c(output$optimal_families$SIRE, output$optimal_families$DAM), 
                       colnames(H) %in% c(output$optimal_families$SIRE, output$optimal_families$DAM)]
  
  return(output)
  
}
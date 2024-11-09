#' allocate.mate.H
#' 
#' @description 
#' This function generates a mating list for a set of parents.  
#' The mating list can be generated i) to minimise the average inbreeding coefficient (F) of families generated or ii) according to assortative mating principles.
#' Inputs include a list of parents and a relationship matrix (H) including these parents.
#' @param H is square additive genetic relationship matrix with individual idenifiers as rownames and colnames.
#' @param parents data frame with the following columns (class in parentheses):
#' \itemize{
#'  \item{'ID' is the individual identifier (character).} 
#'  \item{'SEX' is the sex of the individual - 'M' or 'F', for male and female respectively (character).} 
#'  \item{'EBV' Estimated breeding value (EBV) of the individual - can not be NA if applying assortative mating as the method (numeric).} 
#'  \item{'N_AS_PARENT' The number of families the indivdiual is to contribute to (integer).} 
#' }
#' @param max_F is the maxiumum inbreeding coefficient allowed in the offspring of parents (numeric between 0 and 1)
#' @param method either 'min_F' (to minimise the average inbreeding in offspring) or 'assortative' (to apply assortative mating) (character)
#' @return 'summary' is a data frame containing a summary of all possible families generated from matings between male parents: 
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
#' @return 'all_families' is a data frame containing details of all possible families able to be generated from matings between parents: 
#' \itemize{
#'  \item{'SIRE' male parent.}           
#'  \item{'DAM' female parent.} 
#'  \item{'F' inbreeding coefficient of family members (i.e. the 'family F').}       
#'  \item{'EBV' mean of parent EBVs (i.e. the 'family EBV').}
#'  \item{'SELECTED' N represents an unselected family (i.e. family is not to be generated), Y represents a selected family (i.e. mate allocated - family to be generated).}           
#' } 
#' @return 'optimal_families' is a data frame containing details of mate allocations (i.e.families to be generated from matings between male parents): 
#' \itemize{
#'  \item{'SIRE' male parent.}           
#'  \item{'DAM' female parent.} 
#'  \item{'F' inbreeding coefficient of family members (i.e. the 'family F').}       
#'  \item{'EBV' mean of parent EBVs (i.e. the 'family EBV').}          
#' } 
#' @examples
#' #Retrieve example data
#' H <- AllocateMate::H
#' parents <- AllocateMate::parents
#' 
#' output <- allocate.mate.H(H, parents, max_F = 0.02, method = "min_F") 
#' output$summary
#' head(output$all_families)
#' head(output$optimal_families)
#' @import lpSolveAPI
#' @import dplyr
#' @export
# @import nadiv

allocate.mate.H <- function(H, parents, max_F = 1, method = "min_F") {
  #mhamilton@cgiar.org
  #Feb 2021
  
  # if("nadiv" %in% installed.packages()[, "Package"] == F) {install.packages("nadiv")}   
  # library(nadiv)
  # library(lpSolveAPI)
  # library(dplyr)
  
  check.H(H)
  check.parents(parents)
  check.max_F(max_F)
  check.method(method)
  
  #add ped 
  ped <- data.frame(ID = rownames(H),
                    DAM  = 1:nrow(H),
                    SIRE = (nrow(H)+1):(nrow(H)*2))
  ped[ped$DAM == 0,"DAM"] <- NA
  ped[ped$SIRE == 0,"SIRE"] <- NA
  #ped <- nadiv::prepPed(ped)
  
  ped$ID <- as.character(ped$ID)
  ped$DAM <- as.character(ped$DAM)
  ped$SIRE <- as.character(ped$SIRE)
  
  families <- generate.fams(H = H, parents = parents, ped = ped, max_F = max_F) 
  
  if(method == "assortative") {
    output <- solve_lp(families = families, parents = parents, n_fam_crosses = 1, max_F = max_F, min_trait = "EBV_dev_squared")
  }
  
  if(method == "min_F") {
    output <- solve_lp(families = families, parents = parents, n_fam_crosses = 1, max_F = max_F, min_trait = "F")
  }
  
  return(output)
}

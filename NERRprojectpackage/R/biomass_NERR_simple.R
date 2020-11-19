#' Internal function to produce SSD and total biomass
#'
#' This is the internal function that does all the work of producing the biomass and SSD for each site that is fed back to run_biomass_NERR_simple.
#' @param Params This function is written such that it only takes one set of Mjuv, M, Linf, and k at a time for a given site.  Params is still a list
#' with Params$Mjuv, Params$M, Params$Linf, Params$k, but they only contain one value each.  This function also requires Params$veclength
#' and Params$x, which are created in run_biomass_NERR_simple.
#' @keywords NERR
#' @export
#' @examples
#' biomass_NERR_simple(Params)

biomass_NERR_simple <- function(Params){

  # Do SSD and biomass calculations using IPM for GTM NERR project
  # Started by JW White, April 2019
  # Edited LS Storch Nov 2019, May 2020 for NERR project
  # WARNING - there are a lot of hardcoded variables/parameters in here

 # library(pracma)
#  source("kernmatsimple.R")
#  source("mkkernsimple.R")

  #PARAMETERS - a lot of these were taken from oyster_PP_params or EPR_NERR
  Params$T <- 1 #time step
  Params$F <- 0 #harvest rate - no fishing currently in model
  Params$dx <- diff(Params$x[1:2])
  Params$LW_afdw <- 5.09e-5*Params$x^2.365 # Ash Free Dry Weight; units = g ; derived from Kimbro field samples in salinity zone 2 (taken from oyster_PP_params)
  Params$Rmean <- 7 # mean spat size in week 1 - taken from EPR_NERR
  Params$Rstd <- 5 # std in spat size - taken from EPR_NERR
  Params$Rvec <- dnorm(Params$x,Params$Rmean,Params$Rstd) #taken from EPR_NERR, initializing recruitment vector
  Params$j_length <- 15 #juvenille mortality size limit, taken from oyster_PP_params"
  # Params$Fec <- 19.86e6*Params$LW_afdw^(1.17) #fecundity parameter, taken from oyster_PP_params - no longer used because looking at biomass instead of EPR
  # Params$Fec = Params$Fec * (Params$x>=35) # Impose size at maturity of 35 mm - no longer used because looking at biomass instead of EPR
  # Params$veclength <- 200 #number of patches in the vector, made 200 so the size bins aren't too coarse grain - this is an input now
  # Params$x <- linspace(0,160,Params$veclength) # size bins vector, max value is 160 because that's double the biggest Linf - this is an input now

  # S <- matrix(data = NA,nrow=1,ncol=Params$T) #S <- nan(1,Params$T) # reef structure, currently unused
  # N <- matrix(data = NA, nrow = Params$veclength, ncol = Params$T) #N <- nan(Params$veclength,Params$T) # oyster population size matrix (semi-annual time steps, based on incoming parameters M and k)
  # R <- matrix(data = NA, nrow = Params$veclength, ncol = Params$T) #R <- nan(Params$veclength,Params$T) # recruitment matrix
  # Harv <- matrix(data = NA, nrow = Params$veclength, ncol = Params$T) #Harv <- nan(Params$veclength,Params$T)# oysters harvested (semiannual), currently unused
  # Params$lambdaTAF <- 0.1/52 # Based on annual rate of 0.1 from Powell et al. (2012)
  # Params$density <- 0.849 # grams/cm^3. From aqua-calc.com

  #create growth/mortality kernel
  kmat <- kernmatsimple(Params$x,Params$F,Params,Params$T)# get kernel

  Rdist <- Params$Rvec # recruit vector

  #% Ensure that Rdist integrates to 1:
  Rdist <- Rdist/(sum(Rdist*Params$dx));

  #% Find stable size distribution:
  TT <- 100 #how many times to iterate?
  N <- repmat(cbind(Rdist),1,TT);
  A <- N

  for (t in 2:TT){

    A[,t] = (kmat%*%N[,t-1])*Params$dx # integrate over growth/mortality kernel for each time step

    N[,t] =  A[,t] + Rdist #add in recruits for final population at a time step
  }

  #% Stable size distribution
  SSD <- N[,TT];

  #% Calcluate total biomass by multiplying SSD by the ash free dry weight
  biomass <- sum(Params$LW_afdw*SSD)

  result <- list(SSD,biomass) # R doesn't allow multiple outputs for a function so have to aggregate all results in a list
  return(result)

}







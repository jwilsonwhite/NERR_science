#' NERR function to produce SSD and total biomass
#'
#' Given a set of parameters per site, Mjuv (juvenile mortality), Madult (post-juvenile mortality), Linf (max size), and k (growth rate), this function will produce
#' the EPR and stable size distribution per site and graph the results with the 5th and 95th percentile uncertainties
#' as upper and lower bounds
#'
#' @param Params Params is a list with Params$SiteNames, Params$Mjuv, Params$Madult, Params$k, Params$Linf for all sites.  These can be uploaded by loading Params.Rdata
#' @param CovMats CovMats is a list with the covariance matrices for Mjuv, Madult, Linf, and k (in that order) for the sites.  These can be uploaded by loading CovMats.Rdata
#' @param n n is the number of random parameter draws to create uncertainty bounds
#' @keywords NERR
#' @import pracma
#' @import MASS
#' @import ggplot2
#' @import egg
#' @import grDevices
#' @return A list. The elements are SizesGG and BiomassGG, both ggplot2 objects will plot the size distributions and biomass-per-recruit for each site. SSD has the expected stable size distribution of each site, and Biomass is a dataframe of the biomass-per-recruit at each site.
#' @export
#' @examples
#' run_biomass_NERR_simple(Params,CovMats,n=1000)
#'
run_biomass_NERR_simple <- function(Params,CovMats,n=1000){


  #THINGS THAT NEED TO BE DONE - REDO CALCULATIONS OF M AND MJUV SO THEY COME FROM THE DATA INSTEAD OF MAKING M 10% OF MJUV.
  #MAKE THE DOCUMENTATION IN THE FUNCTIONS CONSISTENT BECAUSE THEY ALL REFER TO M WHEN THAT IS ACTUALLY NOW MJUV AND THE ADULT MORTALITY IS M, ESTIMATED AS 10% OF MJUV
  #ONCE M AND MJUV ARE BOTH CALCULATED FROM THE DATA YOU HAVE TO REDO THE COVARIANCE MATRICES AND THE RANDPARAMS STUFF TO INCLUDE BOTH M AND MJUV AS OUTPUTS TO RANDPARAMS
  #ONCE M AND MJUV ARE REDONE AND BOTH INPUTS YOU HAVE TO COMB OVER THIS FUNCTION AND MAKE SURE YOU'RE NOT DOING, E.G., M <- MJUV*0.1 ANYWHERE

  #LSS December 2019
  #making another function that can deal with the struct of Params for the 7 sites since I restructured EPR_NERR_LS to only take in the parameters for one site at a time
  #here is also where the uncertainty analysis is done

  #the data for this has to be a struct (list) of params: Params$k, Params$Linf, Params$M , Params$sites and also the covariance matrices for each site, CovMats
  #NEW - Params$veclength is now also an input!!!!!


 # source("biomass_NERR_simple.R")
#  source("kernmatsimple.R")
#  source("mkkernsimple.R")
  # The following should be moved to an input option with default = 200
  Params$veclength <- 200 #moved this to the outside of the function because lower down I have to specify the length of the vectors storing the SSD outputs
  Params$x <- seq(from=0,to=160,length.out=Params$veclength) # space vector, max value is 160 because that's double the biggest Linf

  #generating n random parameter combos to run in the model to get upper and lower uncertainty bounds
  RandParams <- vector('list',length(Params$sites))
  for (i in 1:length(Params$sites)){ #need n number of random combos for each site
    invec <- c(Params$Mjuv[i],Params$Linf[i],Params$k[i]) #order of params in the covariance matrices is M, Linf, k, so mimic this in the order of params going into mvrnorm (M is now Mjuv)
    out <- mvrnorm(n, invec, CovMats[[i]][[1]], tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

    #out <- mvrnorm(n, invec, matrix(0,3,3), tol = 1e-6, empirical = FALSE, EISPACK = FALSE) #testing, delete later
    RandParams[[i]] <- out
    dims <- dim(RandParams[[i]])
    for (ii in dims[1]){ #go through all of the rows (dims[1] is number of rows, dims[2] is number of cols)
      if (sum(RandParams[[i]][ii,]>0) > 0){ #if there is a negative value for a parameter in a row, overwrite with NaNs
        RandParams[[i]][ii,] <- NaN

      }
    }

  }

  # Note that at present we have good estimates of juvenile mortality (1st 3 months), but poor estimates of later mortality.
  # As a first approximately, older mortality rates appear to be ~10% of juvenile mortality rates.

  #Now that we have all of the parameters and random parameter combinations, have to multiply M and k by 30 because we have a daily rate and want a monthly
  #Column 1 of RandParams[[i]] is M and column 3 is k
  for (i in 1:length(Params$sites)){
    Params$M[i] <- Params$M[i]*30
    Params$k[i] <- Params$k[i]*30
    Params$Mjuv[i] <- Params$Mjuv[i]*30
    RandParams[[i]][,1] <- RandParams[[i]][,1]*30
    RandParams[[i]][,3] <- RandParams[[i]][,3]*30
  }

  #get the outputs for the 1000 random parameter combos for each site and find the upper and lower 5 and 95 percentile
  RandSSD <- vector('list',length(Params$sites)) #save all the random SSDs
  RandBio <- vector('list',length(Params$sites)) #save all the random biomasses
  SSDmeans <- vector('list',length(Params$sites)) #save the means for all of the random SSDs, using for finding 5th and 95th percentile
  Bio_highlow <- vector('list',length(Params$sites)) #save 5th and 95th percentile biomass
  SSD_highlow <- vector('list',length(Params$sites)) #save 5th and 95th percentile SSD
  SSD_means_highlow <- vector('list',length(Params$sites)) #save 5th and 95th percentile SSD
  for (i in 1:length(Params$sites)){ #outer loop goes through the sites

    RandBio[[i]] <- vector('numeric',n)
    SSDmeans[[i]] <- vector('numeric',n)
    Bio_index <- vector('numeric',2)
    RandSSD[[i]] <- matrix(0,Params$veclength,n) #column-wise is each SSD output for each of the random parameter combos
    SSD_highlow[[i]] <- matrix(0,Params$veclength,2) #save 5th and 95th percentile SSDs
    miniParams <- list()

    for (j in 1:n){ #inner loop goes through the n random parameter combos for each site

      miniParams$Mjuv <- RandParams[[i]][j,1]
      miniParams$M <- RandParams[[i]][j,1]*0.1 #adult mortality approximately 10% of juvenile mortality
      miniParams$Linf <- RandParams[[i]][j,2]
      miniParams$k <- RandParams[[i]][j,3]
      miniParams$veclength <- Params$veclength
      miniParams$x <- Params$x
      out <- biomass_NERR_simple(miniParams)
      if (!is.na(out[[2]])){
        RandBio[[i]][j] <- out[[2]]
      }
      if (!anyNA(out[[1]])){
        RandSSD[[i]][,j] <- cbind(out[[1]])
        SSDmeans[[i]][j] <- apply(cbind(RandSSD[[i]][,j]),2,mean)
      }

    }

    #upper and lower bounds on confidence, grab 5th and 95th percentile biomass
    Bio_highlow[[i]][1] <- quantile(RandBio[[i]],0.95,na.rm=TRUE)
    Bio_highlow[[i]][2] <- quantile(RandBio[[i]],0.05,na.rm=TRUE)
    Bio_diff <- abs(RandBio[[i]] - Bio_highlow[[i]][1])
    Bio_index[1] <- which(Bio_diff == min(Bio_diff))
    Bio_diff <- abs(RandBio[[i]] - Bio_highlow[[i]][2])
    Bio_index[2] <- which(Bio_diff == min(Bio_diff))

    SSD_highlow[[i]][,1] = RandSSD[[i]][,Bio_index[1]] #grab the SSD from the 95th
    SSD_highlow[[i]][,2] = RandSSD[[i]][,Bio_index[2]] #grab the SSD from the 5th
    SSD_means_highlow[[i]][1] <- SSDmeans[[i]][Bio_index[1]]
    SSD_means_highlow[[i]][2] <- SSDmeans[[i]][Bio_index[2]]

  }

  #now that we have the 5th and 95th percentile, lets get the actual SSD and biomasses using the actual parameters so we can plot them with the 5th and 95th
  Biomass <- vector('numeric',length(Params$sites))
  SSD <- vector('list',length(Params$sites))
  maxSSD_highlow <- vector(length=length(Params$sites))
  for (i in 1:length(Params$sites)){ #outer loop goes through the sites
    miniParams$Mjuv <- Params$Mjuv[i]
    miniParams$M <- Params$M[i]
    miniParams$Linf <- Params$Linf[i]
    miniParams$k <- Params$k[i]
    miniParams$veclength <- Params$veclength
    out <- biomass_NERR_simple(miniParams)
    Biomass[i] <- out[[2]]
    SSD[[i]] <- out[[1]]

    # also find max value of SSD_highlow, for plotting
    maxSSD_highlow[i] = max(SSD_highlow[[i]][,1])
  }

  #Plotting - we will want to plot a given SSD[[i]] against SSD_highlow[[i]][,1] and SSD_highlow[[i]][,2] and compare EPR[i] against EPR_highlow[[i]][1] and EPR_highlow[[i]][2]
  # Loop over each site, create a ggobject, then plot them up in a grid
  GGp <- list() # pre-allocate
  Ymax <- max(maxSSD_highlow)*1.001
  for (i in 1:length(Params$sites)){

    # create a data frame for plotting SSD
    Data.sub <- data.frame(SSD=SSD[[i]],Length=Params$x,Hi=SSD_highlow[[i]][,1],Lo=SSD_highlow[[i]][,2])

    GGp[[i]] <- ggplot(data=Data.sub,aes(y=SSD,x=Length))+
      geom_line()+
      geom_line(aes(x=Length,y=Lo))+
      geom_line(aes(x=Length,y=Hi))+
      geom_ribbon(aes(x=Length,ymin=Lo,ymax=Hi),color='blue',fill='blue',alpha=0.5)+
      xlab('Length (mm)')+
      scale_x_continuous(limits=c(0,100))+
      scale_y_continuous(breaks=NULL,limits=c(0,Ymax))+
      ggtitle(Params$sites[i])+
      theme_bw()


  } # end loop over sites



  # Separately, plot biomass:
  Biomass.sub <- data.frame(Site=Params$sites,Biomass=Biomass,Biomassmin = NA,Biomassmax=NA)
  for (i in 1:length(Params$sites)){
    Biomass.sub$Biomassmax[i] = Bio_highlow[[i]][1]
    Biomass.sub$Biomassmin[i] = Bio_highlow[[i]][2]
  }
  # rescale
  Max <- max(Biomass.sub$Biomassmax)
  Biomass.sub$Biomass = Biomass.sub$Biomass/Max
  Biomass.sub$Biomassmin = Biomass.sub$Biomassmin/Max
  Biomass.sub$Biomassmax = Biomass.sub$Biomassmax/Max


  Biomassgg <- ggplot(data=Biomass.sub,aes(x=Site,y=Biomass))+
    geom_linerange(aes(ymin=Biomassmin,ymax=Biomassmax))+
    geom_point()+
    xlab(NULL)+
    scale_x_discrete(guide = guide_axis(angle = 45))+
    ylab('Relative total biomass')+
    theme_bw()



  #scale_x_discrete(limits=c('T','G','St. A','SR','B','M','P'))+

  # Display all results
#  ggarrange(GGp[[1]],GGp[[2]],GGp[[3]],GGp[[4]],GGp[[5]],GGp[[6]],GGp[[7]],ncol = 2)
 # ggarrange(GGp[[1]],GGp[[2]],GGp[[3]],GGp[[4]],GGp[[5]],GGp[[6]],GGp[[7]], ncol = 2)

return(list("SizesGG" = GGp,"BiomassGG" = Biomassgg, "SSD" = SSD, "Biomass" = Biomass.sub))


}



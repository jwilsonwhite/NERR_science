#' Internal middle man function to help create kernel
#'
#' This is an internal middle man function that creates x and y matrices that feed into mkkernsimple.  Not meant to be run on its own.
#' @param x This is Params$x from EPR_NERR_simple
#' @param F This is the fishing rate, set to 0 in EPR_NERR_simple
#' @param T This is the "step size", set to 1 in EPR_NERR_simple
#' @param fixparm This is Params that has Params$Linf, Params$k (only one value for each) but there are other components added in EPR_NERR_simple before Params goes into this function
#' @keywords NERR kernel helper function
#' @export
#' @examples
#' kernmatsimple(x,F,fixparm,T)

kernmatsimple <- function(x,F,fixparm,T){

  #source("mkkernsimple.R")

  #% Set up the integration mesh kernel using Simpsons
  #% For Oyster IPM

  #%adapted from Will from Easterling. Evenly spaced grid, now add weights so
  #%use Simpson's rule using Marissa's code.

  y = x;
  #%this creates a vector (y) that is equal to x

  tmp <- meshgrid(x,y); #% Matlab built in function
  x <- tmp[[1]]
  y <- tmp[[2]]
  #%x is an array (original x by original x) with each row equal to original
  #%vector x
  #%y is an array (original y by original y) with each column equal to
  #%original vector y
  #%X corresponds to size at time t
  #%Y corresponds to size at time t+1

  #% Make the kernel from the grid
  kmat <- mkkernsimple(x,y,F,fixparm,T);

  kmat <- pmax(kmat,0);
}

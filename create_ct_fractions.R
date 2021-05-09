
#' Simulate cell-type fractions
#'
#' This function simulates a data frame of cell-type fractions for 
#' myloid (mean = 0.72) and lyphoid cells (mean = 0.28) using the uniform 
#' distribution, sensible probability cutoffs (i.e cannot exceed 0 to 1), 
#' and known average fractions for each cell type. 
#' Note, the cell type fractions must sum to 1.
#'
#' @param n_samp number of samples to generate (will be doubled cases/controls), default = 100
#' @param seed random seed to set, default = 215
#' @return dataframe of cell-type fractions for each sample
#' @importFrom stats "runif"
#' @export

create_ct_fractions = function(n_samp = 100, seed = 215)
{
  set.seed(seed)
  
  # make sure n_samp is a positive number (greater than zero)
  try(log(n_samp))
  
  myeloid_frac = runif(2*n_samp, 0.44, 1) # mean = .72
  lymphoid_frac = 1 - myeloid_frac # mean = .28
  return(data.frame(myeloid_frac, lymphoid_frac))
}


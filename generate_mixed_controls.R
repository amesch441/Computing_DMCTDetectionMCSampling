
#' Generate Mixture Samples of Cell-types for Control Subjects
#' 
#' This function generates a mixture of the DNAm profiles using the fraction
#' of each cell type that is provided from the function: create_ct_fractions()
#'
#' @param n_samp number of samples to generate, default = 100
#' @param n_dmcts number of DMCTs to select, default = 1000
#' @param seed random seed to set, default = 215
#' @param mc number of monte-carlo runs, default = 10
#' @param lymphoid_frac vector of fractions for lymphoid cells
#' @param myeloid_frac vector of fractions for myeloid cells
#' @param sample_controls_lymphoid list of matrices for the lymphoid specific controls
#' @param sample_controls_myeloid list of matrices for the myeloid specific controls
#' @return returns list of matrices containing mixed samples of cell-types for the control subjects
#' @export

generate_mixed_controls = function(n_dmcts = 1000, n_samp = 100, seed = 215, mc = 10, lymphoid_frac, myeloid_frac, sample_controls_lymphoid, sample_controls_myeloid)
{
  # make sure n_samp is a positive number (greater than zero)
  try(log(n_samp))
  
  sample_controls_mixed = matrix(NA,n_dmcts, n_samp)
  sample_controls_mixed = rep(list(sample_controls_mixed),mc) #contains 1000 x n samples for each MC sample
  
  lymphoid_frac_case = lymphoid_frac[1:n_samp]
  lymphoid_frac_control = lymphoid_frac[(n_samp+1):(2*n_samp)]
  
  myeloid_frac_case = myeloid_frac[1:n_samp]
  myeloid_frac_control = myeloid_frac[(n_samp+1):(2*n_samp)]
  
  for(j in 1:mc) #  MC samples
  {
    for(i in 1:n_samp)
    {
      sample_controls_mixed[[j]][,i] = lymphoid_frac_control[i]*sample_controls_lymphoid[[j]][,i] + myeloid_frac_control[i]*sample_controls_myeloid[[j]][,i] 
      
    }
  }
  return(sample_controls_mixed)
}








#' Generate n Samples for Control Subjects
#' 
#' This function generates DNAm profile samples for the control subjects
#'
#' @param n_samp number of samples to generate, default = 100
#' @param n_dmcts number of DMCTs to select, default = 1000
#' @param seed random seed to set, default = 215
#' @param mc number of monte-carlo runs, default = 10
#' @param beta_cellTypes matrix of beta measures for the given cell-type (myeloid or lymphoid)
#' @param DMCTs list of DMCTs generated from function: sample_DMCTs()
#' @return returns list of matrices containing samples for the control subjects
#' @export

generate_controls = function(n_dmcts = 1000, n_samp = 100, seed = 215, mc = 10, beta_cellTypes, DMCTs)
{
  set.seed(seed)
  replace_q = F
  
  # make sure n_samp is a positive number (greater than zero)
  try(log(n_samp))
  
  if(n_samp > ncol(beta_cellTypes))
  {
    replace_q = T
  }
  
  draw1_control_cellType_index = sample(1:ncol(beta_cellTypes), n_samp, replace = replace_q)
  sample_controls_cellTypes_n_samp = matrix(NA,n_dmcts, n_samp)
  sample_controls_cellTypes_n_samp = rep(list(sample_controls_cellTypes_n_samp),mc) #contains n_dmcts x n samples for each MC sample
  
  
  for(j in 1:mc) # MC samples
  {
    for(i in 1:n_samp)
    {
      sample_controls_cellTypes_n_samp[[j]][,i] = beta_cellTypes[DMCTs,draw1_control_cellType_index[i]]
    }
  }
  return(sample_controls_cellTypes_n_samp)
}

# NOTE: this function will need to be applied twice -> once for myeloid, once for lymphoid








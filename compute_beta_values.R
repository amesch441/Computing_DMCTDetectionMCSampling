#' Compute Beta Values 
#' 
#' This function computes the beta values using the DNA methylated and unmethylated
#' values. Specifically beta = m/(m+u+100) here.
#'
#' @param methylated matrix containing DNA methylated data
#' @param unmethylated matrix containing DNA unmethylated data
#' @return data matrix containing the computed beta values
#' @export

compute_beta_values = function(methylated, unmethylated)
{
  beta_cellTypes_denom = (methylated + unmethylated + 100) 
  
  # if there was incorrected input for methylated or unmethylated
  try(log(beta_cellTypes_denom)) # this should NEVER be 0 or negative for proper beta-values
  
  beta_cellTypes = methylated / beta_cellTypes_denom

 return(beta_cellTypes)
}

# NOTE: this function will need to be applied twice -> once for myeloid, once for lymphoid

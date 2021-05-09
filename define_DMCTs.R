
#' Define DMCT's from the DNAm profiles for the control subjects
#' 
#' This function uses the threshold of beta < 0.1 to 
#' find CpGs that were unmethylated. 
#'
#' @param type string indicating the type of DMCT, default = "non-specific", options: "non-specific", "myeloid", or "lymphoid"
#' @param beta_myeloid data frame containing computed beta values for myeloid cells
#' @param beta_lymphoid data frame containing computed beta values for lyphoid cells
#' @return returns boolean matrix for CpGs 
#' @export


define_DMCTs = function(type = "non-specific", beta_myeloid, beta_lymphoid)
{
  # if not one of these inputs, this function does not work
  if(!try({type == "non-specific" | type == "myeloid" | type == "lymphoid"}))
  {
    type = "non-specific" # make type default
    print("Warning: incorrect input for type, defaulting to non-specific DMCTs")
  }
  
  indicator_t_cells = apply(beta_lymphoid < 0.1, 1, all) 
  indicator_monocyte_cells = apply(beta_myeloid < 0.1, 1, all) 
  
  if(type == "non-specific")
  {
    indicator_both = which(indicator_monocyte_cells == T)[which(indicator_monocyte_cells == T) %in% which(indicator_t_cells == T)]
    return(indicator_both)
  }
  
  if(type == "myeloid")
  {
    indicator_monocyte_only = which(indicator_monocyte_cells == T)[!(which(indicator_monocyte_cells == T) %in% which(indicator_t_cells == T))]
    return(indicator_monocyte_only)
  }
  
  if(type == "lymphoid")
  {
    indicator_tcell_only = which(indicator_t_cells == T)[!(which(indicator_t_cells == T) %in% which(indicator_monocyte_cells == T))]
    return(indicator_tcell_only)
  }

}








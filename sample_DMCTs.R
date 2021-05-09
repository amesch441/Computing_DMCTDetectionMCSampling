

#' Sample DMCTs from CpGs
#' 
#' This function uses random sampling to obtain a list of length n_dmcts 
#' of the DMCTs from the function: define_DMCTs()
#'
#' @param n_dmcts number of DMCTs to select 
#' @param seed random seed to set, default = 215
#' @param DMCTs_indicator boolean matrix for CpGs (generated from define_DMCTs() function)
#' @return returns list of DMCTs of length n_dmcts
#' @export


sample_DMCTs = function(n_dmcts = 1000, seed = 215, DMCTs_indicator)
{
  # randomly select 1000 CpGs to be DMCTs - 3 cases
  set.seed(seed)
  
  try ({length(length(DMCTs_indicator)) > n_dmcts})

  DMCTs_both = names(DMCTs_indicator[sample(1:length(DMCTs_indicator),n_dmcts, replace = FALSE)])
  return(DMCTs_both)
}



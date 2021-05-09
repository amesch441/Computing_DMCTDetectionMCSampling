
#' Generate Plots For Non-Specific DMCTs Results
#' 
#' This function generates the tables comparing the approximate effect size in manuscript
#' measured in the paper, and the average effect size found from the monte-carlo
#' simulation
#'
#' @param n_samp number of samples to generate (will be doubled cases/controls)
#' @param ef_tcells effect size measured from simulation for tcells (lymphoid)
#' @param ef_monocytes effect size measured from simulation for monocytes (myeloid)
#' @param a vector of the first parameter to use for the beta distribution in MC sampling for generating DNAm pofile of case subjects
#' @param b vector of the second parameter to use for the beta distribution in MC sampling for generating DNAm pofile of case subjects
#' @return Table of parameter values, approximate effect size in manuscript, and average effect size for each cell type
#' @importFrom magrittr "%>%"
#' @importFrom kableExtra "kbl"
#' @importFrom kableExtra "kable_classic"
#' @importFrom kableExtra "add_header_above"
#' @export


calc_effect_size = function(a, b, ef_tcells, ef_monocytes, n_samp)
{
  # make sure n_samp is a positive number (greater than zero)
  try(log(n_samp))
  
  effectSize_tcell_mean1 = rowMeans(ef_tcells)
  effectSize_monocyte_mean1 = rowMeans(ef_monocytes)
  Estimate = c(1, 1.1, 1.4, 1.9, 2.7, 3.5, 4.4, 6.8)
  
  ef = data.frame(a, b, Estimate, Lymphoid = effectSize_tcell_mean1, Myeloid  = effectSize_monocyte_mean1)
  
  table = ef %>%
    kbl(caption = paste0("Effect Size Chart non-specific DMCTs, n = ", n_samp),  align = 'c', col.names = c("a", "b", "Estimate","Lymphoid "," Myeloid ")) %>%
    kable_classic("striped", full_width = T,  html_font = "Cambria", bootstrap_options = "striped")
  add_header_above(table, c("Beta Parameters" = 2, "Effect Size" = 3))
  
  return(table)

}
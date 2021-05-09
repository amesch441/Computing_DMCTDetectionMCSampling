
#' Runs the Simulation for MC iterations for Detecting Lymphoid and Myeloid DMCTS
#' 
#' This function runs the monote-carlo sampleing over multiple iterataions
#' using the simulation_setup() function
#'
#' @param a vector of the first parameter to use for the beta distribution in MC sampling for generating DNAm pofile of case subjects
#' @param b vector of the second parameter to use for the beta distribution in MC sampling for generating DNAm pofile of case subjects
#' @param n_samp number of samples to generate, default = 100
#' @param n_dmcts number of DMCTs to select, default = 1000
#' @param seed random seed to set, default = 215
#' @param mc number of monte-carlo runs, default = 10
#' @param DMCTs vector of DMCTs to use
#' @param sample_controls_mixed list of matrices containing mixed samples of cell-types for the control subjects
#' @param sample_controls_tcells list of matrices containing samples of only tcell for the control subjects
#' @param sample_controls_monocytes list of matrices containing samples of of only monocytes for the control subjects
#' @param lymphoid_frac vector of fractions for lymphoid cells
#' @param myeloid_frac vector of fractions for myeloid cells
#' @param spec indicator if cell-type specific
#' @param type indicator of DMCT type, default = "lymphoid"
#' @param parallelOption boolean for running simulation in parellel, default = FALSE
#' @return returns list of results including (1) average sensitivity, 
#' (2) sensitivity values for each MC samples for lymphoid cellss, 
#' (3) sensitivity values for each MC samples for myeloid cells, 
#' (4) average effect size for lymphoid cells across MC samples,
#' (5) average effect size for myeloid cells across MC samples
#' @importFrom stats "rbeta"
#' @importFrom stats "var"
#' @import foreach
#' @import doParallel
#' @importFrom foreach "%dopar%"
#' @importFrom EpiDISH "CellDMC"
#' @importFrom parallel "makeCluster"
#' @import dplyr
#' @export

simulation_run = function(a, b, parallelOption = FALSE, mc = 10, n_dmcts = 1000, type = "lymphoid", spec = FALSE, n_samp = 100, seed = 215, DMCTs,  sample_controls_mixed, lymphoid_frac, myeloid_frac,sample_controls_tcells, sample_controls_monocytes)
{

  num_cores = 1
  if(parallelOption == TRUE)
  {
    num_cores = detectCores()
  }
  
    # in order to run on windows or mac 
    clus_cores =  makeCluster(num_cores-1,type="SOCK")
    registerDoParallel(clus_cores)
    
    sim_results = foreach(i=1:mc, .combine=cbind) %dopar% {
      
      out = simulation_setup(a, b,  n_dmcts = n_dmcts, type = type, spec = spec, n_samp = n_samp, seed = seed, DMCTs = DMCTs,  sample_controls_mixed, lymphoid_frac, myeloid_frac,sample_controls_tcells, sample_controls_monocytes)
      out
      }
    
  return(sim_results)
}
  


#' Compares Parallel Simulation Times 
#' 
#' Performs a comparison of the run time for the monte-carlo simulation run times
#' using parallelOption = False vs. parallelOption = True
#'
#' @param a vector of the first parameter to use for the beta distribution in MC sampling for generating DNAm pofile of case subjects
#' @param b vector of the second parameter to use for the beta distribution in MC sampling for generating DNAm pofile of case subjects
#' @param times number of times to evaluate parallel programing
#' @param DMCTs vector of DMCTs to use
#' @param sample_controls_mixed list of matrices containing mixed samples of cell-types for the control subjects
#' @param sample_controls_tcells list of matrices containing samples of only tcell for the control subjects
#' @param sample_controls_monocytes list of matrices containing samples of of only monocytes for the control subjects
#' @param lymphoid_frac vector of fractions for lymphoid cells
#' @param myeloid_frac vector of fractions for myeloid cells
#' @return vector comaring run-times
#' @importFrom stats "rbeta"
#' @importFrom stats "var"
#' @import foreach
#' @import parallel
##' @import doParallel
#' @importFrom foreach "%dopar%"
#' @importFrom EpiDISH "CellDMC"
#' @importFrom parallel "makeCluster"
#' @import dplyr
#' @export

compare_parallel = function(a, b,times, DMCTs, sample_controls_mixed, lymphoid_frac, myeloid_frac, sample_controls_tcells, sample_controls_monocytes)
{
  time_nonparallel = rep(NA, times)
  time_parallel = rep(NA, times)
  
  for(i in 1:times)
  {
    start_time <-Sys.time()
    sim_results = simulation_run(a = a, b = b, parallelOption = FALSE, mc = 10, n_dmcts = 1000, n_samp = 100,  seed = 215, DMCTs = DMCTs, sample_controls_mixed = sample_controls_mixed, lymphoid_frac = lymphoid_frac, myeloid_frac = myeloid_frac, sample_controls_tcells = sample_controls_tcells, sample_controls_monocytes = sample_controls_monocytes)
    end_time <-Sys.time() 
    time_nonparallel[i] = start_time - end_time
  }
  
  for(i in 1:times)
  {
    start_time <-Sys.time()
    sim_results = simulation_run(a = a, b = b, parallelOption = TRUE, mc = 10, n_dmcts = 1000, n_samp = 100,  seed = 215, DMCTs = DMCTs, sample_controls_mixed = sample_controls_mixed, lymphoid_frac = lymphoid_frac, myeloid_frac = myeloid_frac, sample_controls_tcells = sample_controls_tcells, sample_controls_monocytes = sample_controls_monocytes)
    end_time <-Sys.time() 
    time_parallel[i] = start_time - end_time
  }
  
  return(cbind(time_nonparallel, time_parallel))
}



#' Sets-up the Simulation for Detecting Lymphoid and Myeloid DMCTS
#' 
#' This function sets up the monote-carlo sampleing frame work. For the case patients it samples from 
#' a Beta(a,b) distribution, where parameters a, b are input from the user. It also
#' generates a mixture of cell-type DNAm profiles for the test patients. This function 
#' calculates the effect size as noted in 'Age-related 
#' variations in the methylome associated with gene expression in human
#' monocytes and T cells'. Lastly, CellDMC() is used to determine the sensitivity
#' of this method on the simulated data.
#'
#' @param a vector of the first parameter to use for the beta distribution in MC sampling for generating DNAm pofile of case subjects
#' @param b vector of the second parameter to use for the beta distribution in MC sampling for generating DNAm pofile of case subjects
#' @param n_samp number of samples to generate, default = 100
#' @param n_dmcts number of DMCTs to select, default = 1000
#' @param seed random seed to set, default = 215
#' @param DMCTs vector of DMCTs to use
#' @param sample_controls_mixed list of matrices containing mixed samples of cell-types for the control subjects
#' @param sample_controls_tcells list of matrices containing samples of only tcell for the control subjects
#' @param sample_controls_monocytes list of matrices containing samples of of only monocytes for the control subjects
#' @param lymphoid_frac vector of fractions for lymphoid cells
#' @param myeloid_frac vector of fractions for myeloid cells
#' @param spec indicator if cell-type specific
#' @param type indicator of DMCT type, default = "lymphoid"
#' @return list of single run output list for the simulation
#' @importFrom stats "rbeta"
#' @importFrom stats "var"
#' @importFrom EpiDISH "CellDMC"
#' @export

simulation_setup = function(a, b,  n_dmcts = 1000, type = "lymphoid", spec = FALSE, n_samp = 100, seed = 215, DMCTs,  sample_controls_mixed, lymphoid_frac, myeloid_frac,sample_controls_tcells, sample_controls_monocytes)
{
  set.seed(seed)
  npar = length(a)
  mc = 1 # note: the loops here are not actually running loops (since mc = 1)
  
  sample_cases_tcells = matrix(NA,n_dmcts, n_samp)
  sample_cases_tcells = rep(list(sample_cases_tcells),mc) #contains n_dmcts x n samples for each MC sample
  
  sample_cases_monocytes = matrix(NA,n_dmcts, n_samp)
  sample_cases_monocytes = rep(list(sample_cases_monocytes),mc) #contains n_dmcts x n samples for each MC sample
  
  sample_cases_mixed = matrix(NA,n_dmcts, n_samp)
  sample_cases_mixed = rep(list(sample_cases_mixed),mc) #contains n_dmcts x n samples for each MC sample
  
  sensitivity_1 = matrix(NA, ncol = 2, nrow = mc)
  specificity_1 = matrix(NA, ncol = 2, nrow = mc)
  ppv_1 = matrix(NA, ncol = 2, nrow = mc)
  
  np = npar/2
  n_t = npar*10^mc # 'npar' parameter sets across 'mc' MC samples
  if(spec){n_t = (npar*mc)^2/(np*n_samp) - np}
  if(type == "myeloid"){n_t = (10^9)*(100/n_samp)^2}
  
  sensitivity_mean_1 = matrix(NA, ncol = 2, nrow = npar)
  colnames(sensitivity_mean_1) = c("Lymphoid", "Myeloid")
  
  sensitivity_points_lymph1 = matrix(NA, ncol = mc, nrow = npar)
  sensitivity_points_myeloid1 = matrix(NA, ncol = mc, nrow = npar)
  
  effectSize_tcell_1 = matrix(NA, ncol = mc, nrow = npar)
  effectSize_monocyte_1 = matrix(NA, ncol = mc, nrow = npar)
  
  pheno.v = c(rep(1,n_samp), rep(0,n_samp))
  frac.m = cbind(Tcell = lymphoid_frac, monocyte= myeloid_frac)
  
  lymphoid_frac_case1 = lymphoid_frac[1:n_samp]
  lymphoid_frac_control1 = lymphoid_frac[(n_samp+1):(2*n_samp)]
  
  myeloid_frac_case1 = myeloid_frac[1:n_samp]
  myeloid_frac_control1 = myeloid_frac[(n_samp+1):(2*n_samp)]
  
  for(p in 1:npar)
  {
    print(paste0("Parameter set: ", p))
    # MC sampling for the cases
    for(j in 1:mc) # MC samples
    {
      for(i in 1:n_samp)
      {
        sample_cases_tcells[[j]][,i] = rbeta(n_dmcts, a[p], b[p])
        sample_cases_monocytes[[j]][,i] = rbeta(n_dmcts, a[p], b[p])
        sample_cases_mixed[[j]][,i] = lymphoid_frac_case1[i]*sample_cases_tcells[[j]][,i] + myeloid_frac_case1[i]*sample_cases_monocytes[[j]][,i] 
      } #end for loop for MC sampling of for each sample 1...n
      
      # Compute effect size for both Myloid and Lymphoid cells 
      v1.t = mean(apply(sample_cases_tcells[[j]], 2,var))
      v2.t = mean(apply(sample_controls_tcells[[j]], 2,var))
      
      v1.m = mean(apply(sample_cases_monocytes[[j]], 2,var))
      v2.m = mean(apply(sample_controls_monocytes[[j]], 2,var))
      
      effectSize_tcell_1[p,j] = mean(abs(colMeans(sample_cases_tcells[[j]]) - colMeans(sample_controls_tcells[[j]]))) / sqrt(0.5 *(v1.t+v2.t))
      effectSize_monocyte_1[p,j] = mean(abs(colMeans(sample_cases_monocytes[[j]]) - colMeans(sample_controls_monocytes[[j]]))) / sqrt(0.5 *(v1.m+v2.m))
      
    } #end for loop for simulating case data using MC sampling (across all MC samples) 1...MC
    
    for(t in 1:mc)
    {
      beta.m = cbind(sample_cases_mixed[[t]], sample_controls_mixed[[t]])
      rownames(beta.m) = DMCTs
      
      # RUN CellDMC method
      celldmc.o = CellDMC(beta.m = beta.m, pheno.v = pheno.v, frac.m = frac.m, mc.cores = 1, adjPMethod = "BH",   adjPThresh = 0.05/n_t) # runs in parallel
      
      # Compute summary statistics: sensitivity
      TP_t = length(which(celldmc.o$dmct[,2] == 0))
      TP_mono = length(which(celldmc.o$dmct[,3] == 0))
      
      sensitivity_1[t,1] = 1 - TP_t / n_dmcts
      sensitivity_1[t,2] = 1 - TP_mono / n_dmcts
      
    } #end for loop for running CellDMC() on each MC sample
    sensitivity_points_lymph1[p,] = sensitivity_1[,1]
    sensitivity_points_myeloid1[p,] = sensitivity_1[,2]
    sensitivity_mean_1[p,] = colMeans(sensitivity_1)
  } #end for loop across p parameter sets
  
  
  return(list("sensitivity_mean" = sensitivity_mean_1, "sensitivity_points_lymph" = sensitivity_points_lymph1, "sensitivity_points_myeloid" = sensitivity_points_myeloid1,  "effect_size_tcell" = effectSize_tcell_1, "effect_size_monocyte" = effectSize_monocyte_1))
  
} # end run_simulation function



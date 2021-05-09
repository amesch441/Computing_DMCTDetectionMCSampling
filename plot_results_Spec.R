#' Generate Plots For Cell-type-Specific DMCTs Results
#' 
#' This function generates the scatter plot (Effect size vs Sensitivity) for the 
#' cell-type-specific DMCTS. 
#'
#' @param sim_results_Spec run_simulation() output using myeloid-specific DMCTs for n = 200 (100 cases, 100 controls) 
#' @param sim_results_Spec_600 run_simulation() output using myeloid-specific DMCTs for n = 600 (300 cases, 300 controls)
#' @param sim_results_Spec_100 run_simulation() output using myeloid-specific DMCTs for n = 100 (50 cases, 50 controls)
#' @param sim_results_Spec_800 run_simulation() output using myeloid-specific DMCTs for n = 800 (400 cases, 400 controls)
#' @param cell_type cell type name being plotted, default = "myeloid" (options: "myeloid", "lymphoid")
#' @return Scatter plot for Effect size vs Sensitivity
#' @importFrom graphics "lines"
#' @importFrom graphics "legend"
#' @export

plot_results_Spec = function(cell_type = "myeloid", sim_results_Spec, sim_results_Spec_600, sim_results_Spec_100, sim_results_Spec_800)
{
  
  if(cell_type == "myeloid")
  {
    sensitivity_mean_200 = sim_results_Spec$sensitivity_mean[,2]
    sensitivity_mean_600 = sim_results_Spec_600$sensitivity_mean[,2]
    sensitivity_mean_100 = sim_results_Spec_100$sensitivity_mean[,2]
    sensitivity_mean_800 = sim_results_Spec_800$sensitivity_mean[,2]
    effect_size_est = c(1, 1.1, 1.4, 1.9, 2.7, 3.5, 4.4, 6.8)
    
  }
  
  if(cell_type == "lymphoid")
  {
    sensitivity_mean_200 = sim_results_Spec$sensitivity_mean[,1]
    sensitivity_mean_600 = sim_results_Spec_600$sensitivity_mean[,1]
    sensitivity_mean_100 = sim_results_Spec_100$sensitivity_mean[,1]
    sensitivity_mean_800 = sim_results_Spec_800$sensitivity_mean[,1]
    effect_size_est = c(1, 1.1, 1.4, 1.9, 2.7, 3.5, 4.4, 6.8)
  }

  plot = plot(effect_size_est, sensitivity_mean_200, type = "b", col = "brown4", lwd = 3, pch = 23, xlab = "Effect size", ylab = "SE(DMCT)", main = paste0(cell_type, "-specific DMCTs"))
  lines(effect_size_est, sensitivity_mean_600, type = "b", col = "purple", lwd = 3, pch = 23)
  lines(effect_size_est, sensitivity_mean_100, type = "b", col = "orange", lwd = 3, pch = 23)
  lines(effect_size_est, sensitivity_mean_800, type = "b", col = "thistle4", lwd = 3, pch = 23)
  
  legend("bottomright", c("n = 100","n = 200","n = 600","n = 800"), cex = 0.85, col = c("orange","brown4", "purple", "thistle4"), lwd = 3)
  
  return(plot)
}
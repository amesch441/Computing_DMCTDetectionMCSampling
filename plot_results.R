#' Generate Plots For Non-Specific DMCTs Results
#' 
#' This function generates the scatter plot (Effect size vs Sensitivity) for the 
#' non-specific DMCTS. 
#'
#' @param sim_results run_simulation() output for n = 200 (100 cases, 100 controls) 
#' @param sim_results_600 run_simulation() output for n = 600 (300 cases, 300 controls)
#' @return Scatter plot for Effect size vs Sensitivity
#' @importFrom graphics "lines"
#' @importFrom graphics "legend"
#' @importFrom graphics "points"
#' @export



plot_results = function(sim_results, sim_results_600)
{
  
  sensitivity_mean_1 = sim_results[[1]]
  sensitivity_points_lymph1 = sim_results[[2]]
  sensitivity_points_myeloid1 = sim_results[[3]]
  
  sensitivity_mean_2 = sim_results_600[[1]]
  sensitivity_points_lymph2 = sim_results_600[[2]]
  sensitivity_points_myeloid2 = sim_results_600[[3]]
  
  effect_size_t = rowMeans(sim_results$effect_size_tcell)
  effect_size_m =  rowMeans(sim_results$effect_size_monocyte)
  
  effect_size_t_600 = rowMeans(sim_results_600$effect_size_tcell)
  effect_size_m_600 =  rowMeans(sim_results_600$effect_size_monocyte)
  
  plot = plot(effect_size_t, sensitivity_mean_1[,1], type = "b", col = "blue", lwd = 3, pch = 23, xlab = "Effect size", ylab = "SE(DMCT)", main = "non-specific DMCTs")
  lines(effect_size_m, sensitivity_mean_1[,2], type = "b", col = "brown4", lwd = 3, pch = 23)
  lines(effect_size_t_600, sensitivity_mean_2[,1], type = "b", col = "darkgreen", lwd = 3, pch = 23)
  lines(effect_size_m_600, sensitivity_mean_2[,2], type = "b", col = "purple", lwd = 3, pch = 23)
  
  for(i in 1:8)
  {
    points(rep(effect_size_t[i],10), sensitivity_points_lymph1[i,], col = "lightblue", pch = 18, cex = 1.5)
    points(rep(effect_size_m[i],10), sensitivity_points_myeloid1[i,], col = "red", pch = 18, cex = 1.5)
    
    points(rep(effect_size_t_600[i],10), sensitivity_points_lymph2[i,], col = "green", pch = 18, cex = 1.5)
    points(rep(effect_size_m_600[i],10), sensitivity_points_myeloid2[i,], col = "pink", pch = 18, cex = 1.5)
  }
  
  legend("bottomright", c("lym(n = 200)","mye(n = 200)","lym(n = 600)","mye(n = 600)"), cex = 0.85, col = c("blue","brown4", "darkgreen", "purple"), lwd = 3)
  
  return(plot)
}
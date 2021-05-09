
#' Import DNAm profiles 
#' 
#' This function can be used to import the real DNAm profiles from 'Age-related 
#' variations in the methylome associated with gene expression in human
#' monocytes and T cells'. The .rds data files are uploaded with the package.
#' The four data files include methylated and unmethylated files for both lymphoid
#' and myeloid cells
#'
#' @param file_name file name
#' @return data matrix containing the DNAm profile (methylated or unmethylated)
#' @export

import_real_DNAm_profiles = function(file_name)
{
  file = load(file_name)
  return(get(file))
}


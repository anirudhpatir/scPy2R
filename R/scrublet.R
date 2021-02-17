#' @title Scrublet
#' @description See preprint: Scrublet: computational identification of cell doublets in single-cell transcriptomic data
#' Samuel L Wolock, Romain Lopez, Allon M Klein.  bioRxiv 357368; doi: https://doi.org/10.1101/357368
#' @param data expression matrix
#' @param python_home The python home directory where Scrublet is installed
#' @param return_results_only bool (optional, default False)
#' @param min_counts, int (optional, default=2), See scrublet reference
#' @param min_cells, int (optional, default=3), See scrublet reference
#' @param expected_doublet_rate, float (optional, default=0.1), See scrublet reference - expected_doublet_rate: the
#' fraction of transcriptomes that are doublets, typically 0.05-0.1. Results are not particularly sensitive to this parameter. For this example, the expected doublet rate comes from the Chromium User Guide: https://support.10xgenomics.com/permalink/3vzDu3zQjY0o2AqkkkI4CC
#' @param min_gene_variability_pctl, int (optional, default=85), See scrublet reference
#' @param n_prin_comps, int (optional, default=30), See scrublet reference  (Number of principal components to use)
#' @param sim_doublet_ratio, int (optional, default=2),  the number of doublets to simulate, relative to the number of observed transcriptomes. This should be high enough that all doublet states are well-represented by simulated doublets. Setting it too high is computationally expensive. The default value is 2, though values as low as 0.5 give very similar results for the datasets that have been tested.
#' @param n_neighbors, int (optional) n_neighbors: Number of neighbors used to construct the KNN classifier of observed transcriptomes and simulated doublets. The default value of round(0.5*sqrt(n_cells)) generally works well.
#' Return only a list containing scrublet output
#' @return The doublet_score output from scrublet,
#' @importFrom reticulate use_python
#' @importFrom reticulate source_python
#' @export
scrublet_R <- function(
  data,
  python_home = system("which python", intern = TRUE),
  return_results_only = FALSE,
  min_counts=2,
  min_cells=3,
  expected_doublet_rate=0.1,
  min_gene_variability_pctl=85,
  n_prin_comps=30,
  sim_doublet_ratio=2,
  n_neighbors=NULL
  )
{
  reticulate::use_python(python_home)

  #####################################
  # Checks
  #####################################

  if(!py_available("scrublet")) stop("python module scrublet does not seem to be installed; - try running 'py_config()'")

  #####################################
  # Run
  #####################################

  reticulate::source_python(paste(system.file(package = "scPy2R"),
                                  "scrublet.py", sep = "/"))

  # Parameters
  X <- as(t(data), "TsparseMatrix")
  i <- as.integer(X@i)
  j <- as.integer(X@j)
  val <- X@x
  dim <- as.integer(X@Dim)
  if(is.null(n_neighbors)){
    n_neighbors<-round(0.5*sqrt(nrow(X)))
  }
  scrublet_py_args<-c(list(i=i, j=j, val=val, dim=dim,
                           expected_doublet_rate=expected_doublet_rate, min_counts=min_counts,
                           min_cells=min_cells, min_gene_variability_pctl=min_gene_variability_pctl, n_prin_comps=n_prin_comps,
                           sim_doublet_ratio=sim_doublet_ratio, n_neighbors=n_neighbors))
  # Call scrublet
  scrublet_res <- do.call(scrublet_py, scrublet_py_args)

  #####################################
  # Save
  #####################################

  names(scrublet_res)<-c("doublet_scores_sim", "threshold","doublet_scores", "predicted_doublets")
  return(scrublet_res)
}

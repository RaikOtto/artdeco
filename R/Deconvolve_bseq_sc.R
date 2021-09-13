#' Deconvolve_bseq_sc
#' 
#' \code{Deconvolve_bseq_sc} Utilizes BSeq-sc for deconvolution. 
#' Wrapper function.
#' 
#' @param deconvolution_data Data to be deconvolved.
#' @param models_list List of models.
#' @param models Models utilized.
#' @param Cibersort_absolute_mode CIBRSORT absolute mode
#' @param nr_permutations Amount of permutations.
#' @return Matrix containing the deconvolution results.
Deconvolve_bseq_sc = function(
    deconvolution_data,
    models_list,
    models,
    Cibersort_absolute_mode,
    nr_permutations
){
    
    prediction_res_coeff_list = list()
    prediction_stats_list = list()
    
    for (pred_model in models){
        
        print(paste0("Deconvolving with model: ",pred_model))
        model_basis = models_list[[pred_model]]
        model_basis = model_basis[[1]]
        
        bseqsc_fit = bseqsc::bseqsc_proportions(
            deconvolution_data,
            model_basis,
            verbose = FALSE,
            absolute = Cibersort_absolute_mode,
            log = FALSE,
            perm = nr_permutations
        )
        
        prediction_res_coeff_list[[pred_model]] = bseqsc_fit$coefficients
        prediction_stats_list[[pred_model]] = bseqsc_fit$stats
    }
    
    # create results matrix called meta_data
    
    deconvolution_results = prepare_result_matrix_bseqsc(
        prediction_res_coeff_list = prediction_res_coeff_list,
        deconvolution_data = deconvolution_data,
        prediction_stats_list = prediction_stats_list,
        models = models
    )
    
    subtype_cands = c()
    for (model in models){
        cands = rownames(as.data.frame(prediction_res_coeff_list[model]))
        subtype_cands = c(subtype_cands, cands)
    }
    subtype_cands = unique(subtype_cands)
    
    col_idx <- match(c("model", subtype_cands, "strength_subtype", "subtype", 
                       "sig_score", "p_value", "correlation", "rmse"), 
                       str_to_lower(colnames(deconvolution_results)))
    deconvolution_results <- deconvolution_results[,col_idx]
    
    return(deconvolution_results)
}

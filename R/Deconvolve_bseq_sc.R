Deconvolve_bseq_sc = function(
    deconvolution_data,
    models_list,
    models,
    nr_permutations
){
    
    prediction_res_coeff_list = list()
    prediction_stats_list = list()
    
    for (pred_model in models){
        
        print(paste0("Deconvolving with model: ",pred_model))
        model_basis = models_list[[pred_model]]
        model_basis = model_basis[[1]]
        
        fit = bseqsc_proportions(
            deconvolution_data,
            model_basis,
            verbose = FALSE,
            absolute = TRUE,
            log = FALSE,
            perm = nr_permutations
        )
        
        prediction_res_coeff_list[[pred_model]] = fit$coefficients
        prediction_stats_list[[pred_model]] = fit$stats
    }
    
    # create results matrix called meta_data
    
    deconvolution_results = prepare_result_matrix(
        prediction_res_coeff_list = prediction_res_coeff_list,
        deconvolution_data = deconvolution_data,
        models = models
    )
    colnames(deconvolution_results) = str_to_lower(colnames(deconvolution_results))
    
    deconvolution_results = prepare_sample_result_matrix(
        deconvolution_results = deconvolution_results,
        prediction_stats_list = prediction_stats_list,
        models_list = models_list,
        transcriptome_file = transcriptome_file
    )
    return(deconvolution_results)
}

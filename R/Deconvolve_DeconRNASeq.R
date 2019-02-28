Deconvolve_DeconRNASeq = function(
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
        training_mat = exprs(model_basis[[1]])
        col_names = colnames(training_mat)
        row_names = rownames(training_mat)
        signature_mat = as.double(as.character(unlist(training_mat)))
        signature_mat = matrix(signature_mat,ncol = length(col_names))
        colnames(signature_mat) = col_names
        rownames(signature_mat) = row_names
        signature_mat = as.data.frame(signature_mat,stringsAsFactors=FALSE)

        res = DeconRNASeq(
            datasets = as.data.frame(exprs(deconvolution_data)),
            signatures = signature_mat,
            proportions = NULL,
            checksig = TRUE,
            known.prop = FALSE,
            use.scale = TRUE,
            fig = FALSE
        )
        
        deco_res = res$out.all
        rownames(deco_res) = colnames(deconvolution_data)
        
        deco_res = as.data.frame(deco_res)
        #deco_res = apply(deco_res, MARGIN = 1, FUN = function(vec){return(log(vec+1))})
        annotation_data = t(deco_res)
        annotation_data[1:5,]
        old_colnames = colnames(annotation_data)

        prediction_res_coeff_list[[pred_model]] = deco_res
        prediction_stats_list[[pred_model]]     = as.data.frame(res$out.pca)
    }
    
    # create results matrix called meta_data
    
    deconvolution_results = prepare_result_matrix_deconrnaseq(
        prediction_res_coeff_list = prediction_res_coeff_list,
        deconvolution_data = deconvolution_data,
        models = models
    )
    colnames(deconvolution_results) = str_to_lower(colnames(deconvolution_results))
    
    deconvolution_results = prepare_sample_result_matrix_deconrnaseq(
        deconvolution_results = deconvolution_results,
        prediction_stats_list = prediction_stats_list,
        deconvolution_data = deconvolution_data,
        models_list = models_list
    )
    return(deconvolution_results)
}

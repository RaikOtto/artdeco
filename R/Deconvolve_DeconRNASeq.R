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
        model_basis = model_basis[[1]]
        pData(model_basis)$sampleID = colnames(exprs(model_basis))
        subtypes = as.character(unique(pData(model_basis)$cellType))
        
        res = DeconRNASeq(
            deconvolution_data,
            signatures,
            proportions = NULL,
            checksig = FALSE,
            known.prop = FALSE,
            use.scale = TRUE,
            fig = FALSE
        )
        
        deco_res = res$out.all
        rownames(deco_res) = colnames(query_data)
        
        deco_res = as.data.frame(deco_res)
        deco_res = apply(deco_res, MARGIN = 1, FUN = function(vec){return(log(vec+1))})
        annotation_data = t(deco_res)
        annotation_data[1:5,]
        old_colnames = colnames(annotation_data)

        prediction_res_coeff_list[[pred_model]] = Est.prop.GSE50244$Est.prop.allgene
        prediction_stats_list[[pred_model]]     = Est.prop.GSE50244$r.squared.full
    }
    
    # create results matrix called meta_data
    
    deconvolution_results = prepare_result_matrix_music(
        prediction_res_coeff_list = prediction_res_coeff_list,
        deconvolution_data = deconvolution_data,
        models = models
    )
    colnames(deconvolution_results) = str_to_lower(colnames(deconvolution_results))
    
    deconvolution_results = prepare_sample_result_matrix_music(
        deconvolution_results = deconvolution_results,
        prediction_stats_list = prediction_stats_list,
        deconvolution_data = deconvolution_data,
        models_list = models_list
    )
    return(deconvolution_results)
}

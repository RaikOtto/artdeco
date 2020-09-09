Deconvolve_music = function(
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
        
        Est.prop.GSE50244 = MuSiC::music_prop(
            bulk.eset = deconvolution_data,
            sc.eset = model_basis,
            clusters = 'cellType',
            samples = 'sampleID',
            select.ct = subtypes,
            verbose = FALSE,
            iter.max = nr_permutations
        )

        prediction_res_coeff_list[[pred_model]] = Est.prop.GSE50244$Est.prop.allgene # cell type proportions
    }
    
    # create results matrix called meta_data
    
    deconvolution_results = prepare_result_matrix_music(
        prediction_res_coeff_list = prediction_res_coeff_list,
        deconvolution_data = deconvolution_data,
        models = models
    )
    colnames(deconvolution_results) = str_to_lower(colnames(deconvolution_results))
    
    return(deconvolution_results)
}

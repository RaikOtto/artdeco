prepare_result_matrix_NMF = function(
    prediction_res_coeff_list,
    deconvolution_data,
    prediction_stats_list,
    models
){
    
    rounding_precision = 1
    subtype_cands = c("alpha","beta","gamma","delta","acinar","ductal","hisc")
    nmf_parameter = c("Score")
    
    result_matrix_template = matrix( rep(0.0, length(subtype_cands) * ncol(deconvolution_data)),ncol = length(subtype_cands) )
    colnames(result_matrix_template) = subtype_cands
    result_matrix_template = as.data.frame(result_matrix_template)
    result_matrix_template$Sample_ID = colnames(deconvolution_data)
    result_matrix_template_ori = result_matrix_template
    
    for ( i in 1:( length(models) - 1)){
        result_matrix_template = rbind(result_matrix_template, result_matrix_template_ori)
    }
    
    model_vec = rep("", nrow(result_matrix_template))
    for(i in 1:length(models)){
        start_index =  (i-1) * nrow(result_matrix_template_ori)  + 1
        end_index   =  (i) * nrow(result_matrix_template_ori)
        model_vec[ start_index : end_index ] = models[i]
    }
    result_matrix_template$Model = model_vec
    
    for( i in 1:length(nmf_parameter)){
        result_matrix_template[,nmf_parameter[i]] = rep("",nrow(result_matrix_template))
    }
    
    for ( model in models ){
        
        res_coeff = prediction_res_coeff_list[[model]]
        colnames(res_coeff) = str_replace_all(colnames(res_coeff) ,"^X","")
        res_coeff[ is.na(res_coeff) ] = 0.0
        
        res_coeff = t(res_coeff)
        colnames(res_coeff) = str_to_lower(colnames(res_coeff))
        
        if ("alpha" %in% rownames(res_coeff)) # sanity check
            res_coeff = t(res_coeff)
        
        subtype_cands_found = subtype_cands[
            subtype_cands %in% colnames(res_coeff)
            ]
        
        result_matrix_template[
            result_matrix_template$Model == model,
            subtype_cands_found
            ] = res_coeff[,subtype_cands_found]
        
        prediction_stats = as.data.frame(prediction_stats_list[model])
        result_matrix_template[
            result_matrix_template$Model == model,
            nmf_parameter
        ] = prediction_stats[,1]
        
        #for (subtype in subtype_cands_found){
        
        #    result_matrix_template[ , subtype] = round((as.double(res_coeff[,subtype]) / row_sums)*100,1)
        #    result_matrix_template[result_matrix_template[ , subtype] > 100,subtype] = 100
        #}
    }
    
    colnames(result_matrix_template)[colnames(result_matrix_template) == "progenitor"] = "Progenitor"
    colnames(result_matrix_template)[colnames(result_matrix_template) == "hisc"] = "HISC"
    colnames(result_matrix_template)[colnames(result_matrix_template) == "alpha"] = "Alpha"
    colnames(result_matrix_template)[colnames(result_matrix_template) == "beta"] = "Beta"
    colnames(result_matrix_template)[colnames(result_matrix_template) == "gamma"] = "Gamma"
    colnames(result_matrix_template)[colnames(result_matrix_template) == "delta"] = "Delta"
    colnames(result_matrix_template)[colnames(result_matrix_template) == "acinar"] = "Acinar"
    colnames(result_matrix_template)[colnames(result_matrix_template) == "ductal"] = "Ductal"
    
    return(result_matrix_template)
}

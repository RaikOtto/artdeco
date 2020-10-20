prepare_result_matrix_NMF = function(
    prediction_res_coeff_list,
    deconvolution_data,
    prediction_stats_list,
    models
){
    
    rounding_precision = 1
    subtype_cands = c("alpha","beta","gamma","delta","acinar","ductal","hisc")
    nmf_parameter = c("score")
    
    result_matrix_template = matrix( rep(0.0, length(subtype_cands) * 
                                             ncol(deconvolution_data)),
                                     ncol = length(subtype_cands) )
    colnames(result_matrix_template) = subtype_cands
    result_matrix_template = as.data.frame(result_matrix_template)
    rownames(result_matrix_template) = colnames(deconvolution_data)
    result_matrix_template_ori = result_matrix_template
    
    if (length(models) > 1){
        for ( i in 1:( length(models) - 1)){
            result_matrix_template = rbind(result_matrix_template, 
                                           result_matrix_template_ori)
        }
    }
        
    model_vec = rep("", nrow(result_matrix_template)) # create "model" column
    for(i in 1:length(models)){
        start_index =  (i-1) * nrow(result_matrix_template_ori)  + 1
        end_index   =  (i) * nrow(result_matrix_template_ori)
        model_vec[ start_index : end_index ] = models[i]
    }

    result_matrix_template$model = model_vec 
    # move last column to first column
    result_matrix_template <- result_matrix_template[,c(ncol(
        result_matrix_template),1:(ncol(result_matrix_template)-1))]
    
    for( i in 1:length(nmf_parameter)){ 
        result_matrix_template[,nmf_parameter[i]] = rep(
            "", nrow(result_matrix_template))
    }
    
    for ( model in models ){
        
        res_coeff = prediction_res_coeff_list[[model]]
        colnames(res_coeff) = str_replace_all(colnames(res_coeff), "^X", "")
        res_coeff[ is.na(res_coeff) ] = 0.0
        
        res_coeff = t(res_coeff)
        colnames(res_coeff) = str_to_lower(colnames(res_coeff))
        
        if ("alpha" %in% rownames(res_coeff)) # sanity check
            res_coeff = t(res_coeff)
        
        subtype_cands_found = subtype_cands[
            subtype_cands %in% colnames(res_coeff)
            ]
        
        result_matrix_template[
            result_matrix_template$model == model,
            subtype_cands_found
            ] = res_coeff[,subtype_cands_found]
        
        prediction_stats = as.data.frame(prediction_stats_list[model])
        result_matrix_template[
            result_matrix_template$model == model,
            nmf_parameter
        ] = prediction_stats[,1]
        
        #for (subtype in subtype_cands_found){
        
        #    result_matrix_template[ , subtype] = round(
        #       (as.double(res_coeff[,subtype]) / row_sums)*100,1)
        #    result_matrix_template[result_matrix_template[ , subtype] > 
        #                              100,subtype] = 100
        #}
    }

    
    # add code from prepare_sample_result_matrix_NMF
    # so that columns subtype and strength_subtype are added
    # can be added; doesn't have to be added
    
    result_matrix_template[,"Strength_subtype"] = rep(
        "",nrow(result_matrix_template))
    result_matrix_template[,"Subtype"] = rep("",nrow(result_matrix_template))

    cands_dif_1 = c("alpha","beta","gamma","delta","acinar","ductal")
    if("hisc" %in% colnames(result_matrix_template)){
        cands_dif_2 = "hisc"
    } else {
        cands_dif_2 = c("acinar","ductal")
    }
    
    cands_dif = cands_dif_1[
        (cands_dif_1 %in% colnames(result_matrix_template)) &
            !(cands_dif_1 %in% cands_dif_2)
        ]
    
    
    for( j in 1:nrow(result_matrix_template)){
        
        max_subtype = colnames(result_matrix_template[cands_dif])[
            which.max(result_matrix_template[j,cands_dif])
            ]
        subtype_strength = result_matrix_template[j,max_subtype] /
            sum(result_matrix_template[j,cands_dif])
        if (sum(result_matrix_template[j,cands_dif]) == 0)
            subtype_strength = result_matrix_template[j,max_subtype]
        subtype_strength = round(subtype_strength * 100,rounding_precision)
        
        result_matrix_template[j,"Strength_subtype"] =
            subtype_strength
        
        if (subtype_strength == 0)
            max_subtype = "not_significant"
        
        result_matrix_template[j,"Subtype"] = max_subtype
    }
    
    result_matrix_template$Strength_subtype <- as.numeric(
        result_matrix_template$Strength_subtype)
    
    return(result_matrix_template)
}

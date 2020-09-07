prepare_sample_result_matrix_music = function(
    deconvolution_results,
    prediction_stats_list,
    deconvolution_data,
    models_list
){
    
    rounding_precision = 1

    deconvolution_results[,"Strength_subtype"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Subtype"] = rep("",nrow(deconvolution_results))
    
    ###
    
    res_cor = prediction_stats_list[[1]]
    res_cor[ is.na(res_cor) ] = 0.0
    
    cands_dif_1 = c("alpha","beta","gamma","delta","acinar","ductal")
    if("hisc" %in% colnames(deconvolution_results)){
        cands_dif_2 = "hisc"
    } else {
        cands_dif_2 = c("acinar","ductal")
    }
    
    deconvolution_results[,"Confidence_score_dif"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Confidence_score_dif"] = round(abs(as.double(res_cor)),rounding_precision)
    
    cands_dif = cands_dif_1[
        (cands_dif_1 %in% colnames(deconvolution_results)) &
        !(cands_dif_1 %in% cands_dif_2)
    ]
    
    for( j in 1:ncol(deconvolution_data)){
        
        max_subtype = colnames(deconvolution_results[cands_dif])[
            which.max(deconvolution_results[j,cands_dif])
            ]
        subtype_strength = deconvolution_results[j,max_subtype] /
            sum(deconvolution_results[j,cands_dif])
        if (sum(deconvolution_results[j,cands_dif]) == 0)
            subtype_strength = deconvolution_results[j,max_subtype]
        subtype_strength = round(subtype_strength * 100,rounding_precision)
        
        deconvolution_results[j,"Strength_subtype"] =
            subtype_strength
        
        if (subtype_strength == 0)
            max_subtype = "not_significant"
        
        deconvolution_results[j,"Subtype"] = max_subtype
    }
    
    deconvolution_results[,"Differentiation_score"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Differentiation_score"] = round(res_cor,rounding_precision)
    
    ### de dif
    
    res_cor = prediction_stats_list[[2]]
    res_cor[ is.na(res_cor) ] = 0.0
    
    deconvolution_results[,"Strength_de_differentiation"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Confidence_score_de_dif"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Confidence_score_de_dif"] = round(abs( as.double(res_cor) ),rounding_precision)
    
    for( j in 1:ncol(deconvolution_data)){
        
        de_strength = log( sum(deconvolution_results[j,cands_dif_2]) /
                               sum(deconvolution_results[j,cands_dif_2]) + 1 )
        if (sum(deconvolution_results[j,cands_dif_2]) == 0)
            de_strength = 0
        deconvolution_results[j,"Strength_de_differentiation"] = 
            round(as.double(de_strength),rounding_precision)
        
    }
    
    return(deconvolution_results)
}

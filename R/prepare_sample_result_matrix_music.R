
prepare_sample_result_matrix_music = function(
    deconvolution_results,
    prediction_stats_list,
    deconvolution_data,
    models_list
){
    
    #models = as.character(unlist(str_split(deconvolution_results$Model[1],pattern = "\\|")))
    deconvolution_results[,"R_sqr_dif"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Strength_subtype"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Subtype"] = rep("",nrow(deconvolution_results))

    ###
    
    res_cor = prediction_stats_list[[1]]
    res_cor[ is.na(res_cor) ] = 0.0

    cands_dif = c("alpha","beta","gamma","delta","acinar","ductal")
    cands_dif = cands_dif[cands_dif %in% colnames(deconvolution_results)]
    deconvolution_results[,"Confidence_score_dif"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Confidence_score_dif"] = abs(as.double(res_cor))
    #deconvolution_results[,"Confidence_score_dif"] = deconvolution_results[,"Confidence_score_dif"] /
    #    max(deconvolution_results[,"Confidence_score_dif"])
    
    for( j in 1:ncol(deconvolution_data)){
        
        max_subtype = colnames(deconvolution_results[cands_dif])[
            which.max(deconvolution_results[j,cands_dif])
        ]
        subtype_strength = deconvolution_results[j,max_subtype] /
            sum(deconvolution_results[j,cands_dif])
        if (sum(deconvolution_results[j,cands_dif]) == 0)
            subtype_strength = deconvolution_results[j,max_subtype]
        subtype_strength = round(subtype_strength * 100,1)
        
        deconvolution_results[j,"Strength_subtype"] =
            subtype_strength
        
        if (subtype_strength == 0)
            max_subtype = "not_significant"
        
        deconvolution_results[j,"Subtype"] = max_subtype
    }
    
    deconvolution_results[,"Differentiation_score"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Differentiation_score"] = res_cor

    ### de dif
    
    res_cor = prediction_stats_list[[2]]
    res_cor[ is.na(res_cor) ] = 0.0

    cands_de_dif = c("progenitor","hisc","hesc")
    deconvolution_results[,"Strength_de_differentiation"] = rep("",nrow(deconvolution_results))
    cands_de_dif = cands_de_dif[cands_de_dif %in% colnames(deconvolution_results)]
    deconvolution_results[,"Confidence_score_de_dif"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Confidence_score_de_dif"] = abs( as.double(res_cor) )
    #deconvolution_results[,"Confidence_score_de_dif"] = deconvolution_results[,"Confidence_score_de_dif"] / 
    #    max(deconvolution_results[,"Confidence_score_de_dif"])
    
    for( j in 1:ncol(deconvolution_data)){

        de_strength = log( sum(deconvolution_results[j,cands_de_dif]) /
            sum(deconvolution_results[j,cands_dif]) +1 )
        if (sum(deconvolution_results[j,cands_dif]) == 0)
            de_strength = 0
        deconvolution_results[j,"Strength_de_differentiation"] = 
            as.double(de_strength)

    }

    return(deconvolution_results)
}

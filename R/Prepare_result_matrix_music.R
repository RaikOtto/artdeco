prepare_result_matrix_music = function(
    prediction_res_coeff_list,
    deconvolution_data,
    models
){
    
    rounding_precision = 1
    result_matrix = data.frame(
        row.names = colnames(deconvolution_data),
        "model" = rep( paste0(c(models),collapse="|"), ncol(deconvolution_data))
    )
    
    subtype_cands = c("alpha","beta","gamma","delta","acinar","ductal","hisc")
    
    for (nr_fit in 1:length(models)){
        
        if (nr_fit == 1){
            
            res_coeff_1 = prediction_res_coeff_list[[nr_fit]]
            colnames(res_coeff_1) = str_replace_all(colnames(res_coeff_1), 
                                                    "^X", "")
            
            res_coeff_1[ is.na(res_coeff_1) ] = 0.0
            
            # important sanity check
            if ("alpha" %in% str_to_lower(rownames(res_coeff_1))) 
                res_coeff_1 = t(res_coeff_1)
            
            colnames(res_coeff_1) = str_to_lower(colnames(res_coeff_1))
            subtype_cands_1 = subtype_cands[subtype_cands %in% 
                                                colnames(res_coeff_1)]
            
            row_sums = rowSums(res_coeff_1[,subtype_cands[
                subtype_cands %in% colnames(res_coeff_1)]])
            
            for (subtype in subtype_cands_1){
                
                if (!(subtype %in% colnames(res_coeff_1)) ) next
                
                result_matrix[ , subtype] = round(
                    (res_coeff_1[,subtype] / row_sums)*100,rounding_precision)
                result_matrix[result_matrix[ , subtype] > 100,subtype] = 100
            }
        }
        ## end fit 1
        
        if ( nr_fit == 2 ){
            
            res_coeff_2 = prediction_res_coeff_list[[nr_fit]]
            colnames(res_coeff_2) = str_replace_all(colnames(res_coeff_2), 
                                                    "^X", "")
            res_coeff_2[ is.na(res_coeff_2) ] = 0.0
            
            res_coeff_2 = t(res_coeff_2)
            colnames(res_coeff_2) = str_to_lower(colnames(res_coeff_2))
            
            if ("alpha" %in% rownames(res_coeff_2)) # sanity check
                res_coeff_2 = t(res_coeff_2)
            
            subtype_cands_2 = subtype_cands[
                subtype_cands %in% colnames(res_coeff_2)
                ]
            subtype_cands_2 = subtype_cands_2[!(subtype_cands_2 %in% 
                                                    subtype_cands_1) ]
            row_sums = rowSums(res_coeff_2[,subtype_cands[
                subtype_cands %in% colnames(res_coeff_2)]])
            
            for (subtype in subtype_cands_2){
                
                result_matrix[ , subtype] = round(
                    (as.double(res_coeff_2[,subtype]) / row_sums)*100,1)
                result_matrix[result_matrix[ , subtype] > 100,subtype] = 100
            }
        }
        ## end fit 2
    }
    

    # add code from prepare_sample_result_matrix_NMF
    # so that columns subtype and strength_subtype are added
    # can be added; doesn't have to be added
    
    result_matrix[,"Strength_subtype"] = rep("",nrow(result_matrix))
    result_matrix[,"Subtype"] = rep("",nrow(result_matrix))
    
    cands_dif_1 = c("alpha","beta","gamma","delta","acinar","ductal")
    if("hisc" %in% colnames(result_matrix)){
        cands_dif_2 = "hisc"
    } else {
        cands_dif_2 = c("acinar","ductal")
    }
    
    cands_dif = cands_dif_1[
        (cands_dif_1 %in% colnames(result_matrix)) &
            !(cands_dif_1 %in% cands_dif_2)
        ]
    
    for( j in 1:nrow(result_matrix)){
        
        max_subtype = colnames(result_matrix[cands_dif])[
            which.max(result_matrix[j,cands_dif])
            ]
        subtype_strength = result_matrix[j,max_subtype] /
            sum(result_matrix[j,cands_dif])
        if (sum(result_matrix[j,cands_dif]) == 0)
            subtype_strength = result_matrix[j,max_subtype]
        subtype_strength = round(subtype_strength * 100,rounding_precision)
        
        result_matrix[j,"Strength_subtype"] = subtype_strength
        
        if (subtype_strength == 0)
            max_subtype = "not_significant"
        
        result_matrix[j,"Subtype"] = max_subtype
    }
    
    result_matrix$Strength_subtype <- as.numeric(
        result_matrix$Strength_subtype)
    
    return(result_matrix)
}

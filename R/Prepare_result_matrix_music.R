prepare_result_matrix_music = function(
    prediction_res_coeff_list,
    deconvolution_data,
    models,
    subtype_cands
){
    
    rounding_precision = 1
    result_matrix = data.frame(
        row.names = colnames(deconvolution_data),
        "model" = rep( paste0(c(models),collapse="|"), ncol(deconvolution_data))
    )
    
    for (nr_fit in 1:length(models)){
        
        if (nr_fit == 1){
            
            res_coeff_1 = prediction_res_coeff_list[[nr_fit]]
            colnames(res_coeff_1) = str_replace_all(colnames(res_coeff_1), 
                                                    "^X", "")
            
            res_coeff_1[ is.na(res_coeff_1) ] = 0.0
            
            # important sanity check
            #if ("alpha" %in% str_to_lower(rownames(res_coeff_1))) 
            #    res_coeff_1 = t(res_coeff_1)
            
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
    
    return(result_matrix)
}

prepare_result_matrix_music = function(
    prediction_res_coeff_list,
    deconvolution_data,
    models
){
    
    rounding_precision = 1
    result_matrix = data.frame(
        row.names = colnames(deconvolution_data),
        "Model" = rep( paste0(c(models),collapse="|"), ncol(deconvolution_data))
    )
    
    subtype_cands = c("alpha","beta","gamma","delta","acinar","ductal","hisc")
    
    for (nr_fit in 1:length(models)){
        
        if (nr_fit == 1){
            
            res_coeff_1 = prediction_res_coeff_list[[nr_fit]]
            colnames(res_coeff_1) = str_replace_all(colnames(res_coeff_1) ,"^X","")
            
            res_coeff_1[ is.na(res_coeff_1) ] = 0.0
            
            if ("alpha" %in% str_to_lower(rownames(res_coeff_1))) # important sanity check
                res_coeff_1 = t(res_coeff_1)
            
            colnames(res_coeff_1) = str_to_lower(colnames(res_coeff_1))
            subtype_cands_1 = subtype_cands[subtype_cands %in% colnames(res_coeff_1)]
            
            row_sums = rowSums(res_coeff_1[,subtype_cands[subtype_cands %in% colnames(res_coeff_1)]])
            
            for (subtype in subtype_cands_1){
                
                if (!(subtype %in% colnames(res_coeff_1)) ) next
                
                result_matrix[ , subtype] = round((res_coeff_1[,subtype] / row_sums)*100,rounding_precision)
                result_matrix[result_matrix[ , subtype] > 100,subtype] = 100
            }
        }
        ## end fit 1
        
        if ( nr_fit == 2 ){
            
            res_coeff_2 = prediction_res_coeff_list[[nr_fit]]
            colnames(res_coeff_2) = str_replace_all(colnames(res_coeff_2) ,"^X","")
            res_coeff_2[ is.na(res_coeff_2) ] = 0.0
            
            res_coeff_2 = t(res_coeff_2)
            colnames(res_coeff_2) = str_to_lower(colnames(res_coeff_2))
            
            if ("alpha" %in% rownames(res_coeff_2)) # sanity check
                res_coeff_2 = t(res_coeff_2)
            
            subtype_cands_2 = subtype_cands[
                subtype_cands %in% colnames(res_coeff_2)
                ]
            subtype_cands_2 = subtype_cands_2[!(subtype_cands_2 %in% subtype_cands_1) ]
            row_sums = rowSums(res_coeff_2[,subtype_cands[subtype_cands %in% colnames(res_coeff_2)]])
            
            for (subtype in subtype_cands_2){
                
                result_matrix[ , subtype] = round((as.double(res_coeff_2[,subtype]) / row_sums)*100,1)
                result_matrix[result_matrix[ , subtype] > 100,subtype] = 100
            }
        }
        ## end fit 2
    }
    
    colnames(result_matrix)[colnames(result_matrix) == "progenitor"] = "Progenitor"
    colnames(result_matrix)[colnames(result_matrix) == "hisc"] = "HISC"
    colnames(result_matrix)[colnames(result_matrix) == "alpha"] = "Alpha"
    colnames(result_matrix)[colnames(result_matrix) == "beta"] = "Beta"
    colnames(result_matrix)[colnames(result_matrix) == "gamma"] = "Gamma"
    colnames(result_matrix)[colnames(result_matrix) == "delta"] = "Delta"
    colnames(result_matrix)[colnames(result_matrix) == "acinar"] = "Acinar"
    colnames(result_matrix)[colnames(result_matrix) == "ductal"] = "Ductal"
    
    return(result_matrix)
}

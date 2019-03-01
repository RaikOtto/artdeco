prepare_result_matrix_bseqsc = function(
    prediction_res_coeff_list,
    deconvolution_data,
    models
){
    
    result_matrix = data.frame(
        row.names = colnames(deconvolution_data),
        "Model" = rep( paste0(c(models),collapse="|"), ncol(deconvolution_data))
    )
    
    subtype_cands = c("alpha","beta","gamma","delta","acinar","ductal","hisc")
    
    for (nr_fit in 1:length(models)){
        
        if (nr_fit == 1){
            
            res_coeff = prediction_res_coeff_list[[nr_fit]]
            colnames(res_coeff) = str_replace_all(colnames(res_coeff) ,"^X","")
            
            res_coeff[ is.na(res_coeff) ] = 0.0
            
            colnames(res_coeff) = str_to_lower(colnames(res_coeff))
            subtype_cands_1 = subtype_cands[subtype_cands %in% colnames(res_coeff)]
            
            if ("alpha" %in% str_to_lower(rownames(res_coeff))) # important sanity check
                res_coeff = t(res_coeff)
            
            for (subtype in subtype_cands){
                
                if (!(subtype %in% colnames(res_coeff)) ) next
                
                result_matrix[ , subtype] = round(res_coeff[,subtype]*100,1)
                result_matrix[result_matrix[ , subtype] > 100,subtype] = 100
            }
        }
        ## end fit 1
        
        if ( nr_fit == 2 ){
            
            res_coeff = prediction_res_coeff_list[[nr_fit]]
            colnames(res_coeff) = str_replace_all(colnames(res_coeff) ,"^X","")
            res_coeff[ is.na(res_coeff) ] = 0.0
            
            
            res_coeff = t(res_coeff)
            colnames(res_coeff) = str_to_lower(colnames(res_coeff))
            
            if ("alpha" %in% rownames(res_coeff)) # sanity check
                res_coeff = t(res_coeff)
            
            subtype_cands_2 = subtype_cands[
                subtype_cands %in% colnames(res_coeff)
                ]
            subtype_cands_2 = subtype_cands_2[!(subtype_cands_2 %in% subtype_cands_1) ]
            
            for (subtype in subtype_cands_2){
                
                result_matrix[ , subtype] = res_coeff[,subtype] * 100
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

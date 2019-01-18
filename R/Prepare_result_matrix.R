
prepare_result_matrix = function(
    prediction_res_coeff_list,
    prediction_stats_list,
    parameter_list,
    models_list,
    p_value_threshold,
    scale_values = TRUE,
    relative,
    deconvolution_data
){
    result_matrix = data.frame(
        row.names = colnames(deconvolution_data),
        "Model" = rep( paste0(c(models),collapse="|"), ncol(deconvolution_data))
    )
    for (nr_fit in 1:2){
        
        if (nr_fit == 1){
            
            res_coeff = prediction_res_coeff_list[[nr_fit]]
            res_cor   = prediction_stats_list[[nr_fit]]
            colnames(res_coeff) = str_replace_all(colnames(res_coeff) ,"^X","")
            rownames(res_cor) = str_replace_all(rownames(res_cor) ,"^X","")
            
            res_coeff[ is.na(res_coeff) ] = 0.0
            res_cor[ is.na(res_cor) ] = 0.0

            res_coeff = t(res_coeff)
            colnames(res_coeff) = str_to_lower(colnames(res_coeff))
            
            if ("alpha" %in% rownames(res_coeff))
                res_coeff = t(res_coeff)
            
            result_matrix$Differentation_score = rep(0,nrow(result_matrix))

            training_baseline = parameter_list[[nr_fit]]
            training_baseline = training_baseline[[1]]
            names(training_baseline) = str_to_lower(names(training_baseline))

            for (subtype in c("alpha","beta","gamma","delta","acinar","ductal")){
                
                if (!(subtype %in% colnames(res_coeff)) ) next
                
                if (scale_values) {
                    result_matrix[ , subtype] = exp( res_coeff[,subtype] )
                    result_matrix$Differentation_score = exp(res_cor[,4])
                } else {
                    result_matrix[ , subtype] = res_coeff[,subtype]
                    result_matrix$Differentation_score = res_cor[,4]
                }
                
                subtype_label = paste(subtype,"similarity",sep  ="_")
                result_matrix[,subtype_label] = rep("low",nrow(result_matrix))
            
                if (relative != TRUE){
                
                    baseline = as.double(as.character(unlist(
                        training_baseline[subtype])))
                    result_matrix[ , subtype] = round(
                        (as.double(result_matrix[ , subtype])/
                             quantile(baseline)[2])*100,0)
                    result_matrix[result_matrix[ , subtype] > 100,subtype] = 100
                    
                    result_matrix[
                        result_matrix[,subtype] > 25,
                        subtype_label] = "high"    
                } else {
                    
                    result_matrix[
                        result_matrix[,subtype] >
                        quantile(result_matrix[,subtype],seq(0,1,.01)[67]),
                        subtype_label
                    ] = "high"    
                }
                
            }

            cands = c("alpha","beta","gamma","delta","acinar","ductal")
            cands = cands[cands %in% colnames(result_matrix)]
            Differentiated_sim_scalar = 
                rowSums(result_matrix[,cands])
            result_matrix$Differentiated_similarity = rep("low",nrow(result_matrix))
            Differentiated_sim_scalar = exp(Differentiated_sim_scalar)
            result_matrix$Differentiated_similarity[
                Differentiated_sim_scalar > quantile(
                    Differentiated_sim_scalar,
                    seq(0,1,.01)[50]
                )
            ] = "high"
            
            # set p-values
            
            p_values = res_cor[,1]
            p_values[is.na(p_values)]  = 1
            p_values[p_values == 9999] = 0
            no_sig_index = 
                p_values >= p_value_threshold
            
            index = colnames(result_matrix)[
                str_to_lower(colnames(result_matrix)) %in%
                    c(
                        "alpha_similarity",
                        "beta_similarity",
                        "gamma_similarity",
                        "delta_similarity",
                        "acinar_similarity",
                        "ductal_similarity"
                    )
            ]
            
            result_matrix[no_sig_index,index] = "not_sig"
            }
        ## end fit 1
        
        if ( nr_fit == 2 ){
            
            res_coeff = prediction_res_coeff_list[[nr_fit]]
            res_cor   = prediction_stats_list[[nr_fit]]
            colnames(res_coeff) = str_replace_all(colnames(res_coeff) ,"^X","")
            rownames(res_cor) = str_replace_all(rownames(res_cor) ,"^X","")
            colnames(res_coeff) = str_to_lower(colnames(res_coeff))
            
            res_coeff[ is.na(res_coeff) ] = 0.0
            res_cor[ is.na(res_cor) ] = 0.0
            
            if (relative != TRUE){
                
                training_baseline = parameter_list[[nr_fit]]
                training_baseline = training_baseline[[1]]
                names(training_baseline) = str_to_lower(names(training_baseline))
                
                for (subtype in c("progenitor","hesc","hisc")){
                    
                    if ( !( subtype %in% colnames(res_coeff) ) )
                        next
                    
                    baseline = as.double(as.character(unlist(
                        training_baseline[subtype])))
                    res_coeff[,subtype] = round(
                        (res_coeff[,subtype]/
                        quantile(baseline)[1])*100,0)
                    res_coeff[res_coeff[,subtype] > 100,subtype] = 100
                }
            }
            
            result_matrix$De_differentation_score = rep(0, nrow(result_matrix))
            result_matrix$De_differentation_score = res_cor[,4]

            for ( label in colnames(res_coeff) ) {
            
                similarity_label = paste(label,"similarity",sep="_")
                
                if(scale_values){
                    sim_scalar = exp( res_coeff[,label] )
                    result_matrix$De_differentation_score = exp(result_matrix$De_differentation_score)
                } else {
                    sim_scalar = res_coeff[,label]
                }
                
                result_matrix[,label] = rep("low",nrow(result_matrix))
                result_matrix[,label] = as.double(sim_scalar)
                result_matrix[,similarity_label]   = rep("low",nrow(result_matrix))
                result_matrix[
                    sim_scalar > quantile(
                        sim_scalar,
                        seq(0,1,.01)[80]
                    ),similarity_label
                ] = "high"
            }
            
            p_values = res_cor[,1]
            p_values[is.na(p_values)]  = 1
            p_values[p_values == 9999] = 0
            no_sig_index = 
                p_values >= p_value_threshold
            
            index = colnames(result_matrix)[
                str_to_lower(colnames(result_matrix)) %in%
                    c(
                        "hisc_similarity",
                        "hesc_similarity",
                        "progenitor_similarity"
                    )
            ]
            result_matrix[no_sig_index,index] = "not_sig"
        }
        ## end fit 2
    }

    result_matrix$Differentation_score[
        result_matrix$Differentation_score >= quantile(result_matrix$Differentation_score)[4]
    ] = quantile(result_matrix$Differentation_score)[4]
    result_matrix$Differentation_score[
        result_matrix$Differentation_score <= quantile(result_matrix$Differentation_score)[1]
        ] = quantile(result_matrix$Differentation_score)[2]
    
    result_matrix$De_differentation_score[
        result_matrix$De_differentation_score >= quantile(result_matrix$De_differentation_score)[4]
        ] = quantile(result_matrix$De_differentation_score)[4]
    result_matrix$De_differentation_score[
        result_matrix$De_differentation_score <= quantile(result_matrix$De_differentation_score)[1]
        ] = quantile(result_matrix$De_differentation_score)[2]
    
    colnames(result_matrix)[colnames(result_matrix) == "progenitor_similarity"] = "Progenitor_similarity"
    colnames(result_matrix)[colnames(result_matrix) == "hisc_similarity"] = "HISC_similarity"
    colnames(result_matrix)[colnames(result_matrix) == "hesc_similarity"] = "HESC_similarity"
    
    colnames(result_matrix)[colnames(result_matrix) == "progenitor"] = "Progenitor"
    colnames(result_matrix)[colnames(result_matrix) == "hisc"] = "HISC"
    colnames(result_matrix)[colnames(result_matrix) == "hesc"] = "HESC"
    colnames(result_matrix)[colnames(result_matrix) == "alpha"] = "Alpha"
    colnames(result_matrix)[colnames(result_matrix) == "beta"] = "Beta"
    colnames(result_matrix)[colnames(result_matrix) == "gamma"] = "Gamma"
    colnames(result_matrix)[colnames(result_matrix) == "delta"] = "Delta"
    colnames(result_matrix)[colnames(result_matrix) == "acinar"] = "acinar"
    colnames(result_matrix)[colnames(result_matrix) == "ductal"] = "Ductal"
    colnames(result_matrix)[colnames(result_matrix) == "alpha_similarity"] = "Alpha_similarity"
    colnames(result_matrix)[colnames(result_matrix) == "beta_similarity"] = "Beta_similarity"
    colnames(result_matrix)[colnames(result_matrix) == "gamma_similarity"] = "Gamma_similarity"
    colnames(result_matrix)[colnames(result_matrix) == "delta_similarity"] = "Delta_similarity"
    colnames(result_matrix)[colnames(result_matrix) == "acinar_similarity"] = "acinar_similarity"
    colnames(result_matrix)[colnames(result_matrix) == "ductal_similarity"] = "Ductal_similarity"
    colnames(result_matrix)[colnames(result_matrix) == "progenitor_similarity"] = "Progenitor_similarity"
    colnames(result_matrix)[colnames(result_matrix) == "hisc_similarity"] = "HISC_similarity"
    colnames(result_matrix)[colnames(result_matrix) == "hesc_similarity"] = "HESC_similarity"
    
    cand_labels = c("alpha","beta","gamma","delta","hisc","hesc","progenitor")
    cand_label_indices = which( str_to_lower(colnames(result_matrix)) %in% cand_labels)
    
    maxi = apply( 
        result_matrix[,
                  cand_label_indices
                  ],
        FUN = which.max,
        MARGIN = 1
        
    )
    
    cand_labels_max = colnames(result_matrix)[cand_label_indices]
    cand_labels_max = cand_labels_max[maxi]
    result_matrix$Differentiation_stage = rep("",nrow(result_matrix))
    result_matrix$Differentiation_stage = cand_labels_max
    
    #p_value differentiation stage
    index = colnames(result_matrix)[
        str_to_lower(colnames(result_matrix)) %in%
            c(
                "Differentiation_stage",
                "Differentiated_similarity"
            )
        ]
    
    no_sig_index = apply(result_matrix,MARGIN =1, FUN = function(vec){
        return("not_sig" %in% vec)
    })
    
    result_matrix[no_sig_index,index] = "not_sig"

    return(result_matrix)
}

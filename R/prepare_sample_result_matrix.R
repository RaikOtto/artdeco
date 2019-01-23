
prepare_sample_result_matrix = function(
    deconvolution_results,
    prediction_stats_list,
    models_list,
    transcriptome_file
){
    
    models = as.character(unlist(str_split(deconvolution_results$Model[1],pattern = "\\|")))
    deconvolution_results[,"P_value_subtype"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Strength_subtype"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Subtype"] = rep("",nrow(deconvolution_results))

    ###
    
    res_cor = prediction_stats_list[[1]]
    rownames(res_cor) = str_replace_all(rownames(res_cor) ,"^X","")
    res_cor[ is.na(res_cor) ] = 0.0
    training_mat_dif = as.data.frame(models_list[[1]][[1]])

    transcriptome_file = transcriptome_file[rownames(training_mat_dif),]
    training_mat_dif = training_mat_dif[!is.na(transcriptome_file[,1]),]
    transcriptome_file = transcriptome_file[!is.na(transcriptome_file[,1]),]
    
    cands_dif = c("alpha","beta","gamma","delta","acinar","ductal")
    cands_dif = cands_dif[cands_dif %in% colnames(deconvolution_results)]
    
    for( j in 1:ncol(transcriptome_file)){
        
        dif_p_values = c()
        for (subtype in cands_dif){
            
            dif_p_values = c(
                cor.test(
                    transcriptome_file[,j],
                    training_mat_dif[,subtype]
                )$p.value,
                dif_p_values
            )
        }
        
        max_subtype = colnames(deconvolution_results[cands_dif])[
            which.max(deconvolution_results[j,cands_dif])
        ]
        subtype_strength = deconvolution_results[j,max_subtype] /
            sum(deconvolution_results[j,cands_dif])
        if (sum(deconvolution_results[j,cands_dif]) == 0)
            subtype_strength = deconvolution_results[j,max_subtype]
        subtype_strength = round(subtype_strength * 100,1)
        
        deconvolution_results[j,"P_value_subtype"] =
            dif_p_values[colnames(training_mat_dif) == max_subtype]
        deconvolution_results[j,"P_value_subtype"] = -1*logb(
            as.double(deconvolution_results[j,"P_value_subtype"])
        )
        deconvolution_results[j,"Strength_subtype"] =
            subtype_strength
        
        if (subtype_strength == 0)
            max_subtype = "not_significant"
        
        deconvolution_results[j,"Subtype"] = max_subtype
    }
    
    deconvolution_results[,"Differentiation_score"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"Differentiation_score"] = res_cor[,4]
    
    ### de dif
    
    res_cor = prediction_stats_list[[2]]
    rownames(res_cor) = str_replace_all(rownames(res_cor) ,"^X","")
    res_cor[ is.na(res_cor) ] = 0.0
    training_mat_de_dif = as.data.frame(models_list[[2]][[1]])
    
    transcriptome_file = transcriptome_file[rownames(training_mat_de_dif),]
    training_mat_de_dif = training_mat_de_dif[!is.na(transcriptome_file[,1]),]
    transcriptome_file = transcriptome_file[!is.na(transcriptome_file[,1]),]
    
    cands_de_dif = c("progenitor","hisc","hesc")
    deconvolution_results[,"Strength_de_differentiation"] = rep("",nrow(deconvolution_results))
    cands_de_dif = cands_de_dif[cands_de_dif %in% colnames(deconvolution_results)]
    
    for( j in 1:ncol(transcriptome_file)){
        
        de_dif_p_values = c()
        for (subtype in cands_de_dif){
            
            de_dif_p_values = c(
                cor.test(
                    transcriptome_file[,j],
                    training_mat_de_dif[,subtype]
                )$p.value,
                de_dif_p_values
            )
        }
        
        de_strength = logb( sum(deconvolution_results[j,cands_de_dif]) /
            sum(deconvolution_results[j,cands_dif]) )
        if (sum(deconvolution_results[j,cands_dif]) == 0)
            de_strength = 0
        deconvolution_results[j,"Strength_de_differentiation"] = 
            as.double(de_strength)

        deconvolution_results[j,"P_value_de_dif"] =
            de_dif_p_values[ which.min(de_dif_p_values) ]
        deconvolution_results[j,"P_value_de_dif"] = deconvolution_results[j,"P_value_de_dif"]
    }
    
    deconvolution_results[,"De_differentiation_score"] = rep("",nrow(deconvolution_results))
    deconvolution_results[,"De_differentiation_score"] = res_cor[,4]

    return(deconvolution_results)
}

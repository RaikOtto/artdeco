prepare_result_matrix_bseqsc = function(
    prediction_res_coeff_list,
    deconvolution_data,
    prediction_stats_list,
    models
){

    rounding_precision = 2
    #subtype_cands = c("alpha","beta","gamma","delta","acinar","ductal","hisc")
    subtype_cands = c()
    for (model in models){
        cands = rownames(as.data.frame(prediction_res_coeff_list[model]))
        subtype_cands = c(subtype_cands, cands)
    }
    subtype_cands = unique(subtype_cands)
    bseq_parameter = c("P_value","Correlation","RMSE")
    
    result_matrix_template = matrix(
        rep(0.0, length(subtype_cands) * ncol(deconvolution_data)
    ),ncol = length(subtype_cands))
    colnames(result_matrix_template) = subtype_cands
    result_matrix_template = as.data.frame(result_matrix_template)
    rownames(result_matrix_template) = colnames(deconvolution_data)
    result_matrix_template_ori = result_matrix_template

    # move last column to first column
    #result_matrix_template <- result_matrix_template[,c(
    #    ncol(result_matrix_template),1:(ncol(result_matrix_template)-1))]
    
    for( i in 1:length(bseq_parameter)){
        result_matrix_template[,bseq_parameter[i]] = 
            rep("",nrow(result_matrix_template))
    }

    res_coeff = prediction_res_coeff_list[[1]]
    colnames(res_coeff) = str_replace_all(colnames(res_coeff), 
                                          "^X", "")
    res_coeff[ is.na(res_coeff) ] = 0.0
    
    res_coeff = t(res_coeff)
    colnames(res_coeff) = str_to_lower(colnames(res_coeff))

    matcher = match(str_to_lower(subtype_cands), str_to_lower(colnames(res_coeff)))
    
    result_matrix_template[,
        subtype_cands
    ] = res_coeff[,matcher]
    
    prediction_stats = as.data.frame(prediction_stats_list[1])
    result_matrix_template[,
        bseq_parameter
    ] = prediction_stats[,1:ncol(prediction_stats)]

    # add code from prepare_sample_result_matrix_NMF
    # so that columns subtype and strength_subtype are added
    # can be added; doesn't have to be added
    
    
    result_matrix_template$RMSE <- as.numeric(
        result_matrix_template$RMSE)
    result_matrix_template$P_value <- as.numeric(
        result_matrix_template$P_value)
    result_matrix_template$Correlation <- as.numeric(
        result_matrix_template$Correlation)
    
    result_matrix_template[,colnames(result_matrix_template) != "P_value"] =
        round(result_matrix_template[,
        colnames(result_matrix_template) != "P_value"],
        rounding_precision)
    
    return(result_matrix_template)
}

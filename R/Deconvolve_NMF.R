#' Deconvolve_NMF
#' @param deconvolution_data Data to be deconvolved
#' @param models_list List of models
#' @param models Models utilized
#' @param nr_permutations Amount of permutations 
#' @import NMF
Deconvolve_NMF = function(
    deconvolution_data,
    models_list,
    models,
    nr_permutations
){
    
    prediction_res_coeff_list = list()
    prediction_stats_list = list()
    
    for (pred_model in models){
        
        print(paste0("Deconvolving with model: ",pred_model))
        model_basis = models_list[[pred_model]]
        W = as.data.frame(model_basis@fit@W)
         
        dec_data = as.data.frame(exprs(deconvolution_data))
        if (rownames(dec_data)[1] == "1"){
            stop("Genenames of input data missing, only found '1' as first gene name.")
        }
        W = W[rownames(W) %in% rownames(dec_data),]
        dec_data = dec_data[rownames(dec_data) %in% rownames(W),]
        dec_data = dec_data[rownames(W),]

        proportions = t(W)[,0]
        residuals   = c()
        
        for (i in 1:ncol(dec_data)){
        
            H = t(1/W)%*%dec_data[,i]
            predictions = rowSums(W*(as.matrix(H)))
            row_residuals = dec_data[,1] - predictions
            residuals = c(
                residuals,
                sum(row_residuals)
            )
            percentages = (H / sum(H))
            proportions = cbind(proportions, percentages)
        }
        colnames(proportions) = colnames(dec_data)
        names(residuals) = colnames(dec_data)
        
        prediction_res_coeff_list[[pred_model]] = proportions
        prediction_stats_list[[pred_model]]     = residuals
    }
    
    # create results matrix called meta_data
    
    deconvolution_results = prepare_result_matrix_NMF(
        prediction_res_coeff_list = prediction_res_coeff_list,
        prediction_stats_list = prediction_stats_list,
        deconvolution_data = deconvolution_data,
        models = models
    )
    colnames(deconvolution_results) = str_to_lower(colnames(deconvolution_results_nmf))
    
    #deconvolution_results = prepare_sample_result_matrix_NMF(
    #    deconvolution_results = deconvolution_results,
    #    prediction_stats_list = prediction_stats_list,
    #    models_list = models_list,
    #    deconvolution_data = exprs(deconvolution_data)
    #)
    return(deconvolution_results)
}

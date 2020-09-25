#' Deconvolve_music
#'
#' \code{Deconvolve_music} Utilizes MuSiC for deconvolution. Wrapper function
#'
#' @param deconvolution_data A data frame that contains the gene expression data.
#' Rows are expected to be HGNC symbols and columns are expected to contain the samples.
#' @param models_list List of models to be used.
#' @param models Which model to use
#' @param nr_permutations Amount perturbations
#' @import xbioc
Deconvolve_music = function(
    deconvolution_data,
    models_list,
    models,
    nr_permutations
){
    
    prediction_res_coeff_list = list()
    
    for (pred_model in models){
        
        print(paste0("Deconvolving with model: ",pred_model))
        model_basis = models_list[[pred_model]]
        model_basis = model_basis[[1]]
        pData(model_basis)$sampleID = colnames(exprs(model_basis))
        subtypes = as.character(unique(pData(model_basis)$cellType))
        
        library("xbioc")
        Est.prop.GSE50244 = MuSiC::music_prop(
            bulk.eset = deconvolution_data,
            sc.eset = model_basis,
            clusters = 'cellType',
            samples = 'sampleID',
            select.ct = subtypes,
            verbose = FALSE,
            iter.max = nr_permutations
        )

        prediction_res_coeff_list[[pred_model]] = Est.prop.GSE50244$Est.prop.allgene # cell type proportions
    }
    
    # create results matrix called meta_data
    
    deconvolution_results = prepare_result_matrix_music(
        prediction_res_coeff_list = prediction_res_coeff_list,
        deconvolution_data = deconvolution_data,
        models = models
    )
    
    selection_candidates = c(
        "model","alpha", "beta", "gamma",
        "delta", "acinar", "ductal", "hisc",
        "Strength_subtype", "Subtype")
    selection_candidates = selection_candidates[
        selection_candidates %in% colnames(deconvolution_results)]
    
    col_idx <- match(c(
        selection_candidates
    ), 
        colnames(deconvolution_results)
    )
    deconvolution_results = deconvolution_results[,col_idx]
    
    return(deconvolution_results)
}

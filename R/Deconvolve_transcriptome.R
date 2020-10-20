#' Deconvolve_transcriptome
#'
#' \code{Deconvolve_transcriptome} measures the similarity
#' of one or many query RNA-seq or microarray samples to samples
#' with known differentiation stage contained in the training models.
#'
#' @param transcriptome_data A data frame that contains the gene 
#' expression data. Rows are expected to be HGNC symbols and columns are 
#' expected to contain the samples.
#' @param deconvolution_algorithm Which deconvolution algorithm to choose
#' from. Options: 'music','bseqsc' (CIBERSORT), 'nmf'.
#' @param models List of models to be used. Use show_models_NMF(),
#' show_models_music() or show_models_bseqsc()
#' to view available models or add new model via
#' add_deconvolution_training_model_*().
#' @param nr_permutations Utilized to calculate p-value
#' Higher amount of permutations generally lead to more
#' precise p-value estimates. Default value 1000.
#' @param output_file Path of output file. If not specified,
#' no hard-disk written output will occur.
#' @import NMF
#' @usage
#' Deconvolve_transcriptome(
#'     transcriptome_data,
#'     deconvolution_algorithm,
#'     models,
#'     nr_permutations,
#'     output_file
#' )
#' @examples
#' data("visualization_data")
#' Deconvolve_transcriptome(
#'     transcriptome_data = visualization_data
#' )
#' @return Similarity measurements of differentiation
#' stages
#' @import NMF Biobase
#' @export
Deconvolve_transcriptome = function(
    transcriptome_data,
    deconvolution_algorithm = "nmf",
    models = c(
        "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",
        "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Hisc_Baron"
    ),
    nr_permutations = 1000,
    output_file = ""
){

    # data cleansing
    colnames(transcriptome_data) =
        stringr::str_replace_all(
            colnames(transcriptome_data),
            pattern = "^X",
            ""
    )
    rownames(transcriptome_data) = str_to_upper(rownames(transcriptome_data))
    transcriptome_data = variance_filtering( expr_raw = transcriptome_data)
    
    deconvolution_data = new(
        "ExpressionSet",
        exprs = as.matrix(
            transcriptome_data
        )
    )

    models_list = list()
    parameter_list = list()
    
    deconvolution_algorithm = str_to_lower(deconvolution_algorithm)

    if (deconvolution_algorithm == "music"){
        
        model_indicator = "Models/music"
    } else if (deconvolution_algorithm == "bseqsc"){
        
        model_indicator = "Models/bseqsc"
        
    } else if (deconvolution_algorithm == "nmf"){
        
        model_indicator = "Models/NMF"
    } else {
        stop("Only bseqsc, NMF and MuSiC implemented as of now.")
    }
    
    models = as.character(unlist(models))
    for (model in models){
        
        model_path = paste0(c(
            system.file(
                model_indicator,
                package="artdeco"),
                "/",
                model,".RDS"
            ),
            collapse = ""
        )
        if( !file.exists(model_path))
            stop(paste0(c("Could not find models ",model_path,", aborting"), 
                        collapse = ""))
        
        model_and_parameter = readRDS(model_path)
        
        if (deconvolution_algorithm != "nmf"){
            
            models_list[[model]] = model_and_parameter[1]
            parameter_list[[model]] = model_and_parameter[2]
        } else {
            models_list[[model]] = model_and_parameter
            parameter_list[[model]] = model_and_parameter
        }
    }
    print("Models loaded")

    ### switch algorithms
    
    if ( deconvolution_algorithm == "music"){
        
        deconvolution_results = Deconvolve_music(
            deconvolution_data = deconvolution_data,
            models_list = models_list,
            models = models,
            nr_permutations = nr_permutations
        )
        
    } else if ( deconvolution_algorithm == "bseqsc"){
        
        deconvolution_results = Deconvolve_bseq_sc(
            deconvolution_data = deconvolution_data,
            models_list = models_list,
            models = models,
            nr_permutations = nr_permutations
        )
        
    } else if ( deconvolution_algorithm == "nmf"){
        
        deconvolution_results = Deconvolve_NMF(
            deconvolution_data = deconvolution_data,
            models_list = models_list,
            models = models,
            nr_permutations = nr_permutations
        )
        
    } else {
        stop("Algorithm type not recognized. Options are music, bseqsc 
              and nmf, aborting.")
    }

    if ( output_file != "" ){
        message(paste("Writing output to file: ",output_file, sep =""))
        write.table(
            deconvolution_results,
            output_file,
            sep ="\t",
            row.names = FALSE,
            quote = FALSE
        )
    }
    message("Finished")
    return(deconvolution_results)
}

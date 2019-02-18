#' Determine_differentiation_stage
#'
#' \code{Determine_differentiation_stage} measures the similarity
#' of one or many query RNA-seq or microarray samples to samples
#' with known differentiation stage contained in the training models.
#'
#' @param transcriptome_file File which contains
#' the transcriptome data. Notice the HGNC row convention:
#' Rownames have to include the unique HGNC identifier, see vignette.
#' @param deconvolution_algorithm Which deconvolution algorithm to choose
#' from. Options: 'music','bseqsc' (CIBERSORT), 'centroid' 
#' @param models List of models to be used. Use show_models()
#' to view available models or add new model via
#' add_deconvolution_training_model()
#' @param nr_permutations Utilized to calculate p-value
#' Higher amount of permutations generally lead to more
#' precise p-value estimates
#' @param output_file Path of output file. If not specified,
#' no hard-disk written output will occur.
#' @import bseqsc MuSiC
#' @usage
#' Determine_differentiation_stage(
#'     transcriptome_file,
#'     models,
#'     nr_permutations,
#'     output_file
#' )
#' @examples
#' transcriptome_file = data("Deco_mat")
#' Determine_differentiation_stage(
#'     transcriptome_file = transcriptome_file
#' )
#' @return Similarity measurements of differentiation
#' stages
#' @export
Determine_differentiation_stage = function(
    transcriptome_file,
    deconvolution_algorithm = "music",
    models = c(
        "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",
        "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron_progenitor_stanescu_hisc_haber"
    ),
    nr_permutations = 1000,
    output_file = ""
){

    # data cleansing
    colnames(transcriptome_file) =
        stringr::str_replace_all(
            colnames(transcriptome_file),
            pattern = "^X",
            ""
    )
    rownames(transcriptome_file) = str_to_upper(rownames(transcriptome_file) )
    deconvolution_data = new("ExpressionSet", exprs=as.matrix(transcriptome_file));

    models_list = list()
    parameter_list = list()

    if (deconvolution_algorithm == "music"){
        model_indicator = "Models/music"
    } else if (deconvolution_algorithm == "bseqsc"){
        model_indicator = "Models/bseqsc"
    } else {
        stop("Only bseqsc and music implemented as of now.")
    }
    
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
            stop(paste0(c("Could not find models ",model_path,", aborting"), collapse = ""))
        model_and_parameter = readRDS(model_path)
        models_list[[model]] = model_and_parameter[1]
        parameter_list[[model]] = model_and_parameter[2]
    }
    print("Model(s) loaded")

    ### switch algorithms
    
    if (deconvolution_algorithm == "music"){
        
        deconvolution_results = Deconvolve_music(
            deconvolution_data = deconvolution_data,
            models_list = models_list,
            models = models,
            nr_permutations = nr_permutations
        )
        
    } else if (deconvolution_algorithm == "bseqsc"){
        
        deconvolution_results = Deconvolve_bseq_sc(
            deconvolution_data = deconvolution_data,
            models_list = models_list,
            models = models,
            nr_permutations = nr_permutations
        )
        
    } else if (deconvolution_algorithm == "centroid"){
        stop("Centroid no implemented as of now")
    } else {
        stop("Algorithm type not recognized, aborting.")
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

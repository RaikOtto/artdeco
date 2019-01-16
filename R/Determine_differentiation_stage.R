#' Determine_differentiation_stage
#'
#' \code{Determine_differentiation_stage} measures the similarity
#' of one or many query RNA-seq or microarray samples to samples
#' with known differentiation stage contained in the training models.
#'
#' @param transcriptome_file_path Path to the file which contains
#' the transcriptome data. Notice the HGNC row convention:
#' Rownames have to include the unique HGNC identifier, see vignette.
#' @param deconvolve_exokrine_tissue If TRUE, will deconvolve as well
#' into ductal and accinar tissue percentages
#' @param models List of models to be used. Use show_models()
#' to view available models or add new model via
#' add_deconvolution_training_model()
#' @param nr_permutations Utilized to calculate p-value
#' Higher amount of permutations generally lead to more
#' precise p-value estimates
#' @param output_file Path of output file. If not specified,
#' no hard-disk written output will occur.
#' @import bseqsc
#' @usage
#' Determine_differentiation_stage(
#'     transcriptome_file_path,
#'     deconvolve_exokrine_tissue
#'     models,
#'     nr_permutations,
#'     output_file
#' )
#' @examples
#' transcriptome_file_path = system.file(
#' "/Data/Expression_data/PANnen_Test_Data.tsv", package = "artdeco")
#' Determine_differentiation_stage(
#'     transcriptome_file_path = transcriptome_file_path
#' )
#' @return Similarity measurements of differentiation
#' stages
#' @export
Determine_differentiation_stage = function(
    transcriptome_file_path,
    deconvolve_exokrine_tissue = FALSE,
    models = c(),
    nr_permutations = 100,
    output_file = ""
){
    # check whether model is available
    if (length(models) > 0){
        message("Will utilize specified models")
    } else {
        models = c(
            "Alpha_Beta_Gamma_Delta_Lawlor",
            "Progenitor_Stanescu_HESC_Yan"
        )
        message(paste(
            "Will utilize the following models for deconvolution: ",
            paste(collapse = ", ",models, sep =",")
        ))
    }

    # check for input data availability
    if (!file.exists(transcriptome_file_path)){
        stop(paste0("Could not find file ",transcriptome_file_path))
    }

    # load transcriptome data
    transcriptome_file = read.table(
        transcriptome_file_path,
        sep="\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        row.names= 1
    )

    # data cleansing
    colnames(transcriptome_file) =
        stringr::str_replace_all(
            colnames(transcriptome_file),
            pattern = "^X",
            ""
    )
    rownames(transcriptome_file) = str_to_upper(rownames(transcriptome_file) )
    deconvolution_data = new("ExpressionSet", exprs=as.matrix(transcriptome_file));

    # check if meta data is available

    meta_data = data.frame(
        row.names = colnames(deconvolution_data),
        "Model" = rep( paste0(c(models),collapse="|"), ncol(deconvolution_data))
    )

    models_list = list()
    parameter_list = list()

    for (model in models){
        model_path = paste0(c(
            system.file(
                "Models/",
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

    prediction_stats_list = list()

    for (pred_model in models){

        print(paste0("Deconvolving with model: ",pred_model))
        model_basis = models_list[pred_model]

        fit = suppressMessages(
            bseqsc_proportions(
            deconvolution_data,
            model_basis,
            verbose = FALSE,
            absolute = TRUE,
            log = FALSE,
            perm = nr_permutations
        ))
        prediction_stats_list[[pred_model]] = fit$stats
    }

    # create results matrix called meta_data
    meta_data = prepare_result_matrix(
        prediction_stats_list = prediction_stats_list,
        parameter_list = parameter_list,
        meta_data = meta_data,
        models = models
    )

    if ( output_file != "" )
        write.table(
            meta_data,
            output_file,
            sep ="\t",
            row.names = FALSE,
            quote = FALSE
        )
    
    meta_data = meta_data[,colnames(meta_data) != "Sample"]
    return(meta_data)
}

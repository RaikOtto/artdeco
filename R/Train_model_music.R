#' add_deconvolution_training_model_music
#'
#' \code{add_deconvolution_training_model_music} adds a new model
#'
#' @param transcriptome_data Path to transcriptomic data to be
#' used for training. Has to contain the cell subtypes to which the
#' similarity has to be calculated. Note that the first column has
#' to contain the HGNC symbols and the header not! not the first
#' sample name but a mere label for this HGNC row.
#' @param model_name Name of the model
#' @param subtype_vector Character vector containing the subtype
#' labels of the training data samples
#' @import stringr
#' @usage
#' add_deconvolution_training_model_music(
#'     transcriptome_data,
#'     model_name,
#'     subtype_vector
#' )
#' @examples
#' data("Lawlor") # Data from Lawlor et al.
#' data(meta_data)
#' 
#' subtype_vector = meta_data$Subtype # extract the training sample subtype labels
#' add_deconvolution_training_model_music(
#'     transcriptome_data = Lawlor,
#'     model_name = "my_model",
#'     subtype_vector = subtype_vector
#' )
#' @return Stores a new model in the package directory
#' @export
add_deconvolution_training_model_music = function(
    transcriptome_data,
    model_name = "my_model",
    subtype_vector
){

    if( model_name == "")
        stop("Require model name, aborting")
    model_path = paste(
        c(system.file("Models/music", package="artdeco"),"/",model_name,".RDS"),
        collapse = ""
    )

    if (model_name == "my_model"){
        
    }else if (file.exists(model_path)){
        stop(paste0( collapse= "",
            c(
                "Modelname ",
                model_name,
                " already exists, please choose different name or delete existing model")
            )
        )
    }

    if (length(subtype_vector) == 0)
        stop(paste0("You have to provide the sample subtypes labels for model training"))
    #subtype_vector = str_to_lower(subtype_vector)

    expression_training_mat = transcriptome_data

    ### Data cleansing

    row_var = apply(expression_training_mat, FUN = var, MARGIN = 1)
    expression_training_mat = expression_training_mat[row_var != 0,]
    expression_training_mat = expression_training_mat[rowSums(expression_training_mat) >= 1,]

    eset = new(
        "ExpressionSet",
        exprs = as.matrix(expression_training_mat)
    );
    fData(eset) = data.frame( as.character(subtype_vector) )
    colnames(fData(eset)) = "cellType"
    pData(eset)$sampleID = colnames(expression_training_mat)
    pData(eset)$cellType = fData(eset)$cellType

    model = list(
        eset
    )

    print(paste0("Storing model: ", model_path))
    saveRDS(model,model_path)

    print(paste0("Finished training model: ", model_name))
}

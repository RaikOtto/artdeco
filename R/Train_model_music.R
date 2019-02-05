#' add_deconvolution_training_model
#'
#' \code{add_deconvolution_training_model} adds a new model
#'
#' @param transcriptome_data_path Path to transcriptomic data to be
#' used for training. Has to contain the cell subtypes to which the
#' similarity has to be calculated. Note that the first column has
#' to contain the HGNC symbols and the header not! not the first
#' sample name but a mere label for this HGNC row.
#' @param model_name Name of the model
#' @param subtype_vector Character vector containing the subtype
#' labels of the training data samples
#' @import stringr
#' @usage
#' add_deconvolution_training_model(
#'     transcriptome_data_path = "",
#'     model_name = "",
#'     subtype_vector
#' )
#' @examples
#' transcriptome_data_path = system.file(
#' "Data/Expression_data/PANnen_Test_Data.tsv",package ="artdeco")
#' meta_data_path = system.file("Data/Meta_information/Meta_information.tsv", package = "artdeco")
#' meta_data      = read.table(
#'     meta_data_path, sep = "\t",
#'     header = TRUE,
#'     stringsAsFactors = FALSE
#' )
#' subtype_vector = meta_data$Subtype # extract the training sample subtype labels
#' add_deconvolution_training_model(
#'     transcriptome_data_path = transcriptome_data_path,
#'     model_name = "Test_model",
#'     subtype_vector
#' )
#' @return Stores a new model in the package directory
#' @export
add_deconvolution_training_model_music = function(
    transcriptome_data_path = "",
    model_name = "",
    subtype_vector
){

    if( model_name == "")
        stop("Require model name, aborting")
    model_path = paste(
        c(system.file("Models/music", package="artdeco"),"/",model_name,".RDS"),
        collapse = ""
    )

    if(file.exists(model_path))
        stop(paste0( collapse= "",
                     c("Modelname ",model_name,
                       " already exists, please choose different name or delete existing model"))
        )

    if( ! file.exists(transcriptome_data_path)){
        stop(paste(
            c("Could not find file ",transcriptome_data_path,", aborting"),
            collapse = ""
        ))
    }

    if (length(subtype_vector) == 0)
        stop(paste0("You have to provide the sample subtypes labels for model training"))
    #subtype_vector = str_to_lower(subtype_vector)

    print("Loading training data")
    
    expression_training_mat = read.table(
        transcriptome_data_path,
        sep ="\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        row.names = 1
    )

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

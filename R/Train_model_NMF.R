#' add_deconvolution_training_model_NMF
#'
#' \code{add_deconvolution_training_model_NMF} adds a new model
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
#' add_deconvolution_training_model_NMF(
#'     transcriptome_data,
#'     model_name,
#'     subtype_vector
#' )
#' @examples
#' data("test_data")
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
#' @import NMF
#' @return Stores a new model in the package directory
#' @export
add_deconvolution_training_model_NMF = function(
    transcriptome_data_path = "",
    model_name = "",
    subtype_vector,
    training_nr_marker_genes = 100
){

    if( model_name == "")
        stop("Require model name, aborting")
    model_path = paste(
        c(system.file("Models/NMF", package="artdeco"),"/",model_name,".RDS"),
        collapse = ""
    )

    if(file.exists(model_path))
        stop(paste0( collapse= "",
                     c("Modelname ",model_name,
                       " already exists, please choose different name or delete existing model"))
        )

    if (length(subtype_vector) == 0)
        stop(paste0("You have to provide the sample subtypes labels for model training"))
    #subtype_vector = str_to_lower(subtype_vector)

    expression_training_mat = transcriptome_data

    ### Data cleansing

    row_var = apply(expression_training_mat, FUN = var, MARGIN = 1)
    expression_training_mat = expression_training_mat[row_var != 0,]
    expression_training_mat = expression_training_mat[rowSums(expression_training_mat) >= 1,]
    
    ### dif genes
    
    markers <<- c()
    Marker_Gene_List = list()
    for( subtype in unique(subtype_vector) ){
        markers = c(
            markers,
            identify_marker_genes(
                expression_training_mat = expression_training_mat,
                subtype_vector = subtype_vector,
                subtype = subtype,
                nr_marker_genes = training_nr_marker_genes
            )
        )
        Marker_Gene_List[[subtype]] = markers
    }
    markers = unique(markers)
    expression_training_mat_reduced = expression_training_mat[markers,]
    print("Finished extracting marker genes for subtypes")

    ### NMF training
    
    rank_estimate = length(unique(subtype_vector))
    library("NMF")

    res = nmf(
        expression_training_mat_reduced,
        rank = rank_estimate,
        method = 'brunet',
        .opt = 'tp3',
        nrun = 10,
        maxIter = 100,
        seed = "random"
    )
    summary(res, class = subtype_vector)
    plot(res)
    H=res@fit@H
    W=res@fit@W
    aprx = W %*% H
    residual = expression_training_mat_reduced - aprx
    residual[1:5,1:5]
    dd = NMF::basis(res)
    round(dd["INS",],0)
    
    aggregate(
        
    )

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


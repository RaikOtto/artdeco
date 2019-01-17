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
#' @param marker_gene_list List the contains the marker genes
#' for each subtype. Has to be in the type of list() with each
#' subtype being an entry.
#' @param training_p_value_threshold P-value at which a training
#' is deemed successfull
#' @param training_nr_permutations Amount of perturbation which
#' result in a p-value. Higher number of perturbation generally
#' improves the p-value estiamtes
#' @param training_nr_marker_genes How many genes should be utilized
#' as list of marker genes
#' @import stringr bseqsc
#' @usage
#' add_deconvolution_training_model(
#'     transcriptome_data_path = "",
#'     model_name = "",
#'     subtype_vector,
#'     marker_gene_list = list(),
#'     training_p_value_threshold = 0.05,
#'     training_nr_permutations = 100,
#'     training_nr_marker_genes = 100
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
#'     subtype_vector,
#'     training_p_value_threshold = 0.05,
#'     training_nr_permutations = 100,
#'     training_nr_marker_genes = 100
#' )
#' @return Stores a new model in the package directory
#' @export
add_deconvolution_training_model = function(
    transcriptome_data_path = "",
    model_name = "",
    subtype_vector,
    marker_gene_list = list(),
    training_p_value_threshold = 0.05,
    training_nr_permutations = 100,
    training_nr_marker_genes = 100
){

    if( model_name == "")
        stop("Require model name, aborting")
    model_path = paste(
        c(system.file("Models/", package="artdeco"),"/",model_name,".RDS"),
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

    if (length(marker_gene_list) <= 1){

        Marker_Gene_List <<- list()
        for( subtype in unique(subtype_vector) ){
            Marker_Gene_List[[subtype]] = identify_marker_genes(
                expression_training_mat = expression_training_mat,
                subtype_vector = subtype_vector,
                subtype = subtype,
                nr_marker_genes = training_nr_marker_genes
            )
        }
        print("Finished extracting marker genes for subtypes")
    } else {
        Marker_Gene_List = marker_gene_list
    }

    # Prepare bseq training
    training_mat_bseq = new(
        "ExpressionSet",
        exprs = as.matrix(expression_training_mat)
    )
    fData(training_mat_bseq) = data.frame( as.character(subtype_vector) )
    colnames(fData(training_mat_bseq)) = "subtype_vector"
    pData(training_mat_bseq) = data.frame( as.character(subtype_vector) )
    colnames(pData(training_mat_bseq)) = "subtype_vector"

    Basis = suppressMessages(
        bseqsc_basis(
            training_mat_bseq,
            Marker_Gene_List,
            clusters = 'subtype_vector',
            samples = colnames(expression_training_mat),
            ct.scale = FALSE
    ))

    print( paste0( c("Basis trained, estimating deconvolution thresholds,",
                     "this may take some time"), collapse = ""))

    test_mat = new(
        "ExpressionSet",
        exprs = as.matrix(expression_training_mat)
    );
    fData(test_mat) = data.frame( as.character(subtype_vector) )
    colnames(fData(test_mat)) = "subtype_vector"
    pData(test_mat) = data.frame( as.character(subtype_vector) )
    colnames(pData(test_mat)) = "subtype_vector"
    
    fit = bseqsc_proportions(
        test_mat,
        Basis,
        verbose = TRUE,
        absolute = TRUE,
        log = F,
        perm = training_nr_permutations,
        
    )

    print("Finished threshold determination")

    res_coeff = t(fit$coefficients)
    res_coeff_mat = as.double(unlist(res_coeff))
    res_coeff_mat = as.data.frame(
        matrix(
            res_coeff_mat,
            ncol = ncol(res_coeff),
            nrow = nrow(res_coeff)
        )
    )
    rownames(res_coeff_mat) = rownames(res_coeff)
    colnames(res_coeff_mat) = colnames(res_coeff)
    res_cor   = fit$stats

    res_coeff[ is.na(res_coeff) ] = 0.0
    res_cor[ is.na(res_cor) ] = 0.0

    self_scores = list()
    for (subtype in unique(subtype_vector)){
        self_scores[[subtype]] = as.double(
            res_coeff[
                which(subtype_vector == subtype),
                subtype
            ]
        )
    }

    model = list(Basis, self_scores, Marker_Gene_List)

    print(paste0("Storing model: ", model_path))
    saveRDS(model,model_path)

    print(paste0("Finished training model: ", model_name))
}

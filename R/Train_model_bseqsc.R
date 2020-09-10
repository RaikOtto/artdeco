#' add_deconvolution_training_model_bseqsc
#'
#' \code{add_deconvolution_training_model_bseqsc} adds a new model
#'
#' @param transcriptome_data Path to transcriptomic data to be
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
#' add_deconvolution_training_model_bseqsc(
#'     transcriptome_data,
#'     model_name,
#'     subtype_vector,
#'     marker_gene_list,
#'     training_p_value_threshold,
#'     training_nr_permutations,
#'     training_nr_marker_genes
#' )
#' @examples
#' data("Lawlor") # Data from Lawlor et al.
#' data(meta_data)
#' 
#' subtype_vector = as.character(meta_data$Subtype) # extract the training sample subtype labels
#' 
#' add_deconvolution_training_model_bseqsc(
#'     transcriptome_data = Lawlor,
#'     model_name = "my_model",
#'     subtype_vector = subtype_vector,
#'     training_nr_permutations = 10,
#'     training_nr_marker_genes = 100
#' )
#' @return Stores a new model in the package directory
#' @export
add_deconvolution_training_model_bseqsc = function(
    transcriptome_data,
    model_name,
    subtype_vector,
    marker_gene_list,
    training_p_value_threshold = 0.05,
    training_nr_permutations = 100,
    training_nr_marker_genes = 100
){

    if( model_name == "")
        stop("Require model name, aborting")
    model_path = paste(
        c(system.file("Models/bseqsc", package="artdeco"),"/",model_name,".RDS"),
        collapse = ""
    )

    if (model_name == "my_model"){
        return("")
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

    markers <<- c()
    Marker_Gene_List = list()
    for( subtype in unique(subtype_vector) ){
        
        markers_subtype = identify_marker_genes(
            expression_training_mat = expression_training_mat,
            subtype_vector = subtype_vector,
            subtype = subtype,
            training_nr_marker_genes = training_nr_marker_genes
        )
        markers = c(
            markers,
            markers_subtype
        )
        Marker_Gene_List[[subtype]] = markers_subtype
    }
    markers = unique(markers)
    expression_training_mat_reduced = expression_training_mat[markers,]
    print("Finished extracting marker genes for subtypes")

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
        bseqsc::bseqsc_basis(
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
    
    if ( model_name != "my_model" ){
        fit = bseqsc_proportions(
            test_mat,
            Basis,
            verbose = TRUE,
            absolute = TRUE,
            log = FALSE,
            perm = training_nr_permutations
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
    }

    print(paste0("Finished training model: ", model_name))
}

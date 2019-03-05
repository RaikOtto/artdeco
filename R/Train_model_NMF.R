#' add_deconvolution_training_model_NMF
#'
#' \code{add_deconvolution_training_model_NMF} adds a new model
#'
#' @param transcriptome_data Transcriptomic data to be
#' used for training. Has to contain the cell subtypes to which the
#' similarity has to be calculated. Note that the first column has
#' to contain the HGNC symbols and the header not! not the first
#' sample name but a mere label for this HGNC row.
#' @param model_name Name of the model
#' @param subtype_vector Character vector containing the subtype
#' labels of the training data samples
#' @param rank_estimate Rank of the NMF model. Will be set to amount of 
#' different subtypes defined in the subtype_vector if not specified manually
#' @param exclude_non_interpretable_NMF_components Boolean parameter that indicates
#' whether trained NMF components that cannot clearly be associated with either an endocrine
#' or acinar & ductal or hisc subtyp shall be excluded
#' @param training_nr_marker_genes Amount of genes to be utilized as marker genes
#' for each cell type
#' @import stringr NMF
#' @usage
#' add_deconvolution_training_model_NMF(
#' transcriptome_data,
#' model_name,
#' subtype_vector,
#' rank_estimate,
#' exclude_non_interpretable_NMF_components,
#' training_nr_marker_genes
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
#'     transcriptome_data = transcriptome_data,
#'     model_name = "Test_model",
#'     subtype_vector
#' )
#' @return Stores a new model in the package directory
#' @export
add_deconvolution_training_model_NMF = function(
    transcriptome_data,
    model_name,
    subtype_vector,
    rank_estimate,
    exclude_non_interpretable_NMF_components = TRUE,
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

    ### NMF training
    
    if ( !  exists("rank_estimate") )
        rank_estimate = length(unique(subtype_vector))
    library("NMF")
    
    print("Commencing NMF training")

    res = nmf(
        expression_training_mat_reduced,
        rank = rank_estimate,
        method = 'brunet',
        .opt = 'tp8',
        nrun = 10,
        maxIter = 100,
        seed = "random"
    )
    
    print("Finished NMF training")
    
    summary(res, class = subtype_vector)
    
    W=res@fit@W
    colnames(W) = 1:rank_estimate

    for (subtype in names(Marker_Gene_List)){
        
        marker_genes = as.character(unlist(Marker_Gene_List[subtype]))
        
        col_var = apply( W[marker_genes,],MARGIN = 2,FUN= var)
        max_index = as.integer(which.max(col_var))
        print(paste0(collapse="",c(subtype," subtype found as component ",max_index),sep=""))
        colnames(W)[max_index] = subtype
    }
    
    if(exclude_non_interpretable_NMF_components){
        exclusion_index = which(!str_to_lower(
            colnames(W)) %in% c("alpha","beta","gamma","delta","acinar","ductal","hisc"
            ))
        print(paste("Excluding component ",exclusion_index,sep=""))
        W = W[,-exclusion_index]
    }
    res@fit@W = W

    print(paste0("Storing model: ", model_path))
    saveRDS(res,model_path)

    print(paste0("Finished training model: ", model_name))
}


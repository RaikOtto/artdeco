#' create_PCA_differentiation_stages
#'
#' \code{create_PCA_differentiation_stages}
#' visualizes the differentiation stage predictions as PCA.
#' Please note that the first column of the expression
#' data matrix has to contain the HGNC identifier
#'
#' @param transcriptome_file_path Path to the transcriptome
#' data that shall be visualized. Notice the convention
#' that the first row has to contain the HGNC identifier
#' @param deconvolution_results The dataframe returned
#' by the deconvolution analysis
#' @param annotation_columns Column names that are to be visualized
#' ontop of the correlation matrix
#' @param Graphics_parameters Pheatmap visualization paramters.
#' You can customize visualization colors.
#' Read the vignette for more information.
#' @param baseline Which measurement represents the baseline
#' of the differentiation similarity: 'absolute' = maximal
#' similarity of the training sample to its subtype e.g.
#' alpha cell similarity to alpha cell =
#' @usage
#' create_PCA_differentiation_stages(
#'     transcriptome_file_path,
#'     deconvolution_results,
#'     annotation_columns,
#'     Graphics_parameters,
#'     baseline
#' )
#' @examples
#' meta_data_path = system.file(
#'     "Data/Meta_information/Meta_information.tsv",
#'     package = "artdeco"
#' )
#' deconvolution_results      = read.table(
#'     meta_data_path, sep ="\t",
#'     header = TRUE,
#'     stringsAsFactors = FALSE
#' )
#' rownames(deconvolution_results) = deconvolution_results$Sample
#'
#' transcriptome_file_path = system.file(
#'     "/Data/Expression_data/Visualization_PANnen.tsv",
#'      package = "artdeco"
#' )
#' create_PCA_differentiation_stages(
#'     transcriptome_file_path = transcriptome_file_path,
#'     deconvolution_results = deconvolution_results,
#'     annotation_columns = c(
#'         "Differentiation_Stages_Subtypes",
#'         "Differentiation_Stages_Aggregated",
#'         "Differentiatedness"
#'      ),
#'      Graphics_parameters = "",
#'      baseline = "relative"
#' )
#' @import stringr ggplot2 pheatmap ggfortify
#' @return Plots
#' @export
create_PCA_differentiation_stages = function(
    transcriptome_file_path,
    deconvolution_results,
    annotation_columns = c(
        "Differentiation_Stages_Subtypes",
        "Differentiation_Stages_Aggregated",
        "Differentiatedness"
    ),
    Graphics_parameters = "",
    baseline = "relative"
){

    # check for input data availability
    if (!file.exists(transcriptome_file_path)){
        stop(paste0("Could not find file ",transcriptome_file_path))
    }

    # load transcriptome data
    transcriptome_mat_vis = read.table(
        transcriptome_file_path,
        sep="\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        row.names = 1
    )
    colnames(transcriptome_mat_vis) = str_replace_all(
        colnames(transcriptome_mat_vis),
        pattern = "^X",
        ""
    )

    correlation_matrix = cor(transcriptome_mat_vis)
    pcr = prcomp(t(correlation_matrix))

    # ensure correct ordering of expression and annotation data
    rownames(deconvolution_results) = deconvolution_results$Sample
    deconvolution_results = deconvolution_results[
        colnames(correlation_matrix),
    ]

    # assert that graphics parameters are available
    if ( head(Graphics_parameters,1) == "" )
        Graphics_parameters = configure_graphics()

    # extract similarity measurements
    sim_index = grep(
        x = colnames(deconvolution_results),
        pattern = paste0( c(
            "(_similarity$)",
            "(Differentiation_Stage)",
            "Differentiatedness"
        ), collapse = "|")
    )
    annotation_data = deconvolution_results[, sim_index]

    # case the results data frame does not have a column
    for (annotation in annotation_columns){
        if (!( annotation %in% colnames(deconvolution_results))){
            print(paste0(c(
                "Did not find ",
                annotation,
                " column in the deconvolution matrix. No Visualization will occur."
            )))
        } else {
            annotation_data[,annotation] = rep("",nrow(annotation_data))
            annotation_data[rownames(deconvolution_results),annotation] =
                deconvolution_results[,annotation]
        }
    }

    # assert uniqueness of results
    annotation_data = annotation_data[,unique(colnames(annotation_data))]

    # ensure correct type cast
    if(length(annotation_data$Differentiatedness) > 0)
        annotation_data$Differentiatedness = as.double(
            as.character(
                annotation_data$Differentiatedness
            )
        )
    #    obs.scale = .75,
    #    groups = annotation_data[,"Differentiation_Stages_Aggregated"],
    #    ellipse = TRUE,
    #    circle = TRUE,
    #    var.axes = FALSE
    #    #,labels = meta_data$Name
    #
    p = autoplot(
        pcr,
    )
    p = p + geom_point(
        aes(
            colour = annotation_data[,"Differentiation_Stages_Aggregated"],
            shape = deconvolution_results[,"Differentiation_Stages_Aggregated"]),
        size = annotation_data$Differentiatedness /
            max(annotation_data$Differentiatedness ) * 10
    )
    p = p + scale_shape_manual(
        values = c(1,2),
        labels = c("Deconvolveable","Not deconvolveable")
    )
    p = p + scale_color_manual(
        values = c("darkgreen","darkred"),
        labels = c("Differentiated","Not differentiated")
    )
    p = p + guides( color=guide_legend(
        title="Differentiation stage",
        size = guide_legend(title="Differentiatedness")),
        shape =  guide_legend(
            title="Deconvolveability"
        )
    )
    p
}

#' create_heatmap_differentiation_stages
#'
#' \code{create_heatmap_differentiation_stages}
#' visualizes the differentiation stage predictions as heatmap.
#' Please note that the first column of the expression
#' data matrix has to contain the HGNC identifier
#'
#' @param transcriptome_file_path Path to the transcriptome
#' data that shall be visualized. Notice the convention
#' that the first row has to contain the HGNC identifier
#' @param deconvolution_results The dataframe returned
#' by the deconvolution analysis
#' @param aggregate_differentiated_stages Show the differentiation stage similarities
#' aggregated over all differentiated stages - alpha, beta, gamma, delta
#' accinar and ductal - or specific for each differentiation stage.
#' Default value FALSE, alternative value TRUE
#' @param p_value_threshold Threshold above which deconvolutions are deemed unsuccessful
#' and corresponding results being masked on the the plots
#' @param Graphics_parameters Pheatmap visualization paramters.
#' You can customize visualization colors.
#' Read the vignette for more information.
#' @param relative_baseline Which measurement represents the baseline
#' of the differentiation similarity: 'absolute' = maximal
#' similarity of the training sample to its subtype. TRUE
#' sets the baseline to maximal similarity of the test samples
#' currently analysed
#' @usage
#' create_heatmap_differentiation_stages(
#'     transcriptome_file_path,
#'     deconvolution_results,
#'     aggregate_differentiated_stages,
#'     p_value_threshold,
#'     Graphics_parameters,
#'     relative_baseline
#' )
#' @examples
#' meta_data_path = system.file(
#'     "Data/Meta_information/Meta_information.tsv",
#'     package = "artdeco"
#' )
#' meta_data      = read.table(
#'     meta_data_path, sep ="\t",
#'     header = TRUE,
#'     stringsAsFactors = FALSE
#' )
#' rownames(meta_data) = meta_data$Sample
#'
#' visualization_data_path = system.file(
#'     "/Data/Expression_data/Visualization_PANnen.tsv",
#'      package = "artdeco"
#' )
#' create_heatmap_differentiation_stages(
#'     visualization_data_path = visualization_data_path,
#'     deconvolution_results = meta_data,
#'     aggregate_differentiated_stages = FALSE,
#'     Graphics_parameters = "",
#'     p_value_threshold = 0.05,
#'     relative_baseline = TRUE
#' )
#' @import stringr ggplot2 pheatmap ggfortify
#' @return Plots
#' @export
create_heatmap_differentiation_stages = function(
    visualization_data_path,
    deconvolution_results,
    aggregate_differentiated_stages = FALSE,
    p_value_threshold = 0.05,
    Graphics_parameters = "",
    high_threshold_diff = 20,
    high_threshold_de_diff = 8.25
){

    # check for input data availability
    if (!file.exists(visualization_data_path)){
        stop(paste0("Could not find file ",visualization_data_path))
    }
    
    # data cleansing
    deconvolution_results$Subtype =
        str_to_lower(deconvolution_results$Subtype)
    
    # load transcriptome data
    transcriptome_mat_vis = read.table(
        visualization_data_path,
        sep="\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        row.names = 1
    )
    colnames(transcriptome_mat_vis) = str_replace_all(
        colnames(transcriptome_mat_vis),
        pattern = "^X",
        ""
    )
    
    # get the index for differentiated state
    
    dif_index = grep(
        colnames(deconvolution_results),
        pattern = paste0( c(
            "alpha",
            "beta",
            "gamma",
            "delta",
            "acinar",
            "ductal"
        ), collapse = "|")
    )
    
    de_dif_index = grep(
        colnames(deconvolution_results),
        pattern = paste0( c(
            "progenitor",
            "hisc",
            "hesc"
        ), collapse = "|")
    )

    correlation_matrix = cor(transcriptome_mat_vis)
    deconvolution_results = deconvolution_results[
        colnames(correlation_matrix),
    ]

    # assert that graphics parameters are available
    if ( head(Graphics_parameters,1) == "" )
        Graphics_parameters = configure_graphics()

    # extract similarity measurements
    subtype_index_vis = c
    (
        grep(
            colnames(deconvolution_results),
            pattern = paste0( c(
                "Cor_differentiated_model",
                "Cor_de_differentiated_model",
                "NEC_NET"
            ), collapse = "|")
        )
    )
    
    aggregated_index_vis = grep(
        colnames(deconvolution_results),
        pattern = paste0( c(
            "Subtype",
            "Strength_de_differentiation",
#            "Differentiation_score",
            "NEC_NET"
        ), collapse = "|")
    )
    
    if (aggregate_differentiated_stages == FALSE){
        vis_mat = deconvolution_results[subtype_index_vis]
    } else {
        vis_mat = deconvolution_results[aggregated_index_vis]
    }
    
    vis_mat$Subtype[
        as.double(deconvolution_results$P_value_subtype) <=
            p_value_threshold
    ] = "not_significant"
    vis_mat$Strength_de_differentiation = as.double(vis_mat$Strength_de_differentiation)
    
    # differentiation score scaling
    
    
    for ( score in c("Differentiation_score","De_differentiation_score")){
        
        if (! ( "Differentiation_score" %in% colnames(vis_mat)))
            next()
        
        vis_mat[,score] = (vis_mat[,score] )
        quantiles = quantile(vis_mat[, score ])
        vis_mat[ vis_mat[,score] <=  quantiles[1], score ] = quantiles[1]
        vis_mat[ vis_mat[,score] >=  quantiles[4], score ] = quantiles[4]
        
    }
    
    # correlation heatmap
    pheatmap::pheatmap(
        correlation_matrix,
        annotation_col = vis_mat,
        annotation_colors = Graphics_parameters,
        annotation_legend = TRUE,
        treeheight_col = 0,
        treeheight_row = 0,
        show_colnames = FALSE,
        show_rownames = FALSE
    )
    
    return(deconvolution_results)
}


#' configure_graphics
#'
#' \code{configure_graphics} configure the graphics paramters
#' @usage
#' configure_graphics()
#' @examples
#' graphical_configuration = configure_graphics()
#' @return Heatmap related graphics configuration
#' @export
configure_graphics = function(){

    Graphics_parameters         = list(
        NEC_NET                 = c(NEC= "red", NET = "blue"),
        Alpha_similarity        = c(Not_significant = "gray", low = "white", high = "blue"),
        Beta_similarity         = c(Not_significant = "gray", low = "white", high = "green"),
        Gamma_similarity        = c(Not_significant = "gray", low = "white", high = "brown"),
        Delta_similarity        = c(Not_significant = "gray", low = "white", high = "Purple"),
        Ductal_similarity     = c(Not_significant = "gray", low = "white", high = "Black"),
        Acinar_similarity     = c(Not_significant = "gray", low = "white", high = "Brown"),
        Progenitor_similarity = c(Not_significant = "gray", low = "white", high = "orange"),
        HISC_similarity   = c(Not_significant = "gray",low = "white", high = "black"),
        HESC_similarity   = c(Not_significant = "gray",low = "white", high = "black"),
        Strength_de_differentiation = c(low = "darkgreen", medium = "white", high = "darkred"),
        Strength_subtype = c(low = "black", high = "darkgreen"),
        De_differentiation_score = c(low = "red", medium = "white", high = "green"),
        Differentiation_score = c(low = "red", medium = "white", high = "green"),
        Subtype = c(
            alpha = "blue",
            beta = "green",
            gamma = "brown",
            delta = "purple",
            acinar = "cyan",
            ductal = "red",
            not_significant = "gray"
        ),
        Differentiation_Stage_Aggregated = c(
            differentiated = "darkgreen",
            hisc = "red",
            progenitor     = "orange",
            Not_significant= "gray"
        ),
        Grading          = c( G1 = "Green",G2 = "Yellow", G3 = "Red", G0 = "white")
    )
    return(Graphics_parameters)
}

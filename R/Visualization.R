#' create_PCA_differentiation_stages
#'
#' \code{create_PCA_differentiation_stages}
#' visualizes the differentiation stage predictions as PCA.
#' Please note that the first column of the expression
#' data matrix has to contain the HGNC identifier
#'
#' @param visualization_data_path Path to the transcriptome
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
#'     visualization_data_path,
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
#' visualization_data_path = system.file(
#'     "/Data/Expression_data/Visualization_PANnen.tsv",
#'      package = "artdeco"
#' )
#' create_PCA_differentiation_stages(
#'     visualization_data_path = visualization_data_path,
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
    visualization_data_path,
    deconvolution_results,
    Graphics_parameters = ""
){

    # check for input data availability
    if (!file.exists(visualization_data_path)){
        stop(paste0("Could not find file ",transcriptome_file_path))
    }

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

    correlation_matrix = cor(transcriptome_mat_vis)
    pcr = prcomp(t(correlation_matrix))

    # ensure correct ordering of expression and annotation data
    deconvolution_results = deconvolution_results[
        colnames(correlation_matrix),
    ]

    # assert that graphics parameters are available
    if ( head(Graphics_parameters,1) == "" )
        Graphics_parameters = configure_graphics()
    
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
            #"progenitor",
            "hisc",
            "hesc"
        ), collapse = "|")
    )
    
    # extract similarity measurements
    subtype_index_vis = c(dif_index,de_dif_index,
          grep(
              colnames(deconvolution_results),
              pattern = paste0( c(
                  "Grading",
                  "NEC_NET",
                  "Study"
              ), collapse = "|")
          )
    )
    
    vis_mat = deconvolution_results[subtype_index_vis]
    
    for(subtype in c(cands_de_dif,cands_dif)){
        vis_mat[
            deconvolution_results[,subtype] <=
                quantile(
                    deconvolution_results[,subtype],
                    probs = seq(0,1,.01)
                )[50],
            subtype
            ] = "low"
        vis_mat[
            deconvolution_results[,subtype] >
                quantile(
                    deconvolution_results[,subtype],
                    probs = seq(0,1,.01)
                )[50],
            subtype
            ] = "high"
    }

    
    p = autoplot(pcr)
    p = p + geom_point(
        aes(
            colour = vis_mat[,"hisc"],
            shape = vis_mat[,"Grading"]),
            size = 5
    )
    p = p + guides(
        color=guide_legend(
            title="Differentiation stage",
            size = guide_legend(
                title="HISC-similarity"
            )
        ),
                shape =  guide_legend(
                title="Grading"
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
#' @param visualization_data Matrix of the transcriptome data
#' data that shall be visualized. Notice the convention
#' that the first row has to contain the HGNC identifier
#' @param deconvolution_results The dataframe returned
#' by the deconvolution analysis
#' @param aggregate_differentiated_stages Show the differentiation stage similarities
#' aggregated over all differentiated stages - alpha, beta, gamma, delta
#' accinar and ductal - or specific for each differentiation stage.
#' Default value FALSE, alternative value TRUE
#' @param show_colnames Whether to show the sample column names
#' @param confidence_threshold Threshold above which deconvolutions are deemed unsuccessful
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
#'     show_colnames
#'     confidence_threshold,
#'     Graphics_parameters,
#'     relative_baseline
#' )
#' @examples
#' visualization_data = data("Vis_data)
#' meta_data      = read.table(
#'     meta_data_path, sep ="\t",
#'     header = TRUE,
#'     stringsAsFactors = FALSE
#' )
#' rownames(meta_data) = meta_data$Sample
#'
#' visualization_data = data("Vis_data")
#' create_heatmap_differentiation_stages(
#'     visualization_data = visualization_data,
#'     deconvolution_results = meta_data,
#'     aggregate_differentiated_stages = FALSE,
#'     show_colnames = FALSE,
#'     Graphics_parameters = "",
#'     confidence_threshold = 0.05,
#'     relative_baseline = TRUE
#' )
#' @import stringr ggplot2 pheatmap ggfortify
#' @return Plots
#' @export
create_heatmap_differentiation_stages = function(
    visualization_data,
    deconvolution_results,
    aggregate_differentiated_stages = FALSE,
    confidence_threshold = 1.1,
    show_colnames = FALSE,
    Graphics_parameters = "",
    high_threshold = 101
){

    # data cleansing
    deconvolution_results$Subtype =
        str_to_lower(deconvolution_results$Subtype)

    transcriptome_mat_vis = visualization_data
    colnames(transcriptome_mat_vis) = str_replace_all(
        colnames(transcriptome_mat_vis),
        pattern = "^X",
        ""
    )

    # get the index for differentiated state
    
    cands_dif = c(
        "alpha",
        "beta",
        "gamma",
        "delta",
        "acinar",
        "ductal"
    )
    cands_dif = cands_dif[cands_dif %in% colnames(deconvolution_results)]
    
    cands_de_dif = c(
        #"progenitor",
        "hisc",
        "hesc"
    )
    cands_de_dif = cands_de_dif[cands_de_dif %in% colnames(deconvolution_results)]
    
    dif_index = grep(
        colnames(deconvolution_results),
        pattern = paste0(
            cands_dif,
            collapse = "|"
        )
    )
    
    de_dif_index = grep(
        colnames(deconvolution_results),
        pattern = paste0( 
            cands_de_dif,
            collapse = "|"
        )
    )

    correlation_matrix = cor(transcriptome_mat_vis)
    deconvolution_results = deconvolution_results[
        colnames(correlation_matrix),
    ]

    # assert that graphics parameters are available
    if ( head(Graphics_parameters,1) == "" )
        Graphics_parameters = configure_graphics()

    # extract similarity measurements
    subtype_index_vis = c(dif_index,de_dif_index,
        grep(
            colnames(deconvolution_results),
            pattern = paste0( c(
                "Grading",
                "NEC_NET",
                "Confidence_score_dif",
                "Confidence_score_de_dif"
            ), collapse = "|")
        )
    )

    vis_mat = deconvolution_results[subtype_index_vis]
    
    # differentiation score scaling
    for ( score in c("Confidence_score_de_dif","Confidence_score_dif")){

        vis_mat[,score] = vis_mat[,score] - min(vis_mat[,score])
        vis_mat[,score] = vis_mat[,score] / max(vis_mat[,score])
    }
    
    if ( confidence_threshold > 1.0 )
        confidence_threshold = 
            mean(vis_mat[,"Confidence_score_dif"]) +
            sd(vis_mat[,"Confidence_score_dif"])
    
    for(subtype in c(cands_de_dif,cands_dif)){
        
        # log the result
        #deconvolution_results[,subtype] = log(deconvolution_results[,subtype]+1)
        #deconvolution_results[,subtype] = scale(deconvolution_results[,subtype])
        deconvolution_results[,subtype] = deconvolution_results[,subtype] - min(deconvolution_results[,subtype])
        deconvolution_results[,subtype] = deconvolution_results[,subtype] / max(deconvolution_results[,subtype])
        deconvolution_results[,subtype] = round(deconvolution_results[,subtype] * 100,1)
        
        if (high_threshold > 100){
            m = mean(deconvolution_results[,subtype])
            sd = sd(deconvolution_results[,subtype])*.5
            high_threshold_subtype = round(m + sd,1)
            
            if (
                ( high_threshold_subtype == 0 ) |
                ( length(high_threshold_subtype) == 0 )
            )
            high_threshold_subtype = 1
            message(
                paste0(
                    c("Setting high threshold for subtype ",subtype," to ",high_threshold_subtype),
                    collapse = ""
                )
            )
        } else {
            high_threshold_subtype = high_threshold
        }
        
        if (! (subtype %in% colnames(deconvolution_results)))
            next()
        vis_mat[
            (deconvolution_results[,subtype] <= high_threshold_subtype) |
            (deconvolution_results[,subtype] <= #mean(deconvolution_results[,subtype]))
            #(mean(deconvolution_results[,subtype])+ sd(deconvolution_results[,subtype]))
            quantile(
                deconvolution_results[,subtype],
                probs = seq(0,1,.01)
            )[high_threshold_subtype]),
            subtype
        ] = "low"
        vis_mat[
            (deconvolution_results[,subtype] > high_threshold_subtype) &
            (deconvolution_results[,subtype] > #mean(deconvolution_results[,subtype]))
            #(mean(deconvolution_results[,subtype])+ sd(deconvolution_results[,subtype]))
            quantile(
                deconvolution_results[,subtype],
                probs = seq(0,1,.01)
            )[high_threshold_subtype]),
            subtype
        ] = "high"
    }
    
    # p_value setting
    
    for(subtype in cands_dif){
        vis_mat[
            as.double(vis_mat$Confidence_score_dif) >=
            confidence_threshold,
            subtype
        ] = "not_significant"
    }
    
    for(subtype in c(cands_de_dif)){
        vis_mat[
            as.double(vis_mat$Confidence_score_de_dif) >=
                confidence_threshold,
            subtype
            ] = "not_significant"
    }

    # differentiation score scaling
    for ( score in c("Confidence_score_de_dif","Confidence_score_dif")){
        
        lower_index = which( vis_mat[ ,score] <= (mean(vis_mat[ ,score])-sd(vis_mat[ ,score])*.5) )
        vis_mat[lower_index,score] = mean(vis_mat[ ,score]) - (sd(vis_mat[ ,score])*.5)
        
        higher_index = which( vis_mat[ ,score] >= (mean(vis_mat[ ,score])+sd(vis_mat[ ,score])*.5) )
        vis_mat[higher_index,score] = mean(vis_mat[ ,score]) + (sd(vis_mat[ ,score])*.5)
    }
    vis_mat$Ratio = ( 
        #scale(vis_mat$Confidence_score_de_dif) / ( scale(vis_mat$Confidence_score_dif))
        scale(vis_mat$Confidence_score_de_dif / vis_mat$Confidence_score_dif)
    )
    vis_mat$Ratio_numeric = round(vis_mat$Ratio, 3)
    
    ### ratio adjustment
    
    #lower_index = which( vis_mat$Ratio <= (mean(vis_mat$Ratio)-sd(vis_mat$Ratio)*.5) )
    #higher_index = which( vis_mat$Ratio >= (mean(vis_mat$Ratio) + sd(vis_mat$Ratio)*.5) )
    
    lower_index = which(vis_mat$Ratio <= -.5)
    higher_index = which(vis_mat$Ratio >=.5)
    medium_index = 1:nrow(vis_mat)
    medium_index = medium_index[!(medium_index %in% c(lower_index,higher_index))]
    
    vis_mat$Ratio[lower_index] = "low"
    vis_mat$Ratio[medium_index] = "medium"
    vis_mat$Ratio[higher_index] = "high"
    
    #vis_mat$Ratio[lower_index] = mean(vis_mat$Ratio) - (sd(vis_mat$Ratio)*.5)
    #vis_mat$Ratio[higher_index] = mean(vis_mat$Ratio) + (sd(vis_mat$Ratio)*.5)
    
    #lower_index = which( vis_mat[,"MKI67"] <= (mean(vis_mat[,"MKI67"])-sd(vis_mat[,"MKI67"])*.5) )
    #higher_index = which( vis_mat[,"MKI67"] >= (mean(vis_mat[,"MKI67"]) + sd(vis_mat[,"MKI67"])*.5) )
    
    #vis_mat$MKI67[lower_index] = mean(vis_mat[,"MKI67"]) - (sd(vis_mat[,"MKI67"])*.5)
    #vis_mat$MKI67[higher_index] = mean(vis_mat[,"MKI67"]) + (sd(vis_mat[,"MKI67"])*.5)
    
    # case aggregated similarity
    if (aggregate_differentiated_stages){
        
        vis_mat[,"aggregated_similarity"] = rep("",nrow(vis_mat))
        for (i in 1:nrow(vis_mat)){

            candidate_list_max = vis_mat[i,cands_dif]
            candidate_list_index = which(candidate_list_max %in% c("low","high"))
            
            if (length(candidate_list_index) == 0){
                vis_mat[i,"aggregated_similarity"] = "undetermined"
                next()
            }
            
            high_types = colnames(candidate_list_max)[candidate_list_index]
            
            if (length(high_types) == 1){
                
                max_type = high_types
            } else {
                
                max_index = which.max(deconvolution_results[i,high_types])
                vis_mat[i,"aggregated_similarity"] = colnames(deconvolution_results[,high_types])[max_index]
            }
        }
        
        vis_mat[
            as.double(vis_mat$Confidence_score_dif) >=
                confidence_threshold,
            "aggregated_similarity"
        ] = "not_significant"

        vis_mat = vis_mat[,!(colnames(vis_mat) %in% cands_dif)]
        
        # rearrange order
        
        aggregated_similarity = vis_mat$aggregated_similarity
        vis_mat = vis_mat[,colnames(vis_mat) != "aggregated_similarity"]
        vis_mat = cbind(aggregated_similarity,vis_mat)
        
        vis_mat[vis_mat[,"aggregated_similarity"] == "","aggregated_similarity"] = "undetermined"
    }

    # add mki67 information
    if ("MKI67" %in% colnames(deconvolution_results)){
        
        vis_mat$MKI67 = rep("",nrow(vis_mat))
        mki_67 = as.double(deconvolution_results[,"MKI67"])
        quantiles = quantile(mki_67,probs = seq(0,1,.01))
        vis_mat[mki_67>=quantiles[67],"MKI67"] = "high"
        vis_mat[mki_67<quantiles[67],"MKI67"] = "medium"
        vis_mat[mki_67<=quantiles[34],"MKI67"] = "low"
    }
    
    # remove confidence scores
    vis_mat_filtered = vis_mat[,
        !(colnames(vis_mat) %in% c(
            "Confidence_score_de_dif",
            "Confidence_score_dif",
            "Ratio_numeric"
        ))
    ]

    # correlation heatmap
    pheatmap::pheatmap(
        correlation_matrix,
        annotation_col = vis_mat_filtered,
        annotation_colors = Graphics_parameters,
        annotation_legend = TRUE,
        treeheight_col = 0,
        treeheight_row = 0,
        show_colnames = show_colnames,
        show_rownames = FALSE
    )

    return(vis_mat)
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
        NEC_NET                 = c(NEC= "red", NET = "lightblue"),
        alpha        = c(not_significant = "gray", low = "white", high = "blue"),
        beta         = c(not_significant = "gray", low = "white", high = "green"),
        gamma        = c(not_significant = "gray", low = "white", high = "brown"),
        delta        = c(not_significant = "gray", low = "white", high = "Purple"),
        ductal     = c(not_significant = "gray", low = "white", high = "orange"),
        acinar     = c(not_significant = "gray", low = "white", high = "cyan"),
        #progenitor = c(not_significant = "gray", low = "white", high = "orange"),
        hisc   = c(not_significant = "gray",low = "white", high = "black"),
        hesc   = c(not_significant = "gray",low = "white", high = "black"),
        Confidence_score_de_dif = c(low = "green", medium = "white", high = "red"),
        Ratio = c(low = "darkgreen", medium = "yellow", high = "darkred"),
        MKI67 = c(low = "green", medium = "yellow", high = "red"),
        Confidence_score_dif = c(low = "green", medium = "white", high = "red"),
        Study = c(Groetzinger = "darkgreen", Scarpa = "darkred"),
        aggregated_similarity = c(
            alpha = "blue",
            beta = "green",
            gamma = "brown",
            delta = "purple",
            acinar = "cyan",
            ductal = "orange",
            not_significant = "gray",
            undetermined = "white"
        ),
        Differentiation_Stage_Aggregated = c(
            differentiated = "darkgreen",
            hisc = "red",
            #progenitor     = "orange",
            Not_significant= "gray"
        ),
        Grading          = c( G1 = "darkgreen",G2 = "Yellow", G3 = "darkred", G0 = "white")
    )
    return(Graphics_parameters)
}

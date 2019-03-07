create_visualization_matrix = function(
    deconvolution_results,
    confidence_threshold,
    high_threshold
){
    
    # data cleansing
    deconvolution_results$Subtype =
        str_to_lower(deconvolution_results$Subtype)

    # init variables
    
    subtype_cands = c("alpha","beta","gamma","delta","acinar","ductal","progenitor","hisc")
    cands_dif_endocrine = c("alpha","beta","gamma","delta")
    cands_dif_exokrine = c("ductal","acinar")
    
    # get the index for differentiated state
    
    model_1 = head(unlist(str_split(as.character(deconvolution_results$model),pattern = "\\|")),1)
    model_2 = head(unlist(str_split(as.character(deconvolution_results$model),pattern = "\\|")),2)[2]
    subtypes_1 = str_to_lower(as.character(unlist(str_split(model_1,pattern="_"))))
    subtypes_2 = str_to_lower(as.character(unlist(str_split(model_2,pattern="_"))))
    subtypes_1 = subtypes_1[str_to_lower(subtypes_1) %in% subtype_cands]
    subtypes_2 = subtypes_2[str_to_lower(subtypes_2) %in% subtype_cands]
    cands_dif_2 = subtypes_2[!(subtypes_2 %in% subtypes_1)]
    cands_dif_1 = subtypes_1[!(subtypes_1 %in% cands_dif_2)]
    
    cands_1_index = grep(
        colnames(deconvolution_results),
        pattern = paste0(
            str_to_lower(cands_dif_1),
            collapse = "|"
        )
    )
    
    cands_2_index = grep(
        colnames(deconvolution_results),
        pattern = paste0( 
            str_to_lower(cands_dif_2),
            collapse = "|"
        )
    )
    colnames(deconvolution_results)[c(cands_1_index,cands_2_index)] =
        str_to_lower(colnames(deconvolution_results)[c(cands_1_index,cands_2_index)])
    
    correlation_matrix = cor(transcriptome_mat_vis)
    deconvolution_results = deconvolution_results[
        colnames(correlation_matrix),
        ]

    # extract similarity measurements
    subtype_index_vis = c(
        cands_1_index,
        cands_2_index,
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
    
    # differentiation score scaling, zero offset
    for ( score in c("Confidence_score_de_dif","Confidence_score_dif")){
        
        vis_mat[,score] = as.double(vis_mat[,score])
        vis_mat[vis_mat[,score] == 0.0,score] = 10^-5
        vis_mat[,score] = vis_mat[,score] - min(vis_mat[,score])
        vis_mat[,score] = vis_mat[,score] / max(vis_mat[,score])
    }
    
    # heuristic for significance
    if ( confidence_threshold > 1.0 )
        confidence_threshold = 
        mean(vis_mat[,"Confidence_score_dif"]) +
        sd(vis_mat[,"Confidence_score_dif"])
    
    for(subtype in c(cands_dif_2,cands_dif_1)){
        
        subtype = str_to_lower(subtype)
        deconvolution_results[,subtype] = deconvolution_results[,subtype] - min(deconvolution_results[,subtype])
        
        if ( max( deconvolution_results[,subtype]) > 0 )
            deconvolution_results[,subtype] = 
            deconvolution_results[,subtype] /
            max(deconvolution_results[,subtype]
            )
        deconvolution_results[,subtype] = round(deconvolution_results[,subtype] * 100,1)
        
        # if no manual definition of high subtype similarity has been set
        if (high_threshold > 100){
            
            # heuristic   
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
        
        quantile_threshold = quantile(
            deconvolution_results[,subtype],
            probs = seq(0,1,.01)
        )[high_threshold_subtype]
        
        vis_mat[
            (deconvolution_results[,subtype] <= high_threshold_subtype) |
                (deconvolution_results[,subtype] <= quantile_threshold ),
            subtype
            ] = "low"
        
        vis_mat[
            (deconvolution_results[,subtype] > high_threshold_subtype) &
                (deconvolution_results[,subtype] > quantile_threshold),
            subtype
            ] = "high"
        
        cands_dif_endocrine = cands_dif_endocrine[
            cands_dif_endocrine %in% colnames(deconvolution_results)
            ]
        cands_dif_exokrine = cands_dif_exokrine[
            cands_dif_exokrine %in% colnames(deconvolution_results)
            ]
        
        # identify the most similary subtype i.e. 'subtype' column
        
        if( subtype %in% cands_dif_endocrine){
            
            max_indices = as.integer(apply( 
                deconvolution_results[,cands_dif_endocrine],
                MARGIN = 1,
                FUN = which.max
            ))
            max_subtypes = colnames(deconvolution_results[,cands_dif_endocrine])[max_indices]
            highest_percentage = max_subtypes == subtype
            
        } else if ( subtype %in% cands_dif_exokrine){
            
            max_indices = as.integer(apply( 
                deconvolution_results[,cands_dif_exokrine],
                MARGIN = 1,
                FUN = which.max
            ))
            max_subtypes = colnames(deconvolution_results[,cands_dif_exokrine])[max_indices]
            highest_percentage = max_subtypes == subtype 
            
        } else {
            highest_percentage = rep(TRUE,nrow(deconvolution_results))
        }
        
        vis_mat[
            highest_percentage != TRUE,
            subtype
            ] = "low"
    }
    
    # p_value setting
    
    for(subtype in cands_dif_1){
        subtype = str_to_lower(subtype)
        vis_mat[
            as.double(vis_mat$Confidence_score_dif) >=
                confidence_threshold,
            subtype
            ] = "not_significant"
    }
    
    for(subtype in c(cands_dif_2)){
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
    vis_mat$Confidence_score_dif[vis_mat$Confidence_score_dif == 0] = 1*10^-5
    vis_mat$Confidence_score_de_dif[vis_mat$Confidence_score_de_dif == 0] = 1*10^-5
    vis_mat$Ratio = ( 
        scale(vis_mat$Confidence_score_de_dif / vis_mat$Confidence_score_dif)
    )
    vis_mat$Ratio_numeric = round(vis_mat$Ratio, 3)
    
    ### ratio adjustment
    
    lower_index = which(vis_mat$Ratio <= -.5)
    higher_index = which(vis_mat$Ratio >=.5)
    medium_index = 1:nrow(vis_mat)
    medium_index = medium_index[!(medium_index %in% c(lower_index,higher_index))]
    
    vis_mat$Ratio[lower_index] = "low"
    vis_mat$Ratio[medium_index] = "medium"
    vis_mat$Ratio[higher_index] = "high"
    
    # calculate aggregated similarity
    
    ###vis_mat[,"aggregated_similarity"] = deconvolution_results$Subtype
    
    # add mki67 information
    if ("MKI67" %in% colnames(deconvolution_results)){
        
        vis_mat$MKI67 = rep("",nrow(vis_mat))
        mki_67 = as.double(deconvolution_results[,"MKI67"])
        quantiles = quantile(mki_67,probs = seq(0,1,.01))
        vis_mat[mki_67>=quantiles[67],"MKI67"] = "high"
        vis_mat[mki_67<quantiles[67],"MKI67"] = "medium"
        vis_mat[mki_67<=quantiles[34],"MKI67"] = "low"
    }

    return(vis_mat)
}


#' create_PCA_differentiation_stages
#'
#' \code{create_PCA_differentiation_stages}
#' visualizes the differentiation stage predictions as PCA.
#' Please note that the first column of the expression
#' data matrix has to contain the HGNC identifier
#'
#' @param visualization_data Transcriptome data that shall be visualized.
#' Notice the convention that the first row has to contain the HGNC identifier
#' @param deconvolution_results The dataframe returned
#' by the deconvolution analysis
#' @param Graphics_parameters Pheatmap visualization paramters.
#' You can customize visualization colors.
#' Read the vignette for more information.
#' @usage
#' create_PCA_differentiation_stages(
#'     visualization_data_path,
#'     deconvolution_results,
#'     annotation_columns,
#'     Graphics_parameters,
#'     baseline
#' )
#' @examples
#' data(deconvolution_results)
#' data(visualization_data)
#' create_PCA_differentiation_stages(
#'     visualization_data = visualization_data,
#'     deconvolution_results = deconvolution_results,
#'     Graphics_parameters = ""
#' )
#' @import stringr ggplot2 pheatmap ggfortify
#' @export
create_PCA_differentiation_stages = function(
    visualization_data,
    deconvolution_results,
    aggregate_differentiated_stages = FALSE,
    confidence_threshold = 1.1,
    show_colnames = FALSE,
    Graphics_parameters = "",
    high_threshold = 101
){
        
    vis_mat = create_visualization_matrix(
        deconvolution_results = deconvolution_results,
        confidence_threshold = confidence_threshold,
        high_threshold = high_threshold
    )

    colnames(transcriptome_mat_vis) = str_replace_all(
        colnames(transcriptome_mat_vis),
        pattern = "^X",
        ""
    )
    
    # init variables
    
    subtype_cands = c("alpha","beta","gamma","delta","acinar","ductal","progenitor","hisc")
    cands_dif_endocrine = c("alpha","beta","gamma","delta")
    cands_dif_exokrine = c("ductal","acinar")
    cands_de_dif = "hisc"
    
    correlation_matrix = cor(visualization_data)
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
        pattern = paste0(c(
            cands_dif_endocrine,
            cands_dif_exokrine)
        , collapse = "|")
    )
    
    de_dif_index = grep(
        colnames(deconvolution_results),
        pattern = "hisc"
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
    
    #vis_mat_filtered = deconvolution_results[subtype_index_vis]

    
    p = ggbiplot::ggbiplot(
        pcr,
        obs.scale = .75,
        groups = vis_mat$Grading,
        ellipse = TRUE,
        circle = TRUE,
        var.axes = F#,labels = meta_data$Name
    )
    p = p + geom_point(
        aes(
            colour = as.character(vis_mat$Ratio),
            shape = vis_mat$Ratio#,
            #size = 5
        )
    )
    p = p + guides(
        color=guide_legend(
            title="Ratio"
        ),
        #size = guide_legend(title="Ratio"),
        shape =  guide_legend(title="Grading")
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
#' @param confidence_threshold Threshold above which deconvolutions are deemed unsuccessful
#' and corresponding results being masked on the the plots
#' @param show_colnames Whether to show the sample column names
#' @param Graphics_parameters Pheatmap visualization paramters.
#' You can customize visualization colors.
#' Read the vignette for more information.
#' @param high_threshold Threshold depending on which a deconvolution result
#' is interpreted as 'high'. If not set, a statistical estimation will approximately
#' identify a signficance threshold for a high similarity.
#' @usage
#' create_heatmap_differentiation_stages(
#'     visualization_data,
#'     deconvolution_results,
#'     aggregate_differentiated_stages,
#'     confidence_threshold,
#'     show_colnames,
#'     Graphics_parameters,
#'     high_threshold
#' )
#' @examples
#' data(visualization_data)
#' data(deconvolution_results)
#' data(meta_data)
#'
#' create_heatmap_differentiation_stages(
#'     visualization_data = visualization_data,
#'     deconvolution_results = deconvolution_results,
#'     aggregate_differentiated_stages = FALSE,
#'     confidence_threshold = 1.1,
#'     show_colnames = FALSE,
#'     Graphics_parameters = "",
#'     high_threshold = 101
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

    vis_mat = create_visualization_matrix(
        deconvolution_results = deconvolution_results,
        confidence_threshold = confidence_threshold,
        high_threshold = high_threshold
    )
    
    transcriptome_mat_vis = visualization_data
    colnames(transcriptome_mat_vis) = str_replace_all(
        colnames(transcriptome_mat_vis),
        pattern = "^X",
        ""
    )
    
    # assert that graphics parameters are available
    if ( head(Graphics_parameters,1) == "" )
        Graphics_parameters = configure_graphics()
    
    # remove confidence scores
    vis_mat_filtered = vis_mat[,
        !(colnames(vis_mat) %in% c(
            "Confidence_score_de_dif",
            "Confidence_score_dif",
            "Ratio_numeric"
        ))
    ]
    if (aggregate_differentiated_stages){
        
        max_indices = apply(deconvolution_results[,
                                                  c(cands_dif_endocrine,cands_dif_exokrine)
                                                  ],FUN=which.max,MARGIN = 1)
        subtype_selection = vis_mat[,c(cands_dif_endocrine,cands_dif_exokrine)]
        vis_mat$Aggregated_similarity = rep("",nrow(vis_mat))
        
        for( j in 1:nrow(subtype_selection)){
            
            subtype_strength = subtype_selection[j,max_indices[j]]
            if (subtype_strength == "low"){
                vis_mat$Aggregated_similarity[j] = "not_significant"
            } else {
                vis_mat$Aggregated_similarity[j] = colnames(subtype_selection)[max_indices[j]]
            }
        }
        
        column_candidates = c("Aggregated_similarity","hisc","Grading","MKI67")
        column_candidates = column_candidates[column_candidates %in% colnames(vis_mat)]
        
        vis_mat[
            as.double(vis_mat$Confidence_score_dif) >=
                confidence_threshold,
            "aggregated_similarity"
            ] = "not_significant"
        
        vis_mat_filtered = vis_mat[,column_candidates]
        
    }

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
        Aggregated_similarity = c(
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

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
#' @param annotation_columns Column names that are to be visualized
#' ontop of the correlation matrix
#' @param Graphics_parameters Pheatmap visualization paramters.
#' You can customize visualization colors.
#' Read the vignette for more information.
#' @param baseline Which measurement represents the baseline
#' of the differentiation similarity: 'absolute' = maximal
#' similarity of the training sample to its subtype. 'relative' 
#' sets the baseline to maximal similarity of the test samples
#' currently analysed
#' @usage
#' create_heatmap_differentiation_stages(
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
#' meta_data      = read.table(
#'     meta_data_path, sep ="\t",
#'     header = TRUE,
#'     stringsAsFactors = FALSE
#' )
#' rownames(meta_data) = meta_data$Sample
#'
#' transcriptome_file_path = system.file(
#'     "/Data/Expression_data/Visualization_PANnen.tsv",
#'      package = "artdeco"
#' )
#' create_heatmap_differentiation_stages(
#'     transcriptome_file_path = transcriptome_file_path,
#'     deconvolution_results = meta_data,
#'     annotation_columns = c(
#'         "Differentiation_Stages_Subtypes",
#'         "Differentiation_Stages_Aggregated",
#'         "Differentiatedness"
#'      ),
#'      Graphics_parameters = "",
#'      baseline = "absolute"
#' )
#' @import stringr ggplot2 pheatmap ggfortify
#' @return Plots
#' @export
create_heatmap_differentiation_stages = function(
    transcriptome_file_path,
    deconvolution_results,
    annotation_columns = c(
        "Differentiation_Stages_Subtypes",
        "Differentiation_Stages_Aggregated",
        "Differentiatedness"
    ),
    Graphics_parameters = "",
    baseline = "absolute"
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

    annotation_data$Differentiation_Stages_Subtypes[
        annotation_data$Differentiation_Stages_Subtypes %in% c("hesc","hisc")
    ] = "stem_cell"

    if (baseline == "relative"){

        # procure subtype information

        sub_index = grep(
            x = colnames(deconvolution_results),
            pattern = paste0( c(
                "(_similarity$)"),
                collapse = "|"
            )
        )
        subtypes = colnames(deconvolution_results)[ sub_index]
        subtypes[
            str_replace_all(
                subtypes,
                "_similarity$",
                ""
            )  %in% c("hesc","hisc")
        ] = "stem_cell_similarity"

        # load the model information

        models = as.character(
            unlist(
                str_split(meta_data$model[1],pattern = "\\|")
            )
        )
        # readjust similarity predictions

        for ( model in models){

            model_path = paste0(c(
                system.file(
                    "Models/",
                    package="artdeco"),
                "/",
                model,".RDS"
            ),
            collapse = ""
            )
            model_and_parameter = readRDS(model_path)
            fit = model_and_parameter[1]
            fit = fit[[1]]
            parameter_list = model_and_parameter[2]
            model_parameter = model_parameter[[1]]
            not_sig_samples = rownames(fit)[fit["P-value"] >= 0.05]

            annotation_data = quantify_similarity(
                meta_data = deconvolution_results,
                subtypes = subtypes,
                model = model,
                fit = fit,
                parameter_list,
                baseline = baseline,
                not_sig_samples = not_sig_samples
            )
        }
    }

    # correlation heatmap
    pheatmap::pheatmap(
        correlation_matrix,
        annotation_col = annotation_data,
        annotation_colors = Graphics_parameters,
        annotation_legend = TRUE,
        treeheight_col = 0,
        treeheight_row = 0,
        show_colnames = TRUE,
        show_rownames = FALSE
    )

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
        alpha_similarity        = c(not_significant = "gray", none = "white", traces = "yellow", significant = "blue"),
        beta_similarity         = c(not_significant = "gray", none = "white", traces = "yellow", significant = "darkgreen"),
        gamma_similarity        = c(not_significant = "gray", none = "white", traces = "yellow", significant = "brown"),
        delta_similarity        = c(not_significant = "gray", none = "white", traces = "yellow", significant = "Purple"),
        ductal_similarity     = c(not_significant = "gray", none = "white", traces = "yellow", significant = "Black"),
        acinar_similarity     = c(not_significant = "gray", none = "white", traces = "yellow", significant = "Brown"),
        progenitor_simimilarity = c(not_significant = "gray", none = "white", traces = "yellow", significant = "orange"),
        stem_cell_similaritry   = c(none = "white", traces = "yellow", significant = "darkred",not_significant = "gray"),
        Differentiatedness      = c(none = "white", traces = "yellow", significant = "darkgreen"),
        Differentiation_Stages_Aggregated = c(
            differentiated   = "darkgreen",
            dedifferentiated = "darkred",
            not_significant  = "gray"
        ),
        Differentiation_Stages_Subtypes = c(
            alpha = "blue",
            beta = "green",
            gamma = "brown",
            delta = "purple",
            acinar = "cyan",
            ductal = "red",
            stem_cell      = "black",
            progenitor     = "orange",
            not_significant= "gray"
        ),
        Grading          = c( G1 = "Green",G2 = "Yellow", G3 = "Red", G0 = "white")
    )
    return(Graphics_parameters)
}

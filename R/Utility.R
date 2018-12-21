#' remove_model
#'
#' \code{remove_model} removes a model
#'
#' @param model_name Name of the model
#' @param test_mode Testrun indicator
#' @usage
#' remove_model(
#'     model_name,
#'     test_mode
#' )
#' @examples
#' model = list.files(system.file("Models/", package = "artdeco"))[1]
#' remove_model(
#'     model_name = model,
#'     test_mode = TRUE
#' )
#' @import stringr
#' @export
#' @return list of models
remove_model = function(
    model_name,
    test_mode = FALSE
){

    # procure path
    model_name = str_replace_all(
        model_name,
        pattern = "(.rds)|(.RDS)",
        ""
    )

    model_path = paste0(
        c(system.file("Models/", package = "artdeco"),
            "/", paste0(model_name, ".RDS")
        ), sep ="", collapse = ""
    )

    # check if model exists
    if (!file.exists(model_path)){
        stop( paste0(c(
            "Cannot delete model ", model_name, ", model not found."
            ),collapse = "" )
        )
    }else{
        if (! test_mode)
            file.remove(model_path)
    }
    print(paste0("Deleted model: ",model_name))
}

#' show_models
#'
#' \code{remove_model} shows which models have been trained
#'
#'@usage
#' show_models()
#'@examples
#' show_models()
#'@import stringr
#'@export
#'@return list of models
show_models = function(
){

    package_path = system.file("", package = "artdeco")
    model_path = paste(package_path,"Models", sep = "/")

    models = stringr::str_replace_all(
        list.files(model_path),
        pattern = ".RDS",
        ""
    )

    return(models)
}

#' prepare_result_matrix
#'
#' Prepares the deconvolution results
#' @param prediction_stats_list list containing the fitting results
#' @param parameter_list List containing the model parameters
#' @param p_value P-value significance threshold
#' @param meta_data Contains the deconvolution results as dataframe
#' @param baseline Specifying the baseline mode
#' @param models List of training models
#' @usage
#' prepare_result_matrix(
#'     prediction_stats_list,
#'     parameter_list,
#'     p_value,
#'     meta_data,
#'     baseline,
#'     models
#' )
#' @import stringr
#' @return list of models
prepare_result_matrix = function(
    prediction_stats_list,
    parameter_list,
    p_value,
    meta_data,
    baseline,
    models
){

    for (model in models){ # for all models

        fit = as.data.frame( prediction_stats_list[[model]] )
        colnames(fit) = str_replace_all(
            colnames(fit),
            pattern = paste0(
                c(eval(model),"."),collapse = ""
            ),
            ""
        )
        fit_parameter = parameter_list[[model]]

        # last four entries are always p-value, correlation, RMSE and score
        subtypes = str_to_lower(colnames(fit)[1:(ncol(fit)-4)])
        res_coeff = fit[,subtypes]
        res_cor   = fit[!(colnames(fit) %in% subtypes)]

        res_coeff[ is.na(res_coeff) ] = 0.0
        res_cor[ is.na(res_cor) ]     = 0.0

        not_sig_samples = rownames(res_cor)[
            which(res_cor[,"P-value"] > p_value)]
        not_sig_samples

        meta_data = quantify_similarity(
            meta_data = meta_data,
            subtypes  = subtypes,
            fit       = fit,
            parameter_list = parameter_list,
            model     = model,
            baseline  = baseline,
            not_sig_samples   = not_sig_samples
        )

        # Calculate aggregated proportions
    }

    # find highest similarity

    max_mat  = meta_data[,grep(colnames(meta_data),pattern = "_similarity_percent")]
    subtypes = str_replace_all(
        grep(colnames(meta_data), pattern = "_similarity$",value = TRUE),
        pattern = "_similarity",
        ""
    )
    maxi = apply(  max_mat, FUN = which.max, MARGIN = 1 )
    meta_data$Differentiation_Stages_Subtypes = subtypes[maxi]

    for (i in 1:nrow(meta_data)){

        signess = (meta_data[i,
            paste(meta_data$Differentiation_Stages_Subtypes[i],"similarity", sep ="_")
        ])

        if( signess %in% c("not_significant","none","traces") )
            meta_data[i,"Differentiation_Stages_Subtypes"] = "not_significant"
    }

    # add differentiation information

    differentiated_index = grep(
        pattern = paste0(c(
            "(alpha_similarity_percent)",
            "(beta_similarity_percent)",
            "(gamma_similarity_percent)",
            "(delta_similarity_percent)",
            "(accinar_similarity_percent)",
            "(ductal_similarity_percent)"),
            collapse = "|"),
        colnames(meta_data),
        value = FALSE
    )
    de_differentiated_index = grep(
        pattern = paste0(c(
            "(progenitor_similarity_percent)",
            "(hsc_similarity_percent)"),
            collapse = "|"),
        colnames(meta_data),
        value = FALSE
    )

    meta_data[,"Differentiatedness"] = rep("",nrow(meta_data))
    meta_data[,"Differentiation_Stages_Aggregated"] = rep("",nrow(meta_data))

    # quantify differentiation stages

    if (length(differentiated_index) > 0 ){ # if we have subtypes that quantify diff.

        for (i in 1:nrow(meta_data)){ # for all samples

            differentiatedness = round(sum(meta_data[i,differentiated_index]),0)
            meta_data[i,"Differentiatedness"] = differentiatedness

            # branch we have at least one dedifferentiated subtype

            if( length(de_differentiated_index) > 0 ){

                dedifferentiatedness = sum(meta_data[i,de_differentiated_index])
                differentiatedness = round(
                    differentiatedness / de_differentiated_index
                , 0 )
                meta_data[i,"Differentiatedness"] = differentiatedness # update

                if(dedifferentiatedness > 50) {

                    meta_data[i,"Differentiation_Stages_Aggregated"] = "dedifferentiated"

                } else if (differentiatedness <= 50){

                    meta_data[i,"Differentiation_Stages_Aggregated"] = "not_significant"
                }

            } else { # branch we do not have any dedifferentiation information

                if (differentiatedness > 50) {
                    meta_data[i,"Differentiation_Stages_Aggregated"] = "differentiated"
                } else {
                    meta_data[i,"Differentiation_Stages_Aggregated"] = "not_significant"
                }
            }
        }
    }
    return(meta_data)
}

#' quantify_similarity
#'
#' \code{quantify_similarity} custom function that allows to
#' alternate between different baselines for similarity
#'
#'@param meta_data deconvolution result dataframe
#'@param subtypes list of cell subtypes
#'@param model Deconvolution model
#'@param fit Deconvolution fit
#'@param parameter_list List of parameters of the models
#'@param baseline Absolut or relative mode
#'@param not_sig_samples list of not significantly deconvolable samples
#'@usage
#' quantify_similarity(
#'     meta_data,
#'     subtypes,
#'     model,
#'     fit,
#'     parameter_list,
#'     baseline,
#'     not_sig_samples
#' )
#'@import stringr
#'@export
#'@return dataframe containing baseline adjusted similarities
quantify_similarity = function(
    meta_data,
    subtypes,
    model,
    fit,
    parameter_list,
    baseline,
    not_sig_samples
){

    subtypes = str_replace_all(
        subtypes,
        pattern = "_similarity",
        ""
    )
    subtypes[subtypes %in% c("hesc","hisc")] = "stem_cell"
    subtypes = paste0(
        subtypes,
        "_similarity"
    )

    sim_index = grep(colnames(meta_data),pattern = "_similarity$", value = F)
    colnames(meta_data)[sim_index] = str_replace_all(
        colnames(meta_data)[sim_index],
        pattern = "_similarity$",
        ""
    )
    colnames(meta_data)[sim_index] = paste(
        colnames(meta_data)[sim_index],
        "similarity",
        sep ="_"
    )

    colnames = colnames(fit)
    rownames = rownames(fit)
    fit = matrix(
        as.double(
            as.character(
                unlist(fit)
            )
        ),
        nrow = length(rownames),
        ncol = length(colnames)
    )

    colnames = str_replace_all(
        colnames,
        pattern = "_similarity$",
        ""
    )
    colnames[ colnames %in% c("hesc","hisc") ] = "stem_cell"
    colnames = paste(
        colnames,
        "similarity",
        sep ="_"
    )

    colnames(fit) = colnames
    rownames(fit) = rownames

    res_coeff = data.frame(fit)[,subtypes]

    for (subtype in subtypes){

        subtype = str_replace_all(
            subtype,
            pattern = "_similarity",
            ""
        )

        subtype_label = paste0(
            subtype,
            "_similarity"
        )

        subtype_label_quant = paste0(
            c(
                subtype,
                "_similarity_percent"
            ),
            collapse =""
        )

        # take the measurement
        colnames(res_coeff) = str_replace_all(
            colnames(res_coeff),
            pattern = "_similarity",
            ""
        )
        subtype_sim_scalar = log( res_coeff[,subtype] + 1 )

        # all quantifications are initialized as empty
        meta_data[,subtype_label] =
            rep("",length(subtype_sim_scalar))

        # decide on what to use as baseline
        if (baseline == "absolute"){

            # obtain measurement from training phase of model
            parameter = parameter_list[[model]]
            parameter = parameter[[1]]
            names(parameter)[names(parameter) %in% c("hesc","hisc")] = "stem_cell"

            parameter = as.double(
                as.character(
                    unlist(
                        parameter[subtype]
                    )
                )
            )
            sub_fit_max = max(parameter)
            sub_fit_max = log(sub_fit_max + 1)
            sub_fit_max = max(sub_fit_max, max(subtype_sim_scalar))

            similarity_quantified = round(subtype_sim_scalar / sub_fit_max * 100,0)

            significant_index   = which( similarity_quantified <= 100 )
            traces_index        = which( similarity_quantified <= 25 )
            none_index          = which( similarity_quantified <= 10 )

        } else { # use use-case specific measurements

            similarity_quantified = round(
                subtype_sim_scalar / max( subtype_sim_scalar ) * 100,
                0
            )

            quants = quantile(
                similarity_quantified,
                seq(0,1, by =.1)
            )

            significant_index   = which( similarity_quantified <= quants[10] )
            traces_index        = which( similarity_quantified <= quants[6] )
            none_index          = which( similarity_quantified <= quants[2] )
        }

        not_sig_index           = which( rownames(meta_data) %in% not_sig_samples )

        # shape measurements relative to baseline
        meta_data[,subtype_label_quant]  = similarity_quantified

        meta_data[ significant_index, subtype_label ] = "significant"
        meta_data[ traces_index     , subtype_label ] = "traces"
        meta_data[ none_index       , subtype_label ] = "none"
        meta_data[ not_sig_index    , subtype_label ] = "not_significant"
        meta_data[
            (meta_data[ ,
                subtype_label ] == "") | ( is.na(meta_data[ , subtype_label
            ]) ),
            subtype_label
        ] = "none"
    }
    return(meta_data)
}

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
#' @param meta_data Contains the deconvolution results as dataframe
#' @param models List of training models
#' @usage
#' prepare_result_matrix(
#'     prediction_stats_list,
#'     parameter_list,
#'     meta_data,
#'     baseline,
#'     models
#' )
#' @import stringr
#' @return list of models
prepare_result_matrix = function(
    prediction_stats_list,
    parameter_list,
    meta_data,
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

        meta_data = quantify_similarity(
            meta_data = meta_data,
            subtypes  = subtypes,
            fit       = fit,
            model     =  model,
            parameter_list = parameter_list
        )

        # Calculate aggregated proportions
    }

    # calculate differentiatedness
    
    meta_data[,"Differentiatedness_absolute_log_odds"] = rep(0,nrow(meta_data))
    meta_data[,"Differentiatedness_relative_log_odds"] = rep(0,nrow(meta_data))
    meta_data[,"Differentiatedness_absolute"] = rep("",nrow(meta_data))
    meta_data[,"Differentiatedness_relative"] = rep("",nrow(meta_data))
    
    differentiatedness_relative =
    rowSums(
        meta_data[,
            grep( colnames(meta_data), pattern = 
                "(alpha_similarity_relative_percent)|(beta_similarity_relative_percent)|(gamma_similarity_relative_percent)|(delta_similarity_relative_percent)|(accinar_similarity_relative_percent)|(ductal_similarity_relative_percent)",
                value = T
            )
        ]
    )
    differentiatedness_absolute =
    rowSums(
        meta_data[,
                  grep( colnames(meta_data), pattern = 
                        "(alpha_similarity_absolute_percent)|(beta_similarity_absolute_percent)|(gamma_similarity_absolute_percent)|(delta_similarity_absolute_percent)|(accinar_similarity_absolute_percent)|(ductal_similarity_absolute_percent)",
                        value = T
                  )
                  ]
    )
    
    de_differentiatedness_relative =
        rowSums(
            meta_data[,
                grep( colnames(meta_data), pattern = 
                        "(progenitor_similarity_relative_percent)|(stem_cell_similarity_relative_percent)",
                    value = T
                )
            ]
        )
    de_differentiatedness_absolute =
        rowSums(
            meta_data[,
                grep( colnames(meta_data), pattern = 
                        "(progenitor_similarity_relative_percent)|(stem_cell_similarity_relative_percent)",
                    value = T
                )
            ]
        )
    
    meta_data[,"Differentiatedness_absolute_log_odds"] = log( 
        (differentiatedness_absolute + 1) / ( de_differentiatedness_absolute + 1)
    )
    meta_data[,"Differentiatedness_relative_log_odds"] = log(
        (differentiatedness_relative + 1) / ( de_differentiatedness_relative + 1)
    )
    
    meta_data$Differentiatedness_absolute = "Not_differentiated"
    meta_data$Differentiatedness_absolute[ 
        (-1 >= meta_data$Differentiatedness_absolute_log_odds) &
        (      meta_data$Differentiatedness_absolute_log_odds <= 1)
    ] = "Unclear"
    meta_data$Differentiatedness_absolute[
        meta_data[,"Differentiatedness_absolute_log_odds"] > 1] = "Differentiated"
    
    meta_data$Differentiatedness_relative = "Not_differentiated"
    meta_data$Differentiatedness_relative[ 
        (-1 >= meta_data$Differentiatedness_relative_log_odds) &
            (      meta_data$Differentiatedness_relative_log_odds <= 1)
        ] = "Unclear"
    meta_data$Differentiatedness_relative[
        meta_data[,"Differentiatedness_relative_log_odds"] > 1] = "Differentiated"

    # find highest similarity
    
    max_mat_relative  = meta_data[,grep(colnames(meta_data),pattern = "_similarity_relative_percent")]
    max_mat_absolute  = meta_data[,grep(colnames(meta_data),pattern = "_similarity_absolute_percent")]
    
    subtypes = unique(str_replace_all(
        grep(colnames(meta_data), pattern = "_similarity_absolute$",value = TRUE),
        pattern = "_similarity_absolute",
        ""
    ))
    maxi_relative = as.integer(apply(  max_mat_relative, FUN = which.max, MARGIN = 1 ))
    maxi_absolute = as.integer(apply(  max_mat_absolute, FUN = which.max, MARGIN = 1 ))
    meta_data$Differentiation_Stage_relative = subtypes[maxi_relative]
    meta_data$Differentiation_Stage_absolute = subtypes[maxi_absolute]
    
    meta_data[,"P_value"] = res_cor[,"P-value"]
    meta_data[,"Sig_score"] = res_cor[,ncol(res_cor)]
 
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
#'@usage
#' quantify_similarity(
#'     meta_data,
#'     subtypes,
#'     model,
#'     fit,
#'     parameter_list,
#' )
#'@import stringr
#'@export
#'@return dataframe containing baseline adjusted similarities
quantify_similarity = function(
    meta_data,
    subtypes,
    model,
    fit,
    parameter_list
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

    sim_index = grep(
        colnames(meta_data),
        pattern = "_similarity$",
        value = FALSE
    )
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

        # initiate labels for absolute and relative case
        subtype = str_replace_all(
            subtype,
            pattern = "_similarity",
            ""
        )

        subtype_label_absolute = paste0(
            subtype,
            "_similarity_absolute"
        )
        subtype_label_relative = paste0(
            subtype,
            "_similarity_relative"
        )

        subtype_label_quant_absolute = paste0(
            c(
                subtype,
                "_similarity_absolute_percent"
            ),
            collapse =""
        )
        subtype_label_quant_relative = paste0(
            c(
                subtype,
                "_similarity_relative_percent"
            ),
            collapse =""
        )

        # take the measurement
        colnames(res_coeff) = str_replace_all(
            colnames(res_coeff),
            pattern = "_similarity",
            ""
        )
        subtype_sim_scalar = res_coeff[,subtype]

        # all quantifications are initialized as empty
        meta_data[,subtype_label_absolute] =
            rep("",length(meta_data$Sample))
        meta_data[,subtype_label_relative] =
            rep("",length(meta_data$Sample))

        # case baseline is absolute

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
            maximum = max(parameter)
            #sub_fit_max = log(sub_fit_max + 1)
            #maximum_val = max(maximum, max(subtype_sim_scalar))

            similarity_quantified_absolute = round(
                exp( subtype_sim_scalar )  /
                    maximum * 100,
                1
            )
            similarity_quantified_absolute = similarity_quantified_absolute /
                max(similarity_quantified_absolute) * 100
            
            # load data

            none_index_absolute          = which( similarity_quantified_absolute >= 0.0 )
            traces_index_absolute        = which( similarity_quantified_absolute >= 10.0 )
            significant_index_absolute   = which( similarity_quantified_absolute >= 20.0 )
            
            meta_data[,subtype_label_quant_absolute] = as.double(similarity_quantified_absolute)
            meta_data[,subtype_label_absolute] = rep("",nrow(meta_data))
            
            meta_data[ none_index_absolute       , subtype_label_absolute ] = "none"
            meta_data[ traces_index_absolute     , subtype_label_absolute ] = "traces"
            meta_data[ significant_index_absolute, subtype_label_absolute ] = "present"
            
            meta_data[
                (meta_data[ ,
                            subtype_label_absolute ] == "") | ( is.na(meta_data[ , subtype_label_absolute
                                                                                 ]) ),
                subtype_label_absolute
            ] = "none"

        #case baseline relative

            similarity_quantified_relative = round(
                subtype_sim_scalar / max( subtype_sim_scalar ) * 100,
                0
            )
            #similarity_quantified_relative = log(similarity_quantified_relative+1)

            quants = quantile(
                similarity_quantified_relative,
                seq(0,1, by =.01)
            )

            # load data
            
            none_index_relative          = which( similarity_quantified_relative >= 0 )
            traces_index_relative        = which( similarity_quantified_relative >= quants[33]  )
            significant_index_relative   = which( similarity_quantified_relative >= quants[67] )

            meta_data[,subtype_label_quant_relative] = as.double(similarity_quantified_relative)
            meta_data[,subtype_label_relative] = rep("",nrow(meta_data))
            
            meta_data[ none_index_relative       , subtype_label_relative ] = "none"
            meta_data[ traces_index_relative     , subtype_label_relative ] = "traces"
            meta_data[ significant_index_relative, subtype_label_relative ] = "present"
            
            meta_data[
                (meta_data[ ,
                            subtype_label_relative ] == "") | ( is.na(meta_data[ , subtype_label_relative
                                                                                 ]) ),
                subtype_label_relative
            ] = "none"

    }

    return(meta_data)
}

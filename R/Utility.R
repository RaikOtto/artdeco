#' remove_model
#'
#' \code{remove_model} removes a model from a library
#'
#' @param model_name Name of the model
#' @param lib_name Name of the library that contains the model (e.g. "NMF", "music" or "bseqsc")
#' @param test_mode Testrun indicator
#' @usage
#' remove_model(
#'     model_name,
#'     lib_name,
#'     test_mode
#' )
#' @examples
#' remove_model(
#'     model_name = "My_model",
#'     lib_name = "bseqsc",
#'     test_mode = TRUE
#' )
#' @import stringr
#' @export
remove_model = function(
    model_name,
    lib_name = "NMF",
    test_mode = FALSE
){
    
    # procure path
    model_name = str_replace_all(
        model_name,
        pattern = "(.rds)|(.RDS)",
        ""
    )
    
    model_path <- paste0(system.file("Models", package = "artdeco"), 
                         "/", lib_name, "/", model_name, ".RDS")
    #model_path = paste0(
    #    c(system.file("Models/bseqsc", package = "artdeco"),
    #      "/", paste0(model_name, ".RDS")
    #    ), sep ="", collapse = ""
    #)
    
    # check if model exists
    if ( 
        file.exists(model_path)
    ){
        file.remove(model_path)
    }else if ( test_mode){
        print("Test mode active")
    }else{
        stop(
            paste0(
                c(
                    "Cannot delete model ",
                    model_name,
                    ", model not found. Check the library name."
                ),
                collapse = ""
            )
        )    
    }
    print(paste0("Deleted model ", model_name, "of library ", lib_name, "."))
}

#' remove_model_bseqsc
#'
#' \code{remove_model_bseqsc} removes a model from the bseqsc library
#'
#' @param model_name Name of the model
#' @param test_mode Testrun indicator
#' @usage
#' remove_model_bseqsc(
#'     model_name,
#'     test_mode
#' )
#' @examples
#' remove_model_bseqsc(
#'     model_name = "My_model",
#'     test_mode = TRUE
#' )
#' @import stringr
#' @export
remove_model_bseqsc = function(
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
        c(system.file("Models/bseqsc", package = "artdeco"),
          "/", paste0(model_name, ".RDS")
        ), sep ="", collapse = ""
    )
    
    # check if model exists
    if ( 
        file.exists(model_path)
    ){
        file.remove(model_path)
    }else if ( test_mode){
        print("Test mode active")
    }else{
        stop(
            paste0(
                c(
                    "Cannot delete model ",
                    model_name,
                    ", model not found."
                ),
                collapse = ""
            )
        )    
    }
    print(paste0("Deleted model: ",model_name))
}

#' remove_model_NMF
#'
#' \code{remove_model_NMF} removes a model from the NMF library
#'
#' @param model_name Name of the model
#' @param test_mode Testrun indicator
#' @usage
#' remove_model_NMF(
#'     model_name,
#'     test_mode
#' )
#' @examples
#' remove_model_NMF(
#'     model_name = "My_model",
#'     test_mode = TRUE
#' )
#' @import stringr
#' @export
remove_model_NMF = function(
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
        c(system.file("Models/NMF", package = "artdeco"),
          "/", paste0(model_name, ".RDS")
        ), sep ="", collapse = ""
    )
    
    # check if model exists
    if ( 
        file.exists(model_path)
    ){
        file.remove(model_path)
    }else if ( test_mode){
        print("Test mode active")
    }else{
        stop(
            paste0(
                c(
                    "Cannot delete model ",
                    model_name,
                    ", model not found."
                ),
                collapse = ""
            )
        )    
    }
    print(paste0("Deleted model: ",model_name))
}

#' remove_model_music
#'
#' \code{remove_model_music} removes a model from the music library
#'
#' @param model_name Name of the model
#' @param test_mode Testrun indicator
#' @usage
#' remove_model_music(
#'     model_name,
#'     test_mode
#' )
#' @examples
#' remove_model_music(
#'     model_name = "My_model",
#'     test_mode = TRUE
#' )
#' @import stringr
#' @export
remove_model_music = function(
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
        c(system.file("Models/music", package = "artdeco"),
          "/", paste0(model_name, ".RDS")
        ), sep ="", collapse = ""
    )
    
    # check if model exists
    if ( 
        file.exists(model_path)
    ){
        file.remove(model_path)
    }else if ( test_mode){
        print("Test mode active")
    }else{
        stop(
            paste0(
                c(
                    "Cannot delete model ",
                    model_name,
                    ", model not found."
                ),
                collapse = ""
            )
        )    
    }
    print(paste0("Deleted model: ",model_name))
}

#' show_models
#'
#' \code{show_models} shows which models have been trained for what library
#'
#'@param lib_name Specify the library (e.g. "NMF", "music", "bseqsc" or "all")
#'@usage
#' show_models(
#' lib_name
#' )
#'@examples
#' show_models(
#' lib_name = "NMF"
#' )
#'@import stringr
#'@export
#'@return list of models
show_models = function(
    lib_name = "NMF"
){
    
    package_path = system.file("", package = "artdeco")
    
    if (lib_name == "all"){
        model_path_nmf <- paste0(package_path,"Models/NMF")
        model_path_music <- paste0(package_path,"Models/music")
        model_path_bseqsc <- paste0(package_path,"Models/bseqsc")
        
        nmf_models <- stringr::str_replace_all(list.files(model_path_nmf), pattern = ".RDS", "")
        music_models <- stringr::str_replace_all(list.files(model_path_music), pattern = ".RDS", "")
        bseqsc_models <- stringr::str_replace_all(list.files(model_path_bseqsc), pattern = ".RDS", "")
        
        models <- list(nmf_models, music_models, bseqsc_models)
        names(models) <- c("NMF", "music", "bseqsc")
        
    } else {
        model_path = paste0(package_path, "/Models/", lib_name)
        models = stringr::str_replace_all(
            list.files(model_path),
            pattern = ".RDS",
            "")
    }
    
    return(models)
}

#' show_models_music
#'
#' \code{show_models_music} shows which models have been trained
#'
#'@usage
#' show_models_music()
#'@examples
#' show_models_music()
#'@import stringr
#'@export
#'@return list of models
show_models_music = function(
){
    
    package_path = system.file("", package = "artdeco")
    model_path = paste(package_path,"Models/music", sep = "/")
    
    models = stringr::str_replace_all(
        list.files(model_path),
        pattern = ".RDS",
        ""
    )
    
    return(models)
}

#' show_models_bseqsc
#'
#' \code{show_models_bseqsc} shows which models have been trained for bseqsc
#'
#'@usage
#' show_models_bseqsc()
#'@examples
#' show_models_bseqsc()
#'@import stringr
#'@export
#'@return list of models
show_models_bseqsc = function(
){

    package_path = system.file("", package = "artdeco")
    model_path = paste(package_path,"Models/bseqsc", sep = "/")

    models = stringr::str_replace_all(
        list.files(model_path),
        pattern = ".RDS",
        ""
    )

    return(models)
}

#' show_models_NMF
#'
#' \code{show_models_NMF} shows which models have been trained for the NMF method
#'
#'@usage
#' show_models_NMF()
#'@examples
#' show_models_NMF()
#'@import stringr
#'@export
#'@return list of models
show_models_NMF = function(
){
    
    package_path = system.file("", package = "artdeco")
    model_path = paste(package_path,"Models/NMF", sep = "/")
    
    models = stringr::str_replace_all(
        list.files(model_path),
        pattern = ".RDS",
        ""
    )
    
    return(models)
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
#'     parameter_list
#' )
#'@import stringr
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
            names(parameter) = str_to_lower(names(parameter))
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
                (subtype_sim_scalar / maximum) * 100,
                3
            )

            # load data

            meta_data[,subtype_label_absolute] = rep("low",nrow(meta_data))
            meta_data[,subtype_label_quant_absolute] = rep(0,nrow(meta_data))
            meta_data[,subtype_label_quant_absolute] = similarity_quantified_absolute
            meta_data[ 
                similarity_quantified_absolute >= mean(similarity_quantified_absolute),
                subtype_label_absolute
            ] = "high"
            
            meta_data[
                ( meta_data[ , subtype_label_absolute ] == "") |
                ( is.na(meta_data[ , subtype_label_absolute ])
                  ),
                subtype_label_absolute
            ] = "none"

        #case baseline relative

            similarity_quantified_relative = round(
                subtype_sim_scalar / max( subtype_sim_scalar ) * 100,
                0
            )
            #similarity_quantified_relative = log(similarity_quantified_relative+1)
            #similarity_quantified_relative = log(similarity_quantified_relative+1)

            # load data
            
            meta_data[,subtype_label_quant_relative] = round(as.double(
                similarity_quantified_relative
            ),0)
            meta_data[,subtype_label_relative] = rep("low",nrow(meta_data))
            meta_data[
                similarity_quantified_relative>mean(similarity_quantified_relative),
                subtype_label_relative
            ] = "high"
            
            meta_data[
                (meta_data[ ,
                            subtype_label_relative ] == "") | ( is.na(meta_data[ , subtype_label_relative
                                                                                 ]) ),
                subtype_label_relative
            ] = "none"

    }
    return(meta_data)
}

context("Sanity checks")

test_that("Assert that data can be analyzed by Deconvolve_transcriptome", {

    library("artdeco")
    
    data("visualization_data")
    data("deconvolution_results")

    expect_equal(dim(visualization_data), c(12997, 97))

    deconvolution_results_test = Deconvolve_transcriptome(
         transcriptome_data = visualization_data
    )

    expect_equal(dim(deconvolution_results), c(97, 14))

    expect_identical(
            colnames(deconvolution_results_test),
            c("model", "alpha", "beta", "gamma", "delta", "acinar", "ductal", "hisc","Strength_subtype", "Subtype", "score")
        )
    
    expect_error(deconvolution_results_test = Deconvolve_transcriptome( # algorithm does not exist
        transcriptome_data = visualization_data,
        deconvolution_algorithm = "DeconRNASeq"
    ))
    expect_error(deconvolution_results_test = Deconvolve_transcriptome( # model does not exist
        transcriptome_data = visualization_data,
        models = "a_model"
    ))
    
})



test_that("Test if visualization works", {

    library("artdeco")
    library("stringr")
    
    deconvolution_results = deconvolution_results[
        colnames(deconvolution_results) != "Strength_subtype"
    ]
    data(visualization_data, envir = environment())

    decon_heatmap <- create_heatmap_deconvolution(
        visualization_data = visualization_data,
        deconvolution_results = deconvolution_results
    )
    
    expect_identical(decon_heatmap$tree_col$labels, colnames(visualization_data))
    
    decon_pca <- create_PCA_deconvolution(
        visualization_data = visualization_data,
        deconvolution_results = deconvolution_results
    )
    
    expect_identical(decon_pca$data$labels, colnames(visualization_data))
    expect_true(str_detect(decon_pca$labels[[2]], "standardized PC1"))
    expect_true(str_detect(decon_pca$labels[[1]], "standardized PC2"))
    
})



test_that("Adding, showing and removing models", {

    data(meta_data)
    
    expect_equal(dim(meta_data), c(638, 15))
    
    rownames(meta_data) = as.character(meta_data$Name)

    subtype_vector = as.character(meta_data$Subtype) # extract the training sample subtype labels

    data("Lawlor")
    expect_equal(dim(Lawlor), c(20655, 638))

    # add new model
    model_name = "Test_model"
    model_path = paste(
        c(system.file("Models/NMF", package="artdeco"),"/", model_name,".RDS"),
        collapse = ""
    )
    
    expect_error(add_deconvolution_training_model_NMF( # no model name given
        transcriptome_data = Lawlor,
        subtype_vector = subtype_vector
    ))
    
    expect_error(add_deconvolution_training_model_NMF( # model already exists
        transcriptome_data = Lawlor,
        model_name = "Alpha_Beta_Gamma_Delta_Acinar_Ductal_Baron",
        subtype_vector = subtype_vector
    ))
    
    expect_error(add_deconvolution_training_model_NMF( # subtype_vector is not character vector
        transcriptome_data = Lawlor,
        model_name = model_name,
        subtype_vector = meta_data$Subtype
    ))
    
    expect_error(add_deconvolution_training_model_NMF( # no subtype labels given
        transcriptome_data = Lawlor,
        model_name = model_name,
        subtype_vector = c()
    ))
    
    add_deconvolution_training_model_NMF(
        transcriptome_data = Lawlor,
        model_name = model_name,
        subtype_vector = subtype_vector
    )
    
    # test existence of newly added model
    expect_true(file.exists(model_path))
    nmf_models <- show_models_NMF()
    nmf_models_too <- show_models(lib_name = "NMF")
    
    expect_identical(tail(nmf_models, 1), model_name)
    expect_identical(tail(nmf_models_too, 1), model_name)

    # test if model is non-empty
    new_model <- readRDS(model_path)
    
    expect_equal(dim(new_model@fit@W), c(595, 8))
    expect_equal(dim(new_model@fit@H), c(8, 638))
    
    expect_identical(colnames(new_model@fit@H), colnames(Lawlor))
    

    # remove model
   remove_model_NMF(model_name = model_name) 
   # remove_model(model_name = model_name, lib_name = "NMF")
    
   expect_true(!file.exists(model_path))

})

context("Sanity checks")

test_that("Assert that data can be analyzed by Deconvolve_transcriptome", {

    library("artdeco")
    
    data("visualization_data")
    data("deconvolution_results")

    expect_equal(nrow(visualization_data), 12997)
    expect_equal(ncol(visualization_data), 97)

    deconvolution_results_test = Deconvolve_transcriptome(
         transcriptome_data = visualization_data
    )

    expect_equal(nrow(deconvolution_results), 97)
    expect_equal(ncol(deconvolution_results), 14)

    expect_identical(
            colnames(deconvolution_results_test),
            c("model", "alpha", "beta", "gamma", "delta", "acinar", "ductal", "hisc","Strength_subtype", "Subtype", "score")
        )
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
    
    expect_equal(nrow(meta_data), 638)
    expect_equal(ncol(meta_data), 15)
    
    rownames(meta_data) = as.character(meta_data$Name)

    subtype_vector = as.character(meta_data$Subtype) # extract the training sample subtype labels

    data("Lawlor")
    expect_equal(nrow(Lawlor), 20655)
    expect_equal(ncol(Lawlor), 638)

    # add new model
    model_name = "My_model"
    model_path = paste(
        c(system.file("Models/NMF", package="artdeco"),"/", model_name,".RDS"),
        collapse = ""
    )
    
    add_deconvolution_training_model_NMF(
        transcriptome_data = Lawlor,
        model_name = model_name,
        subtype_vector = subtype_vector,
        rank_estimate = 0,
        exclude_non_interpretable_NMF_components = FALSE,
        training_nr_marker_genes = 100,
        parallel_processes = 1,
        nrun = 1
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
    if (file.exists(model_path))
        remove_model_NMF(
            model_name = model_name
        ) # or remove_model(model_name = model_name, lib_name = "NMF")
    
    expect_true(!file.exists(model_path))

})

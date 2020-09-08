context("Sanity checks")

test_that("Assert that data can be analyzed by Deconvolve_transcriptome", {

    library("artdeco")
    data("visualization_data")
    data("deconvolution_results")

    expect_true( nrow(visualization_data) == 12997 ) # do same for deconvolution_results
    expect_true( ncol(visualization_data) == 97 )

    deconvolution_results_test = Deconvolve_transcriptome(
         transcriptome_data = visualization_data
    )

    #expect_true(
    #    deconvolution_results_test == deconvolution_results)
    expect_true(
        nrow(deconvolution_results) == 194)

    expect_true(
        ncol(deconvolution_results) == 11)

    expect_true(
        colnames(deconvolution_results_test)==c("model", "alpha", "beta", "gamma", "delta", "acinar", "ductal", "hisc", 
                                                "Strength_subtype", "Subtype", "score"))
})



test_that("Test if visualization works", {

    library("artdeco")
    deconvolution_results = deconvolution_results[
        colnames(deconvolution_results) != "Strength_subtype"
    ]
    data(visualization_data, envir = environment())

    create_heatmap_deconvolution(
        visualization_data,
        deconvolution_results = deconvolution_results
    )
})



test_that("Adding models", {
    ### adding models

    data(meta_data)
    
    expect_true( nrow(meta_data) == 638 )
    expect_true( ncol(meta_data) == 15 )
    
    rownames(meta_data) = as.character(meta_data$Name)

    subtype_vector = meta_data$Subtype # extract the training sample subtype labels

    data("Lawlor")
    expect_true( nrow(Lawlor) == 20655)
    expect_true( ncol(Lawlor) == 638)

    model_name = "my_model"
    model_path = paste(
        c(system.file("Models/NMF", package="artdeco"),"/",model_name,".RDS"),
        collapse = ""
    )

    if (file.exists(model_path))
        file.remove(model_path)

    data(meta_data)
     
    subtype_vector = as.character(meta_data$Subtype) # extract the training sample subtype labels
     
    add_deconvolution_training_model_NMF(
         transcriptome_data = Lawlor,
         model_name = "my_model",
         subtype_vector = subtype_vector,
         rank_estimate = 0,
         exclude_non_interpretable_NMF_components = FALSE,
         training_nr_marker_genes = 100,
         parallel_processes = 1,
         nrun = 1
    )
    
    expect_true(file.exists(model_path))

    # remove model
    remove_model_NMF(
        model_name = "My_model",
        test_mode = TRUE
    )
})

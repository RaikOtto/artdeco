context("Sanity checks")

test_that("Assert that data can be analyzed by Determine_differentiation_stage", {

    library("artdeco")
    path_transcriptome_file = system.file(
        "/Data/Expression_data/PANnen_Test_Data.tsv",
        package="artdeco"
    )

    expect_true( length(path_transcriptome_file) > 0)

    # testrun
    deconvolution_results = Determine_differentiation_stage(
        transcriptome_file_path = path_transcriptome_file
    )

    expect_true(
        nrow(deconvolution_results) == 57)

    expect_true(
        ncol(deconvolution_results) == 12)

    deconvolution_results_path = system.file(
        "/Data/Expression_data/",
        package="artdeco"
    )
    deconvolution_results_path = paste0(
        deconvolution_results_path,
        "/Deconvolution_test_data.tsv"
    )
})

test_that("Test if visualization works", {

    library("artdeco")
    visualization_data_path = system.file(
        "/Data/Expression_data/Visualization_PANnen.tsv",
        package="artdeco")
    expect_true( length(visualization_data_path) > 0)

    meta_data_path = system.file(
        "Data/Meta_information/Meta_information.tsv",
        package = "artdeco"
    )
    meta_data      = read.table(
        meta_data_path, sep ="\t",
        header = TRUE,
        stringsAsFactors = FALSE
    )
    rownames(meta_data) = meta_data$Sample

    create_heatmap_differentiation_stages(
        visualization_data_path,
        deconvolution_results = meta_data
    )
})

test_that("Adding models", {
    ### adding models

    meta_data_path = system.file(
        "Data/Meta_information/Meta_information.tsv",
        package = "artdeco"
    )
    meta_data      = read.table(
        meta_data_path, sep ="\t", header = TRUE,
        stringsAsFactors = FALSE)
    expect_true( length(meta_data_path) > 0)

    subtype_vector = meta_data$Subtype # extract the training sample subtype labels
    expect_true( length(subtype_vector) > 0)

    training_data_path = system.file(
        "Data/Expression_data/PANnen_Test_Data.tsv", package = "artdeco")
    expect_true( length(training_data_path) > 0)

    model_name = "Test_model"
    model_path = paste(
        c(system.file("Models/", package="artdeco"),"/",model_name,".RDS"),
        collapse = ""
    )

    if (file.exists(model_path))
        file.remove(model_path)

    add_deconvolution_training_model(
        transcriptome_data_path = training_data_path,
        model_name = "Test_model",
        subtype_vector,
        training_p_value_threshold = 0.05,
        training_nr_permutations = 100,
        training_nr_marker_genes = 100
    )
    expect_true(file.exists(model_path))

    # remove model
    remove_model(model_name)
    expect_false(file.exists(model_path))
})

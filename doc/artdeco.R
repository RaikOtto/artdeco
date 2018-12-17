## ------------------------------------------------------------------------
library("artdeco")

# obtain path to testdata
path_transcriptome_file = system.file(
    "/Data/Expression_data/PANnen_Test_Data.tsv",
    package="artdeco")

# testrun
deconvolution_results = Determine_differentiation_stage(
    transcriptome_file_path = path_transcriptome_file
)

# show results
deconvolution_results[1:5,1:5]

## ------------------------------------------------------------------------
visualization_data_path = system.file(
    "/Data/Expression_data/Visualization_PANnen.tsv",
    package="artdeco")

create_heatmap_differentiation_stages(
    visualization_data_path,
    deconvolution_results
)

## ------------------------------------------------------------------------
meta_data_path = system.file("Data/Meta_Data.tsv", package = "artdeco")
meta_data      = read.table(
    meta_data_path, sep ="\t", header = TRUE,
    stringsAsFactors = FALSE)
subtype_vector = meta_data$Subtype # extract the training sample subtype labels

subtype_vector[1:6] # show subtype definition

training_data_path = system.file(
    "Data/Expression_data/PANnen_Test_Data.tsv", package = "artdeco")

add_deconvolution_training_model(
    training_data = training_data_path,
    model_name = "My_model",
    subtype_vector = subtype_vector,
    training_nr_marker_genes = 5
)

## ------------------------------------------------------------------------
models = list.files(system.file("Models/", package = "artdeco"))
print(models)
model_to_be_removed = models[length(models)]
remove_model(
    model_name = model_to_be_removed
)

## ------------------------------------------------------------------------
training_data_path = system.file(
    "Data/Expression_data/PANnen_Test_Data.tsv", package = "artdeco")
expression_matrix = read.table(
    training_data_path,
    sep ="\t",
    header = TRUE,
    row.names = 1)

expression_matrix[1:5,1:5]


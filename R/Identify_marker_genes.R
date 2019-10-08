#' identify_marker_genes
#' @param expression_training_mat matrix that contains the expression data
#' @param subtype_vector Vector containing the subtype of the training samples
#' @param subtype Subtype for which to determine the differentially expressed
#' @param training_nr_marker_genes Number of desired marker genes to describe a subtype
#' @import stringr limma
#' @return Returns marker genes
identify_marker_genes = function(
    expression_training_mat,
    subtype_vector,
    subtype,
    training_nr_marker_genes
){

    print(paste("Calculating marker genes for subtype: ",subtype, sep =""))

    groups = subtype_vector
    groups[groups == subtype] = "CASE"
    groups[groups != "CASE"] = "CTRL"
    design <- model.matrix(~0 + groups)
    colnames(design) = c("Case","Ctrl")

    vfit = lmFit(expression_training_mat,design)
    contr.matrix = makeContrasts(
        contrast = Case - Ctrl,
        levels = design
    )

    vfit = contrasts.fit( vfit, contrasts = contr.matrix)
    efit = eBayes(vfit)

    result_t = topTable(
        efit,
        coef     = "contrast",
        number   = nrow(expression_training_mat),
        genelist = efit$genes,
        "none",
        sort.by  = "B",
        resort.by= NULL,
        p.value  = 1,
        lfc      = 0,
        confint  = FALSE
    )
    result_t$hgnc_symbol = rownames(result_t)
    colnames(result_t) = c("Log_FC","Average_Expr","t","P_value","adj_P_value","B","HGNC")

    result_t = result_t[c("HGNC","Log_FC","Average_Expr","P_value","adj_P_value")]
    result_t = result_t[order(result_t$P_value, decreasing = FALSE),]
    result_t$Log_FC = round(result_t$Log_FC, 1)
    result_t$Average_Expr = round(result_t$Average_Expr, 1)
    result_t = result_t[order(result_t$Log_FC,decreasing = TRUE),]

    result_t = result_t[order(result_t$Log_FC, decreasing = TRUE),]
    marker_genes = as.character(result_t$HGNC)
    marker_genes = marker_genes[marker_genes %in% rownames(expression_training_mat)]
    marker_genes = marker_genes[1:training_nr_marker_genes]

    return(result_t)
}

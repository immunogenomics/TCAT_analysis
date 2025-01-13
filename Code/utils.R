library(SingleCellExperiment)
library(Matrix)

export_sce <- function(sce, filepath_prefix) {
    # Export expression matrix to Matrix Market format
    expr_matrix <- assay(sce, "counts") # Adjust as needed for the desired assay
    if(!inherits(expr_matrix, "sparseMatrix")) {
        expr_matrix <- Matrix(expr_matrix, sparse = TRUE)
    }
    writeMM(expr_matrix, paste0(filepath_prefix, "_expression.mtx"))
    
    # Export cell metadata to CSV
    cell_meta <- as.data.frame(colData(sce))
    write.csv(cell_meta, file=paste0(filepath_prefix, "_cell_metadata.csv"), quote=FALSE, row.names=TRUE)
    
    # Export gene metadata to CSV
    gene_meta <- as.data.frame(rowData(sce))
    write.csv(gene_meta, file=paste0(filepath_prefix, "_gene_metadata.csv"), quote=FALSE, row.names=TRUE)
}


export_dgcmatrix <- function(data, filepath_prefix) {
    # Export expression matrix to Matrix Market format
    writeMM(obj = data, paste0(filepath_prefix, "_counts.mtx"))
    
    # Export cell metadata to CSV
    cell_meta <- as.data.frame(colnames(data))
    write.csv(cell_meta, file=paste0(filepath_prefix, "_cell_metadata.csv"), quote=FALSE, row.names=TRUE)
    
    # Export gene metadata to CSV
    gene_meta <- as.data.frame(rownames(data))
    write.csv(gene_meta, file=paste0(filepath_prefix, "_gene_metadata.csv"), quote=FALSE, row.names=TRUE)
}
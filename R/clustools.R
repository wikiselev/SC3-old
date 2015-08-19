#' Default parameters for the first filtering step
#'
#' Returns default min.cells, max.cells, min.reads parameters used in
#' gene_filter1()
#'
#' @param dataset Name of the toy dataset (either "quake", "sandberg",
#' "linnarsson" or "bernstein")
#' @return min.cells, max.cells, min.reads required to be able to filter genes
#' using filter1.
#' @examples
#' filter1_params("quake")
filter1_params <- function(dataset) {
    n.cells <- dim(dataset)[2]

    min.cells <- ceiling(0.06*n.cells)
    max.cells <- ceiling(0.06*n.cells)
    min.reads <- 2

    return(list(min.cells = min.cells, max.cells = max.cells, min.reads = min.reads))
}

#' First filtering step
#'
#' Filter genes that would not contribute to clustering, because they are either
#' expreseed or not expressed in almost all cells.
#'
#' @param d Expression matrix with rows as genes and columns as cells
#' @param min.cells Minimum number of cells in which a given gene is expressed
#' @param max.cells Maximum number of cells in which a given gene is expressed
#' @param min.reads Minimum number of reads per gene per cell
#' @return Filtered expression matrix in which only genes that are expressed in
#' more than \code{min.cells} with more than \code{min.reads} reads and also are
#' expressed in less than [total number of cells - \code{max.cells}].
#' @examples
#' gene_filter1(quake, 3, 3, 2)
gene_filter1 <- function(d, min.cells, max.cells, min.reads) {
    d <- d[rowSums(d > min.reads) >= min.cells & rowSums(d > 0) <= dim(d)[2] - max.cells, ]
    d <- unique(d)
    return(d)
}

support_vector_machines1 <- function(teach, study, kern) {
    cat("Dimensions of teacher:\n")
    cat(dim(teach))
    cat("\n")
    cat("Dimensions of study:\n")
    cat(dim(study))
    cat("\n")

    teach <- t(teach)
    labs <- factor(rownames(teach))
    rownames(teach) <- NULL
    # length(unique(colnames(teach)))
    cat("Performing svm...\n")
    model <- tryCatch(svm(teach, labs, kernel = kern),
                      error = function(cond) return(NA))
    cat("Performing prediction...\n")
    # if(!is.na(model)) {
    pred <- predict(model, t(study))
    return(pred = pred)
    # } else {
    # return(NA)
    # }
}

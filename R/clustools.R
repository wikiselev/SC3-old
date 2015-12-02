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
filter1_params <- function(dataset, fraction) {
    n.cells <- dim(dataset)[2]

    min.cells <- ceiling(fraction*n.cells)
    max.cells <- ceiling(fraction*n.cells)
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

#' Calculate a distance matrix
#'
#' Calculate a distance between column vectors of the input dataset using a specified
#' distance metrics.
#'
#' @param d Expression matrix with rows as genes and columns as cells
#' @param method Distance metrics: "spearman", "pearson", "euclidean", "maximum",
#' "manhattan", "canberra", "binary" or "minkowski"
#' @return A distance matrix
#' @examples
#' calculate_distance(quake, "spearman")
calculate_distance <- function(d, method) {
    return(if (method == "spearman") {
        # there is no spearman distance method in 'proxy' package - have to define manually
        as.matrix(1 - cor(d, method = "spearman"))
        # the output is bit different from Martin's results - need to figure out
    } else if (method == "pearson") {
        as.matrix(1 - cor(d, method = "pearson"))
    } else {
        as.matrix(dist(t(d), method = method))
    })
}

#' Dimensionality reduction of a distance matrix
#'
#' Transform a distance matrix to a new basis
#'
#' @param dists A distance matrix
#' @param method Dimensionality reduction method: "pca", "spectral", "spectral_reg"
#' or "mds"
#' @return A transformed distance matrix
#' @examples
#' transformation(d, "spectral")
transformation <- function(dists, method) {
    if (method == "pca") {
        t <- prcomp(dists, center = TRUE, scale. = TRUE)
        list(t$rotation, t$sdev)
    } else if (method == "spectral") {
        L <- norm_laplacian(exp(-dists/max(dists)), 0)
        # here need to sort eigenvectors by their eigenvalues in increasing order!
        list(eigen(L)$vectors[, order(eigen(L)$values)], eigen(L)$values[order(eigen(L)$values)])
    } else if (method == "spectral_reg") {
        L <- norm_laplacian(exp(-dists/max(dists)), 1000)
        # here need to sort eigenvectors by their eigenvalues in increasing order!
        list(eigen(L)$vectors[, order(eigen(L)$values)], eigen(L)$values[order(eigen(L)$values)])
    } else if (method == "mds") {
        t <- cmdscale(dists, k = ncol(dists) - 1)
        list(t)
    }
}

support_vector_machines1 <- function(teach, study, kern) {
#     cat("Dimensions of teacher:\n")
#     cat(dim(teach))
#     cat("\n")
#     cat("Dimensions of study:\n")
#     cat(dim(study))
#     cat("\n")

    teach <- t(teach)
    labs <- factor(rownames(teach))
    rownames(teach) <- NULL
    # length(unique(colnames(teach)))
    # cat("Performing svm...\n")
    model <- tryCatch(svm(teach, labs, kernel = kern),
                      error = function(cond) return(NA))
    # cat("Performing prediction...\n")
    # if(!is.na(model)) {
    pred <- predict(model, t(study))
    return(pred = pred)
    # } else {
    # return(NA)
    # }
}

#' Laplacian computation
#'
#' @param x Adjacency/distance matrix
#' @param tau Regularization term
#' @return Laplacian of the adjacency/distance matrix
#' @examples
#' L <- norm_laplacian(dists, 0)
norm_laplacian <- function(x, tau) {
    x <- x + tau * matrix(1, dim(x)[1], dim(x)[2])
    D <- diag(colSums(x))
    D1 <- D^(-0.5)
    D1[D1 == Inf] <- 0
    return(diag(dim(D)[1]) - D1 %*% x %*% D1)
}

# consensus clustering analysis Cluster-based similarity partitioning algorithm
consensus_clustering <- function(clusts) {
    n.cells <- length(unlist(strsplit(clusts[1], " ")))
    res <- matrix(0, nrow = n.cells, ncol = n.cells)
    for (i in 1:length(clusts)) {
        t <- clusts[i]
        t <- as.numeric(unlist(strsplit(t, " ")))
        t <- as.matrix(dist(t))
        t[t != 0] <- -1
        t[t == 0] <- 1
        t[t == -1] <- 0
        res <- res + t
    }
    res <- res/i
    return(res)
}

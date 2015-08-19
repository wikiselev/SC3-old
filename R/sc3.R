run_sc3 <- function(filename, ks) {
    if(filename %in% c("quake", "quake_all_fpkm", "quake_all_read", "sandberg", "sandberg_all_read", "sandberg_all_rpkm",
                       "bernstein", "linnarsson", "zhong", "kirschner")) {
        dataset <- get(filename)
    } else {
        dataset <- read.table(filename, header = T)
        if(filename == "../clustering/tallulah/data/Bergiers_Exp1_FeatureCounts_BaseGenomeAnnotations_No_Multimapping.out" |
           filename == "../clustering/tallulah/data/Bergiers_Exp2_FeatureCounts_BaseGenomeAnnotations_No_Multimapping.out" |
           filename == "../clustering/tallulah/data/Bergiers_Exp1_FeatureCounts_BaseGenomeAnnotations_No_Multimapping_RUVnorm.txt" |
           filename == "../clustering/tallulah/data/Bergiers_Exp2_FeatureCounts_BaseGenomeAnnotations_No_Multimapping_RUVnorm.txt" |
           filename == "../clustering/tallulah/data/Bergiers_Exp1_FeatureCounts_BaseGenomeAnnotations_No_Multimapping_SFnorm.txt" |
           filename == "../clustering/tallulah/data/Bergiers_Exp2_FeatureCounts_BaseGenomeAnnotations_No_Multimapping_SFnorm.txt") {
            rownames(dataset) <- dataset[, 1]
            dataset <- dataset[ , 2:dim(dataset)[2]]
        }
    }

    # hard cell filter
    # more than 2000 genes have to be expressed in each cell
    dataset <- dataset[ , colSums(dataset > 1e-2) > 2000]

    if(dim(dataset)[2] == 0) {
        cat("Your dataset did not pass cell filter (more than 2000 genes have to be expressed in each cell)! Stopping now...")
        return()
    }

    svm.num.cells <- 1000
    distances <- c("euclidean", "pearson", "spearman")
    dimensionality.reductions <- c("pca", "spectral")

    cat("1. Preliminary gene filtering...\n")
    filter1.params <- filter1_params(dataset)
    min.cells <- filter1.params$min.cells
    max.cells <- filter1.params$max.cells
    min.reads <- filter1.params$min.reads
    dataset <- gene_filter1(dataset, min.cells, max.cells, min.reads)

    if(dim(dataset)[1] == 0) {
        cat("Your dataset did not pass gene filter! Stopping now...")
        return()
    }

    cat("2. Log2-transforming data...\n")
    if(filename != "bernstein") {
        dataset <- log2(1 + dataset)
    }

    study.dataset <- data.frame()
    if(dim(dataset)[2] > svm.num.cells) {
        cat("\n")
        cat(paste0("Your dataset contains more than ", svm.num.cells, " cells, therefore clustering wil be performed on a random sample of ", svm.num.cells, " cells, the rest of the cells will be predicted using SVM."))
        cat("\n")
        cat("\n")
        working.sample <- sample(1:dim(dataset)[2], svm.num.cells)
        study.dataset <- dataset[ , setdiff(1:dim(dataset)[2], working.sample)]
        dataset <- dataset[, working.sample]
    }

    n.cells <- dim(dataset)[2]
    n.dim <- floor(0.05 * n.cells) : ceiling(0.08 * n.cells)

    if(length(n.dim) > 15) {
        n.dim <- sample(n.dim, 15)
    }

    hash.table <- expand.grid(distan = distances,
                              dim.red = dimensionality.reductions,
                              k = c(min(ks) - 1, ks),
                              n.dim = n.dim, stringsAsFactors = F)

    # register local cluster
    cl <- makeCluster(detectCores() - 1, outfile="")
    registerDoParallel(cl, cores = detectCores() - 1)

    cat("3. Calculating distance matrices...\n")
    dists = foreach(i = distances) %dopar% {
        try({
                if (i == "spearman") {
                    # there is no spearman distance method in 'proxy' package - have to define manually
                    as.matrix(1 - cor(dataset, method = "spearman"))
                    # the output is bit different from Martin's results - need to figure out
                } else if (i == "pearson") {
                    as.matrix(1 - cor(dataset, method = "pearson"))
                } else {
                    as.matrix(dist(t(dataset), method = i))
                }
        })
    }
    names(dists) <- distances

    pb <- txtProgressBar(min = 1, max = dim(hash.table)[1], style = 3)

    cat("4. Performing dimensionality reduction and kmeans clusterings...\n")
    labs = foreach(i = 1:dim(hash.table)[1],
                   .combine = rbind) %dopar% {
        try({
            norm_laplacian <- function(x, tau) {
                x <- x + tau * matrix(1, dim(x)[1], dim(x)[2])
                D <- diag(colSums(x))
                D1 <- D^(-0.5)
                D1[D1 == Inf] <- 0
                return(diag(dim(D)[1]) - D1 %*% x %*% D1)
            }

            dist <- get(hash.table[i, 1], dists)
            method <- hash.table[i, 2]
            if (method == "pca") {
                t <- prcomp(dist, center = TRUE, scale. = TRUE)
                t <- t$rotation
            } else if (method == "spectral") {
                L <- norm_laplacian(exp(-dist/max(dist)), 0)
                # here need to sort eigenvectors by their eigenvalues in increasing order!
                t <- eigen(L)$vectors[, order(eigen(L)$values)]
            } else if (method == "spectral_reg") {
                L <- norm_laplacian(exp(-dist/max(dist)), 1000)
                # here need to sort eigenvectors by their eigenvalues in increasing order!
                t <- eigen(L)$vectors[, order(eigen(L)$values)]
            } else if (method == "mds") {
                t <- cmdscale(dist, k = ncol(dist) - 1)
            }

            s <- paste(kmeans(t[, 1:hash.table[i, 4]],
                         hash.table[i, 3],
                         iter.max = 1e+09,
                         nstart = 1000)$cluster,
                  collapse = " ")
            setTxtProgressBar(pb, i)
            return(s)
        })
    }

    close(pb)

    res <- cbind(hash.table, labs)
    res$labs <- as.character(res$labs)
    rownames(res) <- NULL

    cat("5. Computing consensus matrix and labels...\n")
    all.combinations <- NULL
    for(k in c(min(ks) - 1, ks)) {
        for(i in 1:length(distances)) {
            for(j in 1:length(dimensionality.reductions)) {
                dist.combs <- combn(distances, i)
                dim.red.combs <- combn(dimensionality.reductions, j)
                for(m in 1:dim(dist.combs)[2]) {
                    for(n in 1:dim(dim.red.combs)[2]) {
                        all.combinations <- rbind(
                            all.combinations,
                            cbind(paste(dist.combs[, m], collapse = " "),
                                  paste(dim.red.combs[, n], collapse = " "),
                                  as.numeric(k)))
                    }
                }
            }
        }
    }

    cons = foreach(i = 1:dim(all.combinations)[1]) %dopar% {
        try({

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

            d <- res[res$distan %in% strsplit(all.combinations[i, 1], " ")[[1]] &
                         res$dim.red %in% strsplit(all.combinations[i, 2], " ")[[1]] &
                         res$k == as.numeric(all.combinations[i, 3]), ]

            dat <- consensus_clustering(d$labs)

            diss <- dist(dat)
            clust <- hclust(diss)
            clusts <- cutree(clust, k = as.numeric(all.combinations[i, 3]))

            labs <- NULL
            for(j in unique(clusts[clust$order])) {
                labs <- rbind(labs, paste(names(clusts[clusts == j]), collapse = " "))
            }

            labs <- as.data.frame(labs)
            colnames(labs) <- "Labels"

            return(list(dat, labs, clust))
        })
    }



    # stop local cluster
    stopCluster(cl)

    run_shiny_app(filename, distances, dimensionality.reductions,
                   cbind(all.combinations, cons),
                   dataset, study.dataset, svm.num.cells)
}

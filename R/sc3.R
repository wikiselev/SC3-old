run_sc3 <- function(filename, ks) {
    if(filename == "quake_all_fpkm") {
        dataset <- get(filename)
    } else {
        if(!grepl("csv", filename)) {
            dataset <- read.table(filename)
        } else if(grepl("csv", filename)) {
            dataset <- read.csv(filename, header = F)
        }
        rownames(dataset) <- dataset[, 1]
        dataset <- dataset[ , 2:dim(dataset)[2]]
    }

    # original.dataset <- dataset

#     # hard cell filter
#     # more than 2000 genes have to be expressed in each cell
#     dataset <- dataset[ , colSums(dataset > 1e-2) > 2000]
#
#     if(dim(dataset)[2] == 0) {
#         cat("Your dataset did not pass cell filter (more than 2000 genes have to be expressed in each cell)! Stopping now...")
#         return()
#     }

    svm.num.cells <- 50
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
        study.sample <- setdiff(1:dim(dataset)[2], working.sample)
        study.dataset <- dataset[ , study.sample]
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
    dists = foreach(i = distances, .packages = "SC3") %dopar% {
        try({
            calculate_distance(dataset, i)
        })
    }
    names(dists) <- distances

    pb <- txtProgressBar(min = 1, max = dim(hash.table)[1], style = 3)

    cat("4. Performing dimensionality reduction and kmeans clusterings...\n")
    labs = foreach(i = 1:dim(hash.table)[1], .packages = "SC3",
                   .combine = rbind) %dopar% {
        try({
            t <- transformation(get(hash.table[i, 1], dists), hash.table[i, 2])[[1]]
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

    cons = foreach(i = 1:dim(all.combinations)[1], .packages = "SC3") %dopar% {
        try({
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
                   dataset, study.dataset, svm.num.cells, working.sample, study.sample)
}

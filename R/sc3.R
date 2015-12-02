
get_data <- function(name) {
    if(!is.character(name)) {
        return(name)
    } else {
        if(!grepl("csv", name)) {
            return(as.matrix(read.table(name, check.names = F)))
        } else if(grepl("csv", name)) {
            return(as.matrix(read.csv(name, check.names = F)))
        }
    }
}

cell_filter <- function(data, cell.filter.genes) {
    # more than cell.filter.genes have to be expressed in each cell
    # by default cell.filter.genes = 2000
    # this criterium is taken from bernstein paper:
    # Patel, A. P. et al. Single-cell RNA-seq highlights intratumoral heterogeneity
    # in primary glioblastoma. Science 344, 1396–1401 (2014).
    cat("Cell filtering...\n")
    data <- data[ , colSums(data > 1e-2) > cell.filter.genes]
    if(dim(data)[2] == 0) {
        cat(paste0("Your dataset did not pass cell filter (more than ", cell.filter.genes, " genes have to be expressed in each cell)! Stopping now..."))
        return()
    } else {
        return(data)
    }
}

gene_filter <- function(data, fraction) {
    cat("Gene filtering and log2-scaling...\n")
    filter1.params <- filter1_params(data, fraction)
    min.cells <- filter1.params$min.cells
    max.cells <- filter1.params$max.cells
    min.reads <- filter1.params$min.reads
    data <- gene_filter1(data, min.cells, max.cells, min.reads)
    data <- log2(1 + data)

    if(dim(data)[1] == 0) {
        cat("All genes were removed after the gene filter! Stopping now...")
        return()
    } else {
        return(data)
    }
}

sc3 <- function(filename, ks = 3:7, cell.filter = F, interactivity = T, svm.num.cells = 1000, cell.filter.genes = 2000, gene.filter.fraction = 0.06, show.original.labels = F, d.region.min = 0.04, d.region.max = 0.07, chisq.quantile = 0.9999) {

    # initial parameters
    set.seed(1)
    distances <- c("euclidean", "pearson", "spearman")
    dimensionality.reductions <- c("pca", "spectral")
    RSelenium::startServer()

    # get input data
    dataset <- get_data(filename)

    # remove duplicated genes
    dataset <- dataset[!duplicated(rownames(dataset)), ]

    # cell filter
    if(cell.filter) {
        dataset <- cell_filter(dataset, cell.filter.genes)
    }

    # gene filter
    if(deparse(substitute(filename)) != "bernstein" & deparse(substitute(filename)) != "kedar.norm") {
        dataset <- gene_filter(dataset, gene.filter.fraction)
    }

    # define the output file basename
    filename <- ifelse(!is.character(filename), deparse(substitute(filename)), basename(filename))

    # define cell names from the input dataset
    cell.names <- c(1:dim(dataset)[2])
    cell.names <- colnames(dataset)

    # prepare for SVM (optional)
    study.cell.names <- NULL
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

        study.cell.names <- study.sample
        study.cell.names <- colnames(study.dataset)

        cell.names <- working.sample
        cell.names <- colnames(dataset)
    }

    # define number of cells and region of dimensions
    n.cells <- dim(dataset)[2]
    n.dim <- floor(d.region.min * n.cells) : ceiling(d.region.max * n.cells)

    # for large datasets restrict the region of dimensions to 15
    if(length(n.dim) > 15) {
        n.dim <- sample(n.dim, 15)
    }

    # create a hash table for running on parallel CPUs
    hash.table <- expand.grid(distan = distances,
                              dim.red = dimensionality.reductions,
                              k = c(min(ks) - 1, ks),
                              n.dim = n.dim, stringsAsFactors = F)

    # register computing cluster (n-1 CPUs) on a local machine
    cl <- makeCluster(detectCores() - 1, outfile="")
    registerDoParallel(cl, cores = detectCores() - 1)

    # calculate distances in parallel
    cat("Calculating distance matrices...\n")
    dists = foreach(i = distances, .packages = "SC3", .options.RNG=1234) %dorng% {
        try({
            calculate_distance(dataset, i)
        })
    }
    names(dists) <- distances

    # perform kmeans in parallel
    # add a progress bar to be able to see the progress
    pb <- txtProgressBar(min = 1, max = dim(hash.table)[1], style = 3)
    cat("Performing dimensionality reduction and kmeans clusterings...\n")
    labs = foreach(i = 1:dim(hash.table)[1], .packages = "SC3",
                   .combine = rbind, .options.RNG=1234) %dorng% {
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

    # perform consensus clustering in parallel
    cat("Computing consensus matrix and labels...\n")
    # first make another hash table for consensus clustering
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

    # run consensus clustering in parallel
    cons = foreach(i = 1:dim(all.combinations)[1], .packages = c("SC3", "cluster"), .options.RNG=1234) %dorng% {
        try({
            d <- res[res$distan %in% strsplit(all.combinations[i, 1], " ")[[1]] &
                         res$dim.red %in% strsplit(all.combinations[i, 2], " ")[[1]] &
                         res$k == as.numeric(all.combinations[i, 3]), ]

            dat <- consensus_clustering(d$labs)

            diss <- dist(dat)
            clust <- hclust(diss)
            clusts <- cutree(clust, k = as.numeric(all.combinations[i, 3]))

            silh <- silhouette(clusts, diss)

            labs <- NULL
            for(j in unique(clusts[clust$order])) {
                labs <- rbind(labs, paste(names(clusts[clusts == j]), collapse = " "))
            }

            labs <- as.data.frame(labs)
            colnames(labs) <- "Labels"

            return(list(dat, labs, clust, silh))
        })
    }

    # stop local cluster
    stopCluster(cl)

    output.param <- list(filename, distances, dimensionality.reductions,
                         cbind(all.combinations, cons),
                         dataset, study.dataset, svm.num.cells, cell.names,
                         study.cell.names, show.original.labels,
                         chisq.quantile)

    if(interactivity) {
        # start a shiny app in a browser window
        sc3_interactive(output.param)
    } else {
        saveRDS(output.param, paste0(filename, ".rds"))
    }
}

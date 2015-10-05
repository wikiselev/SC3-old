reindex_clusters <- function(ordering) {
    new.index <- NULL
    j <- 1
    for(i in unique(ordering)) {
        tmp <- rep(j, length(ordering[ordering == i]))
        names(tmp) <- names(ordering[ordering == i])
        new.index <- c(new.index, tmp)
        j <- j + 1
    }
    return(new.index)
}

de_gene_heatmap_param <- function(res) {
    row.ann <- data.frame("minus.log10.p.value" = -log10(res))
    rownames(row.ann) <- names(res)

    return(list(row.ann = row.ann))
}

mark_gene_heatmap_param <- function(mark.res, labs) {
    mark.res.plot <- NULL
    for(i in labs) {
        tmp <- mark.res[mark.res[,2] == i, ]
        if(dim(tmp)[1] > 10) {
            mark.res.plot <- rbind(mark.res.plot, tmp[1:10, ])
        }
    }

    row.ann <- data.frame(Cluster = factor(mark.res.plot$clusts, levels = unique(mark.res.plot$clusts)))
    rownames(row.ann) <- rownames(mark.res.plot)

    row.gaps <- as.numeric(mark.res.plot$clusts)
    row.gaps <- which(diff(row.gaps) != 0)

    return(list(mark.res.plot = mark.res.plot, row.ann = row.ann, row.gaps = row.gaps))
}

outl_cells_main <- function(d) {
    outl.res <- list()
    for(i in unique(colnames(d))) {
        # reduce p dimensions by using robust PCA
        t <- tryCatch({
            PcaHubert(d[ , colnames(d) == i])
        }, warning = function(cond) {
            message(cond)
        }, error = function(cond) {
            message(paste0("No outliers detected in cluster ", i, ". Distribution of gene expression in cells is too skewed towards 0."))
            return(NULL)
        })
        if(class(t) != "NULL") {
            # degrees of freedom used in mcd and chisquare distribution
            if(dim(t@loadings)[1] <= 6) {
                message(paste0("No outliers detected in cluster ", i, ". Small number of cells in the cluster."))
                out <- rep(0, dim(d[ , colnames(d) == i])[2])
                names(out) <- rep(i, dim(d[ , colnames(d) == i])[2])
                outl.res[[i]] <- out
            } else {
                df <- ifelse(dim(t@loadings)[2] > 3, 3, dim(t@loadings)[2])

                mcd <- NULL
                if(df != 1) {
                    mcd <- tryCatch({
                        covMcd(t@loadings[ , 1:df])
                    }, warning = function(cond) {
                        message(cond)
                    }, error = function(cond) {
                        message("No outliers detected in the cluster. Error in MCD.")
                        return(NULL)
                    })
                }

                if(class(mcd) != "NULL") {
                    # sqrt(mcd$mah) - sqrt of robust distance
                    # sqrt(qchisq(.95, df = length(mcd$best))) - sqrt of 97.5% quantile of a
                    # chi-squared distribution with p degrees of freedom
                    outliers <- sqrt(mcd$mah) - sqrt(qchisq(.9999, df = df))
                    outliers[which(outliers < 0)] <- 0
                    outl.res[[i]] <- outliers
                } else {
                    out <- rep(0, dim(d[ , colnames(d) == i])[2])
                    names(out) <- rep(i, dim(d[ , colnames(d) == i])[2])
                    outl.res[[i]] <- out
                }
            }
        } else {
            out <- rep(0, dim(d[ , colnames(d) == i])[2])
            names(out) <- rep(i, dim(d[ , colnames(d) == i])[2])
            outl.res[[i]] <- out
        }
    }

    nams <- NULL
    vals <- NULL
    for(i in 1:length(outl.res)) {
        vals <- c(vals, outl.res[[i]])
        nams <- c(nams, names(outl.res[[i]]))
    }
    names(vals) <- nams

    return(vals)
}
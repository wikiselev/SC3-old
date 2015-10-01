prepare_dataset <- function(dataset, clusts, cell.order, study.dataset) {
    colnames(dataset) <- clusts
    if(with_svm) {
        tmp <- cbind(dataset, study.dataset)
        cols <- colnames(tmp)
        inds <- NULL
        for(i in unique(clusts[cell.order])) {
            inds <- c(inds, which(cols == i))
        }
        return(tmp[ , inds])
    } else {
        return(dataset[ , cell.order])
    }
}

mark_genes_main <- function(d, cluster.order) {
    mark.res <<- get_marker_genes(d, as.numeric(colnames(d)))

    mark.res.plot <- NULL
    for(i in cluster.order) {
        tmp <- mark.res[mark.res[,2] == i, ]
        if(dim(tmp)[1] > 10) {
            mark.res.plot <- rbind(mark.res.plot, tmp[1:10, ])
        }
    }

    return(mark.res.plot)
}

mark_gene_heatmap_param <- function(d, mark.res.plot) {
    row.ann <- data.frame(Cluster = factor(mark.res.plot$clusts, levels = unique(mark.res.plot$clusts)))
    rownames(row.ann) <- rownames(mark.res.plot)

    row.gaps <- as.numeric(mark.res.plot$clusts)
    row.gaps <- which(diff(row.gaps) != 0)

    col.gaps <- as.numeric(colnames(d))
    col.gaps <- which(diff(col.gaps) != 0)

    return(list(row.ann = row.ann, row.gaps = row.gaps, col.gaps = col.gaps))
}

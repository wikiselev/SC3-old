getAUC <- function(gene, labels) {
    score <- rank(gene)
    # Get average score for each cluster
    ms <- aggregate(score ~ labels, FUN = mean)
    # Get cluster with highest average score
    posgroup <- ms[ms$score == max(ms$score), ]$labels
    # Return negatives if there is a tie for cluster with highest average score
    # (by definition this is not cluster specific)
    if(length(posgroup) > 1) {
        return (c(-1,-1,1))
    }
    # Create 1/0 vector of truths for predictions, cluster with highest average score vs everything else
    truth <- as.numeric(labels == posgroup)
    #Make predictions & get auc using RCOR package.
    pred <- prediction(score,truth)
    val <- unlist(performance(pred,"auc")@y.values)
    pval <- suppressWarnings(wilcox.test(score[truth == 1],score[truth == 0])$p.value)
    return(c(val,posgroup,pval))
}

get_marker_genes <- function(dataset, labels) {
    geneAUCs <- apply(dataset, 1, getAUC, labels = labels)
    geneAUCsdf <- data.frame(matrix(unlist(geneAUCs), nrow=length(geneAUCs)/3, byrow=T))
    rownames(geneAUCsdf) <- rownames(dataset)
    colnames(geneAUCsdf) <- c("AUC","clusts", "p.value")
    geneAUCsdf$AUC <- as.numeric(as.character(geneAUCsdf$AUC))
    geneAUCsdf$clusts <- as.numeric(as.character(geneAUCsdf$clusts))
    geneAUCsdf$p.value <- as.numeric(as.character(geneAUCsdf$p.value))

    geneAUCsdf$p.value <- p.adjust(geneAUCsdf$p.value)
    geneAUCsdf <- geneAUCsdf[geneAUCsdf$p.value < 0.01 & !is.na(geneAUCsdf$p.value), ]

    geneAUCsdf <- geneAUCsdf[geneAUCsdf$AUC > 0.85, ]

    d <- NULL
    for(i in sort(unique(geneAUCsdf$clusts))) {
        tmp <- geneAUCsdf[geneAUCsdf$clusts == i, ]
        tmp <- tmp[order(tmp$AUC, decreasing = T),]
        d <- rbind(d, tmp)
    }

    colnames(d) <- c("AUC","clusts","p.value")
    return(d)
}

kruskal_statistics <- function(dataset, labels) {
    t <- apply(dataset, 1, kruskal.test, g = factor(labels))
    ps <- unlist(lapply(t, "[[", "p.value"))
    ps <- p.adjust(ps)
    ps <- ps[!is.na(ps)]
    ps <- ps[ps < 0.05]
    return(ps[order(ps)])
}

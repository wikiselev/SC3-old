getAUC <- function(gene, labels) {
    score <- rank(gene)
    # Get average score for each cluster
    ms <- aggregate(score ~ labels, FUN = mean)
    # Get cluster with highest average score
    posgroup <- ms[ms$score == max(ms$score), ]$labels
    # Return negatives if there is a tie for cluster with highest average score
    # (by definition this is not cluster specific)
    if(length(posgroup) > 1) {
        return (c(-1,-1))
    }
    # Create 1/0 vector of truths for predictions, cluster with highest average score vs everything else
    truth <- as.numeric(labels == posgroup)
    #Make predictions & get auc using RCOR package.
    pred <- prediction(score,truth)
    val <- unlist(performance(pred,"auc")@y.values)
    return(c(val,posgroup))
}

get_marker_genes <- function(dataset, labels) {
    geneAUCs <- apply(dataset, 1, getAUC, labels = labels)
    geneAUCsdf <- data.frame(matrix(unlist(geneAUCs), nrow=length(geneAUCs)/2, byrow=T))
    rownames(geneAUCsdf) <- rownames(dataset)
    colnames(geneAUCsdf) <- c("AUC","Group")
    geneAUCsdf <- geneAUCsdf[geneAUCsdf$AUC > 0.8,]
    geneAUCsdf <- geneAUCsdf[order(geneAUCsdf[,1], decreasing = T),]
    return(geneAUCsdf)
}

kruskal_statistics <- function(dataset, labels) {
    t <- apply(dataset, 1, kruskal.test, g = factor(labels))
    ps <- unlist(lapply(t, "[[", "p.value"))
    ps <- p.adjust(ps, "bonferroni")
    ps <- ps[ps < 0.05]
    return(ps[order(ps)])
}

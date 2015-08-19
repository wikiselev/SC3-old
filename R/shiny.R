run_shiny_app <- function(filename, distances, dimensionality.reductions, cons.table, dataset, study.dataset, svm.num.cells, working.sample, study.sample) {

    dist.opts <- strsplit(unlist(cons.table[,1]), " ")
    dim.red.opts <- strsplit(unlist(cons.table[,2]), " ")

    distances <- as.list(distances)
    names(distances) <- distances

    dimensionality.reductions <- as.list(dimensionality.reductions)
    names(dimensionality.reductions) <- dimensionality.reductions

    colour.pallete <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(median(as.numeric(unlist(cons.table[,3]))))
    plot.width <- 600
    plot.height <- 600

    shinyApp(
        ui = fluidPage(
            headerPanel(
                HTML("SC<sup>3</sup> - Single-Cell Consensus Clustering")
            ),
            sidebarPanel(
                h4("1. Clustering"),
                sliderInput("clusters", label = "k",
                            min = min(as.numeric(unlist(cons.table[,3]))) + 1,
                            max = max(as.numeric(unlist(cons.table[,3]))),
                            value = min(as.numeric(unlist(cons.table[,3])))),

                checkboxGroupInput("distance", label = "Distance metrics",
                                   choices = distances,
                                   selected = distances[1]),

                checkboxGroupInput("dimRed", label = "Dimensionality reduction",
                                   choices = dimensionality.reductions,
                                   selected = dimensionality.reductions[1]),

                if(dim(study.dataset)[2] > 0) {
                    h4("1+. SVM")},
                if(dim(study.dataset)[2] > 0) {
                    p("Press this button when you have found the best clustering\n\n")},
                if(dim(study.dataset)[2] > 0) {
                    actionButton("svm", label = "Run SVM")},

                h4("2. Gene identificatiion"),
                p("\n\n"),
                actionButton("get_de_genes", label = "Get DE genes"),
                p("\n\n"),
                actionButton("get_mark_genes", label = "Get Marker genes"),

                h4("3. Save results"),
                p("\n\n"),
                downloadLink('labs', label = "Save cell labels"),
                p("\n\n"),
                downloadLink('markers', label = "Save cluster markers"),
                p("\n\n"),
                downloadLink('de', label = "Save de genes")
            ),
            mainPanel(
                uiOutput('mytabs')
            )
        ),
        server = function(input, output) {
            output$mytabs = renderUI({
                if(dim(study.dataset)[2] > 0) {
                    myTabs <- list(tabPanel("Consensus Matrix (1)", plotOutput('plot')),
                                   tabPanel("Expression Matrix (1)", plotOutput('matrix')),
                                   tabPanel("Cell Labels (1)", div(htmlOutput('labels'), style = "font-size:80%")),
                                   tabPanel("SVM (1+)", textOutput('svm_panel')),
                                   tabPanel("DE genes (2)", plotOutput('de_genes')),
                                   tabPanel("Marker genes (2)", plotOutput('mark_genes')))
                } else {
                    myTabs <- list(tabPanel("Consensus Matrix (1)", plotOutput('plot')),
                                   tabPanel("Expression Matrix (1)", plotOutput('matrix')),
                                   tabPanel("Cell Labels (1)", div(htmlOutput('labels'), style = "font-size:80%")),
                                   tabPanel("DE genes (2)", plotOutput('de_genes')),
                                   tabPanel("Marker genes (2)", plotOutput('mark_genes')))
                }
                do.call(tabsetPanel, myTabs)
            })
            get_consensus <- reactive({
                res <- cons.table[unlist(lapply(dist.opts, function(x){setequal(x, input$distance)})) &
                                      unlist(lapply(dim.red.opts, function(x){setequal(x, input$dimRed)})) &
                                      as.numeric(cons.table[ , 3]) == input$clusters, 4]
                return(res[[1]])
            })

            get_consensus_1 <- reactive({
                res <- cons.table[unlist(lapply(dist.opts, function(x){setequal(x, input$distance)})) &
                                      unlist(lapply(dim.red.opts, function(x){setequal(x, input$dimRed)})) &
                                      as.numeric(cons.table[ , 3]) == (input$clusters - 1), 4]
                return(res[[1]])
            })

            output$plot <- renderPlot({
                d <- get_consensus()
                hc <- d[[3]]
                d <- d[[1]]
                show_labs <- TRUE
                if(dim(d)[1] > 80) {
                    show_labs <- FALSE
                }
                withProgress(message = 'Plotting...', value = 0, {
                    pheatmap(d,
                             color = colour.pallete,
                             cluster_rows = hc,
                             cutree_rows = input$clusters, cutree_cols = input$clusters,
                             show_rownames = show_labs, show_colnames = show_labs)
                })
            }, height = plot.height, width = plot.width)

            output$matrix <- renderPlot({
                hc <- get_consensus()[[3]]
                withProgress(message = 'Plotting...', value = 0, {
                    pheatmap(dataset, cluster_cols = hc,
                             cutree_cols = input$clusters,
                             color = colour.pallete,
                             kmeans_k = 100, show_rownames = F, show_colnames = F,
                             treeheight_col = 0)
                })
            }, height = plot.height, width = plot.width)

            get_de_genes <- eventReactive(input$get_de_genes, {
                validate(
                    need(try(!is.null(rownames(dataset))), "\nNo gene names provided in the input expression matrix!")
                )
                withProgress(message = 'Calculating DE genes...', value = 0, {
                    hc <- get_consensus()[[3]]
                    clusts <- cutree(hc, input$clusters)
                    if(dim(study.dataset)[2] > 0) {
                        d <- cbind(dataset, study.dataset)
                        colnames(d) <- c(clusts, colnames(study.dataset))
                    } else {
                        d <- dataset
                        colnames(d) <- clusts
                    }

                    res <- kruskal_statistics(d, colnames(d))
                    de.res <<- res

                    res <- head(res, 68)
                    d <- d[names(res), ]

                    p.value.ann <- split(res, ceiling(seq_along(res)/17))
                    p.value.ranges <- as.vector(unlist(lapply(p.value.ann, function(x){rep(max(x), length(x))})))
                    p.value.ranges <- format(p.value.ranges, scientific = T, digits = 2)
                    p.value.ranges <- paste("<", p.value.ranges, sep=" ")

                    p.value.ann <- data.frame(p.value = factor(p.value.ranges, levels = unique(p.value.ranges)))
                    rownames(p.value.ann) <- names(res)

                    if(dim(study.dataset)[2] > 0) {
                        col.gaps <- as.numeric(colnames(d))
                        col.gaps <- col.gaps[order(col.gaps)]
                        col.gaps <- which(diff(col.gaps) != 0)
                        pheatmap(d[, order(colnames(d))], color = colour.pallete,
                                 show_colnames = F,
                                 cluster_rows = F, cluster_cols = F, annotation_row = p.value.ann,
                                 annotation_names_row = F,
                                 treeheight_col = 0, gaps_col = col.gaps)
                    } else {
                        pheatmap(d,
                                 color = colour.pallete, show_colnames = F,
                                 cluster_cols = hc,
                                 cutree_cols = input$clusters, cluster_rows = F,
                                 annotation_row = p.value.ann, annotation_names_row = F,
                                 treeheight_col = 0)
                    }
                })

            })

            get_mark_genes <- eventReactive(input$get_mark_genes, {
                validate(
                    need(try(!is.null(rownames(dataset))), "\nNo gene names provided in the input expression matrix!")
                )
                withProgress(message = 'Calculating Marker genes...', value = 0, {
                    hc <- get_consensus()[[3]]
                    clusts <- cutree(hc, input$clusters)

                    if(dim(study.dataset)[2] > 0) {
                        d <- cbind(dataset, study.dataset)
                        colnames(d) <- c(clusts, colnames(study.dataset))
                    } else {
                        d <- dataset
                        colnames(d) <- clusts
                    }

                    res <- get_marker_genes(d, as.numeric(colnames(d)))
                    mark.res <<- res
                    colnames(mark.res) <<- c("AUC", "clusts")

                    res1 <- NULL
                    for(i in unique(res$Group)) {
                        tmp <- res[res[,2] == i, ]
                        if(dim(tmp)[1] > 10) {
                            tmp <- tmp[1:10, ]
                        }
                        res1 <- rbind(res1, tmp)
                    }

                    d <- d[rownames(res1), ]
                    # colnames(d) <- names(clusts)

                    row.ann <- data.frame(Cluster = factor(res1$Group, levels = unique(res1$Group)))
                    rownames(row.ann) <- rownames(res1)

                    # col.ann <- data.frame(Cluster1 = factor(clusts, levels = unique(clusts[hc$order])))
                    # rownames(col.ann) <- clusts

                    # col.labs <- clusts

                    row.gaps <- res1$Group
                    row.gaps <- which(diff(row.gaps) != 0)

                    if(dim(study.dataset)[2] > 0) {
                        col.gaps <- as.numeric(colnames(d))
                        col.gaps <- col.gaps[order(col.gaps)]
                        col.gaps <- which(diff(col.gaps) != 0)
                        pheatmap(d[, order(colnames(d))], color = colour.pallete,
                                 show_colnames = F,
                                 cluster_rows = F, cluster_cols = F, annotation_row = row.ann,
                                 annotation_names_row = F,
                                 treeheight_col = 0, gaps_row = row.gaps, gaps_col = col.gaps)
                    } else {
                        pheatmap(d,
                                 color = colour.pallete, show_colnames = F,
                                 cluster_cols = hc,
                                 cutree_cols = input$clusters, cluster_rows = F,
                                 annotation_row = row.ann,
                                 annotation_names_row = F,
                                 treeheight_col = 0, gaps_row = row.gaps)
                    }
                })
            })

            get_svm <- eventReactive(input$svm, {
                withProgress(message = 'Running SVM...', value = 0, {
                    hc <- get_consensus()[[3]]
                    clusts <- cutree(hc, input$clusters)
                    colnames(dataset) <- clusts
                    prediction <- support_vector_machines1(dataset, study.dataset, "linear")
                    colnames(study.dataset) <<- prediction
                    return(prediction)
                })
            })

            output$svm_panel <- renderText({
                svm.prediction <<- get_svm()
                "SVM finished!"
                # c("SVM finished!", "\n\n", svm.prediction)
            })

            output$de_genes <- renderPlot({
                get_de_genes()
            }, height = plot.height, width = plot.width)

            output$mark_genes <- renderPlot({
                get_mark_genes()
            }, height = plot.height, width = plot.width)

            output$labels <- renderUI({
                d <- get_consensus()[[2]]
                d1 <- get_consensus_1()[[2]]

                labs1 <- list()
                cols <- brewer.pal(input$clusters - 1, "Paired")
                for(i in 1:(input$clusters - 1)) {
                    col <- cols[i]
                    ind <- unlist(strsplit(as.character(d1[i, ]), " "))
                    for(j in ind) {
                        labs1[[j]] <- paste0("<font color=\"", col, "\">", j, "</font>")
                    }
                }
                labs <- "<br/>"
                for(i in 1:input$clusters) {
                    ind <- unlist(strsplit(as.character(d[i, ]), " "))
                    for(j in ind) {
                        labs <- c(labs, labs1[[j]])
                    }
                    labs <- c(labs, c("<br/>", "<hr>"))
                }

                HTML(paste0(labs))
            })

            output$labs <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-labels.csv")
                },
                content = function(file) {
                    hc <- get_consensus()[[3]]
                    clusts <- cutree(hc, k = input$clusters)
                    if(dim(study.dataset)[2] > 0) {
                        names(clusts) <- working.sample
                        names(svm.prediction) <- study.sample
                        clusts <- c(clusts, svm.prediction)
                        clusts <- clusts[order(as.numeric(names(clusts)))]
                    }
                    write.table(t(data.frame(clusts)), file = file)
                }
            )

            output$markers <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-markers.csv")
                },
                content = function(file) {
                    write.csv(mark.res, file = file)
                }
            )

            output$de <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-de-genes.csv")
                },
                content = function(file) {
                    write.csv(de.res, file = file)
                }
            )

        }
    )
}

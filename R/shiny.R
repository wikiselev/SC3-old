sc3_interactive <- function(input.param) {
    filename <- input.param[[1]]
    distances <- input.param[[2]]
    dimensionality.reductions <- input.param[[3]]
    cons.table <- input.param[[4]]
    dataset <- input.param[[5]]
    study.dataset <- input.param[[6]]
    svm.num.cells <- input.param[[7]]
    cell.names <- input.param[[8]]
    study.cell.names <- input.param[[9]]
    show.original.labels <- input.param[[10]]
    chisq.quantile <- input.param[[11]]

    ## define UI parameters
    dist.opts <- strsplit(unlist(cons.table[,1]), " ")
    dim.red.opts <- strsplit(unlist(cons.table[,2]), " ")

    distances <- as.list(distances)
    names(distances) <- distances

    dimensionality.reductions <- as.list(dimensionality.reductions)
    names(dimensionality.reductions) <- dimensionality.reductions

    colour.pallete <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(median(as.numeric(unlist(cons.table[,3]))))
    plot.width <- 800
    plot.height <- 800
    plot.height.small <- 400

    ## define server global variables

    values <- reactiveValues()

    if(dim(study.dataset)[2] > 0) {
        with_svm <- TRUE
        values$svm <- FALSE
        values$svm.ready <- FALSE
    } else {
        with_svm <- FALSE
    }

    shinyApp(
        ui = fluidPage(
            headerPanel(
                paste0("Clustering of ", filename)
            ),
            fluidRow(
                column(3,
                       wellPanel(
                            h4("1. Clustering"),
                            sliderInput("clusters", label = "Number of clusters k",
                                        min = min(as.numeric(unlist(cons.table[,3]))) + 1,
                                        max = max(as.numeric(unlist(cons.table[,3]))),
                                        value = median(as.numeric(unlist(cons.table[,3]))),
                                        step = 1,
                                        animate = animationOptions(interval = 2000, loop = F)),

                            checkboxGroupInput("distance", label = "Distance metrics",
                                               choices = distances,
                                               selected = distances),

                            checkboxGroupInput("dimRed", label = "Dimensionality reduction",
                                               choices = dimensionality.reductions,
                                               selected = dimensionality.reductions),

                            if(with_svm) {
                                h4("1+. SVM")},
                            if(with_svm) {
                                p("Press this button when you have found the best clustering\n\n")},
                            if(with_svm) {
                                actionButton("svm", label = "Run SVM")},

                            h4("2. Gene/Cell Analysis"),
                            p("\n\nOpen a corresponding panel first, then press a button:"),
                            actionButton("get_de_genes", label = "Get DE genes"),
                            p("\n\n"),
                            actionButton("get_mark_genes", label = "Get Marker genes"),
                            p("\n\n"),
                            actionButton("get_outliers", label = "Get Cells outliers"),

                            h4("3. GO Analysis"),
                            p("\n\nRun Marker genes analysis first."),
                            selectInput("cluster", "Choose a cluster:",
                                        c("None" = "NULL")),
                            actionButton("GO", label = "Go to Webgestalt"),
                            p(" (will open in Firefox)"),

                            h4("4. Save results"),
                            p("\n\n"),
                            downloadLink('labs', label = "Save cell labels"),
                            p("\n\n"),
                            downloadLink('de', label = "Save DE genes"),
                            p("\n\n"),
                            downloadLink('markers', label = "Save Marker genes"),
                            p("\n\n"),
                            downloadLink('outl', label = "Save cell outliers")
                       )
                ),
                column(9,
                       uiOutput('mytabs')
                )
            )
        ),
        server = function(input, output, session) {
            output$mytabs = renderUI({
                myTabs <- list(tabPanel("Consensus Matrix", plotOutput('consensus')),
                               tabPanel("Silhouette", plotOutput('silh')),
                               tabPanel("Cell Labels", div(htmlOutput('labels'), style = "font-size:80%")),
                               tabPanel("Expression Matrix", plotOutput('matrix')),
                               tabPanel("DE genes", plotOutput('de_genes')),
                               tabPanel("Marker genes", plotOutput('mark_genes')),
                               tabPanel("Cells outliers", plotOutput('outliers')))
                do.call(tabsetPanel, myTabs)
            })

            ## main reactive function for extraction of precalculated variables

            observe({
                res <- cons.table[unlist(lapply(dist.opts, function(x){setequal(x, input$distance)})) &
                                      unlist(lapply(dim.red.opts, function(x){setequal(x, input$dimRed)})) &
                                      as.numeric(cons.table[ , 3]) == input$clusters, 4][[1]]
                res1 <- cons.table[unlist(lapply(dist.opts, function(x){setequal(x, input$distance)})) &
                                      unlist(lapply(dim.red.opts, function(x){setequal(x, input$dimRed)})) &
                                      as.numeric(cons.table[ , 3]) == (input$clusters - 1), 4][[1]]

                values$consensus <- res[[1]]
                values$labels <- res[[2]]
                values$labels1 <- res1[[2]]
                values$hc <- res[[3]]
                values$silh <- res[[4]]

                clusts <- cutree(values$hc, input$clusters)
                cell.order <- order.dendrogram(as.dendrogram(values$hc))

                d <- dataset
                colnames(d) <- clusts
                d <- d[ , cell.order]
                values$original.labels <- cell.names[cell.order]
                values$new.labels <- reindex_clusters(colnames(d))
                colnames(d) <- values$new.labels
                values$dataset <- d

                values$col.gaps <- which(diff(as.numeric(colnames(d))) != 0)
            })

            observe({
                if(with_svm) {
                    values$svm.ready <- values$svm &
                        values$svm.clusters == paste(input$clusters, collapse = "_") &
                        values$svm.distance == paste(input$distance, collapse = "_") &
                        values$svm.dimRed == paste(input$dimRed, collapse = "_")
                }
            })

            observe({
                if(!is.null(values$mark.res)) {
                    clusts <- unique(colnames(values$dataset))
                    updateSelectInput(session, "cluster", choices = clusts)
                }

            })

            ## REACTIVE PANELS

            output$consensus <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    if(show.original.labels) {
                        ann <- data.frame(Input.Labels = factor(cell.names))
                        pheatmap(values$consensus,
                                 color = colour.pallete,
                                 cluster_rows = values$hc,
                                 cluster_cols = values$hc,
                                 cutree_rows = input$clusters,
                                 cutree_cols = input$clusters,
                                 annotation_col = ann,
                                 show_rownames = F,
                                 show_colnames = F)
                    } else {
                        pheatmap(values$consensus,
                                 color = colour.pallete,
                                 cluster_rows = values$hc,
                                 cluster_cols = values$hc,
                                 cutree_rows = input$clusters,
                                 cutree_cols = input$clusters,
                                 show_rownames = F,
                                 show_colnames = F)
                    }
                })
            }, height = plot.height, width = plot.width)

            output$silh <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    plot(values$silh, col = "black")
                })
            }, height = plot.height, width = plot.width)

            output$labels <- renderUI({
                labs1 <- list()
                cols <- iwanthue(input$clusters - 1)
                for(i in 1:(input$clusters - 1)) {
                    col <- cols[i]
                    ind <- unlist(strsplit(as.character(values$labels1[i, ]), " "))
                    for(j in ind) {
                        labs1[[j]] <- paste0("<font color=\"", col, "\">", j, "</font>")
                    }
                }

                labs <- paste0("<br/><font size=\"3\">Colours correspond to clusters obtained by clustering the data by <b>",
                               input$clusters - 1, "</b> clusters</font><br/>")
                labs <- c(labs, "<br/>")
                for(i in 1:input$clusters) {
                    ind <- unlist(strsplit(as.character(values$labels[i, ]), " "))
                    for(j in ind) {
                        labs <- c(labs, labs1[[j]])
                    }
                    labs <- c(labs, c("<br/>", "<hr>"))
                }

                HTML(paste0(labs))
            })

            output$matrix <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    if(show.original.labels) {
                        ann <- data.frame(Input.Labels = factor(values$original.labels))
                        t <- values$dataset
                        colnames(t) <- rownames(ann)
                        pheatmap(t,
                                 # color = colour.pallete,
                                 kmeans_k = 100,
                                 cluster_cols = F,
                                 show_rownames = F,
                                 show_colnames = F,
                                 annotation_col = ann,
                                 gaps_col = values$col.gaps)
#                                  main = "Expression matrix is log2 scaled and clustered in 100 clusters by kmeans.
#                                          Only the values of the cluster centers are shown.")
                    } else {
                        pheatmap(values$dataset,
                                 # color = colour.pallete,
                                 kmeans_k = 100,
                                 cluster_cols = F,
                                 show_rownames = F,
                                 show_colnames = F,
                                 gaps_col = values$col.gaps)
#                                  main = "Expression matrix is log2 scaled and clustered in 100 clusters by kmeans.
#                                  Only the values of the cluster centers are shown.")
                    }
                })
            }, height = plot.height, width = plot.width)

            ## REACTIVE BUTTONS

            observeEvent(input$svm, {
                withProgress(message = 'Running SVM...', value = 0, {
                    values$dataset.svm <- study.dataset

                    original.labels <- c(values$original.labels, study.cell.names)
                    colnames(values$dataset.svm) <- support_vector_machines1(values$dataset, study.dataset, "linear")

                    tmp <- cbind(values$dataset, values$dataset.svm)
                    cols <- colnames(tmp)
                    inds <- NULL
                    for(i in unique(colnames(tmp))) {
                        inds <- c(inds, which(cols == i))
                    }
                    tmp <- tmp[ , inds]

                    values$original.labels.svm <- original.labels[inds]
                    colnames(tmp) <- reindex_clusters(colnames(tmp))
                    values$new.labels.svm <- colnames(tmp)
                    values$dataset.svm <- tmp

                    values$col.gaps.svm <- which(diff(as.numeric(colnames(values$dataset.svm))) != 0)

                    values$svm <- TRUE
                    values$svm.clusters <- paste(input$clusters, collapse = "_")
                    values$svm.distance <- paste(input$distance, collapse = "_")
                    values$svm.dimRed <- paste(input$dimRed, collapse = "_")
                })
            })

            get_de_genes <- eventReactive(input$get_de_genes, {
                if(with_svm) {
                    validate(need(try(values$svm.ready), "\nPlease run SVM prediction first!"))
                }
                validate(
                    need(try(!is.null(rownames(dataset))), "\nNo gene names provided in the input expression matrix!")
                )
                withProgress(message = 'Calculating DE genes...', value = 0, {
                    # prepare dataset for plotting
                    if(with_svm) {
                        d <- values$dataset.svm
                        col.gaps <- values$col.gaps.svm
                    } else {
                        d <- values$dataset
                        col.gaps <- values$col.gaps
                    }
                    # define de genes
                    values$de.res <- kruskal_statistics(d, colnames(d))
                    # check the results of de_genes_main:
                    # global variable de.res
                    validate(
                        need(try(length(values$de.res) != 0), "\nUnable to find significantly (p-value < 0.05) differentially expressed genes from obtained clusters! Please try to change the number of clusters k and run DE analysis again.")
                    )

                    d.param <- de_gene_heatmap_param(head(values$de.res, 70))

                    pheatmap(d[names(head(values$de.res, 70)), ],
                             # color = colour.pallete,
                             show_colnames = F,
                             cluster_rows = F,
                             cluster_cols = F,
                             annotation_row = d.param$row.ann,
                             annotation_names_row = F,
                             gaps_col = col.gaps)
                })

            })

            get_mark_genes <- eventReactive(input$get_mark_genes, {
                if(with_svm) {
                    validate(need(try(values$svm.ready), "\nPlease run SVM prediction first!"))
                }
                validate(
                    need(try(!is.null(rownames(dataset))), "\nNo gene names provided in the input expression matrix!")
                )
                withProgress(message = 'Calculating Marker genes...', value = 0, {
                    # prepare dataset for plotting
                    if(with_svm) {
                        d <- values$dataset.svm
                        col.gaps <- values$col.gaps.svm
                    } else {
                        d <- values$dataset
                        col.gaps <- values$col.gaps
                    }
                    # define marker genes
                    values$mark.res <- get_marker_genes(d, as.numeric(colnames(d)))
                    # check the results of mark_genes_main:
                    # global variable mark.res
                    validate(
                        need(try(dim(values$mark.res)[1] != 0), "\nUnable to find significant marker genes from obtained clusters! Please try to change the number of clusters k and run marker analysis again.")
                    )
                    d.param <- mark_gene_heatmap_param(values$mark.res, unique(colnames(d)))
                    pheatmap(d[rownames(d.param$mark.res.plot), ],
                             # color = colour.pallete,
                             show_colnames = F,
                             cluster_rows = F,
                             cluster_cols = F,
                             annotation_row = d.param$row.ann,
                             annotation_names_row = F,
                             gaps_row = d.param$row.gaps,
                             gaps_col = col.gaps)
                })
            })

            get_outl <- eventReactive(input$get_outliers, {
                if(with_svm) {
                    validate(need(try(values$svm.ready), "\nPlease run SVM prediction first!"))
                }
                withProgress(message = 'Calculating cell outliers...', value = 0, {
                    # prepare dataset for plotting
                    if(with_svm) {
                        d <- values$dataset.svm
                    } else {
                        d <- values$dataset
                    }
                    # compute outlier cells
                    values$outl.res <- outl_cells_main(d, chisq.quantile)

                    t <- as.data.frame(values$outl.res)
                    colnames(t)[1] <- "outl"
                    t$Cluster <- names(values$outl.res)
                    t$Cells <- 1:dim(t)[1]
#                     figure4a <- t
#                     save(figure4a, file = "figure4a.rda")
                    t$Cluster <- factor(t$Cluster, levels = unique(as.character(sort(as.numeric(t$Cluster)))))
                    cols <- iwanthue(length(unique(t$Cluster)))
                    ggplot(t, aes(x = Cells, y = outl, fill = Cluster, color = Cluster)) +
                        geom_bar(stat = "identity") +
                        geom_point() +
                        scale_fill_manual(values = cols) +
                        scale_color_manual(values = cols) +
                        guides(color = FALSE, fill = FALSE) +
                        labs(y = "Outlier score") +
                        # coord_cartesian(xlim = c(0, length(t$Cells))) +
                        theme_bw()

#                     plot(values$outl.res,
#                          col = names(values$outl.res),
#                          type = "p", ylab = "Outliers", xlab = "Cells",
#                          pch = 16, cex = 1.1)
                })
            })

            observeEvent(input$GO, {
                validate(
                    need(try(!is.null(values$mark.res)), "\nPlease run marker genes analysis by clicking on \"Get Marker genes\" button!")
                )
                open_webgestalt_go(rownames(values$mark.res[values$mark.res[,2] == input$cluster, ]))
            })

            ## PANELS REACTIVE ON BUTTON CLICK

            output$de_genes <- renderPlot({
                get_de_genes()
            }, height = plot.height, width = plot.width)

            output$mark_genes <- renderPlot({
                get_mark_genes()
            }, height = plot.height, width = plot.width)

            output$outliers <- renderPlot({
                get_outl()
            }, height = plot.height.small, width = plot.width)

            ## DOWNLOAD LINKS

            output$labs <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-labels-", filename, ".xls")
                },

                content = function(file) {
                    if(with_svm) {
                        validate(need(try(values$svm.ready), "\nPlease run SVM prediction first!"))
                    }
                    if(with_svm) {
                        write.table(data.frame(new.labels = values$new.labels.svm, original.labels = values$original.labels.svm),
                                    file = file, row.names = F, quote = F, sep = "\t")}
                    else{
                        write.table(data.frame(new.labels = values$new.labels, original.labels = values$original.labels),
                                    file = file, row.names = F, quote = F, sep = "\t")
                    }
                }
            )

            output$de <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-de-genes-", filename, ".xls")
                },
                content = function(file) {
                    if(with_svm) {
                        validate(need(try(values$svm.ready), "\nPlease run SVM prediction first!"))
                    }
                    validate(
                        need(try(!is.null(values$de.res)), "\nPlease run differential expression analysis by clicking on \"Get DE genes\" button!")
                    )
                    write.table(data.frame(gene = names(values$de.res), p.value = values$de.res),
                                file = file, row.names = F, quote = F, sep = "\t")
                }
            )

            output$markers <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-markers-", filename, ".xls")
                },
                content = function(file) {
                    if(with_svm) {
                        validate(need(try(values$svm.ready), "\nPlease run SVM prediction first!"))
                    }
                    validate(
                        need(try(!is.null(values$mark.res)), "\nPlease run marker genes analysis by clicking on \"Get Marker genes\" button!")
                    )
                    write.table(data.frame(new.labels = values$mark.res[,2],
                                           gene = rownames(values$mark.res),
                                           AUROC = values$mark.res[,1],
                                           p.value = values$mark.res[,3]),
                                file = file, row.names = F, quote = F, sep = "\t")
                }
            )

            output$outl <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-outliers-", filename, ".xls")
                },
                content = function(file) {
                    if(with_svm) {
                        validate(need(try(values$svm.ready), "\nPlease run SVM prediction first!"))
                    }
                    validate(
                        need(try(!is.null(values$outl.res)), "\nPlease run cell outlier analysis by using \"Get Cells outliers\" button!")
                    )
                    if(with_svm) {
                        write.table(data.frame(new.labels = names(values$outl.res), original.labels = values$original.labels.svm, MCD.dist = values$outl.res),
                                    file = file, row.names = F, quote = F, sep = "\t")
                    } else {
                        write.table(data.frame(new.labels = names(values$outl.res), original.labels = values$original.labels, MCD.dist = values$outl.res),
                                    file = file, row.names = F, quote = F, sep = "\t")
                    }
                }
            )

            session$onSessionEnded(function() { stopApp() } )
        },
        options = list(launch.browser = T)
    )

}

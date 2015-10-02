run_shiny_app <- function(filename, distances, dimensionality.reductions, cons.table, dataset, study.dataset, svm.num.cells, working.sample, study.sample, cell.names, study.cell.names) {

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
    plot.height.small <- 300

    ## define server global variables

    values <- reactiveValues()

    if(dim(study.dataset)[2] > 0) {
        with_svm <- TRUE
    } else {
        with_svm <- FALSE
    }

    shinyApp(
        ui = fluidPage(
            headerPanel(
                paste0("Clustering of ", filename)
            ),
            sidebarPanel(
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

                h4("2. Analysis"),
                p("\n\n"),
                actionButton("get_de_genes", label = "Get DE genes"),
                p("\n\n"),
                actionButton("get_mark_genes", label = "Get Marker genes"),
                p("\n\n"),
                actionButton("get_outliers", label = "Get Cells outliers"),

                h4("3. Save results"),
                p("\n\n"),
                downloadLink('labs', label = "Save cell labels"),
                p("\n\n"),
                downloadLink('de', label = "Save DE genes"),
                p("\n\n"),
                downloadLink('markers', label = "Save cluster markers"),
                p("\n\n"),
                downloadLink('outl', label = "Save cell outliers")
            ),
            mainPanel(
                uiOutput('mytabs')
            )
        ),
        server = function(input, output, session) {
            output$mytabs = renderUI({
                if(with_svm) {
                    myTabs <- list(tabPanel("Consensus Matrix (1)", plotOutput('consensus')),
                                   tabPanel("Silhouette (1)", plotOutput('silh')),
                                   tabPanel("Cell Labels (1)", div(htmlOutput('labels'), style = "font-size:80%")),
                                   tabPanel("Expression Matrix (1)", plotOutput('matrix')),
                                   tabPanel("SVM (1+)", textOutput('svm_panel')),
                                   tabPanel("DE genes (2)", plotOutput('de_genes')),
                                   tabPanel("Marker genes (2)", plotOutput('mark_genes')),
                                   tabPanel("Cells outliers (2)", plotOutput('outliers')))
                } else {
                    myTabs <- list(tabPanel("Consensus Matrix (1)", plotOutput('consensus')),
                                   tabPanel("Silhouette (1)", plotOutput('silh')),
                                   tabPanel("Cell Labels (1)", div(htmlOutput('labels'), style = "font-size:80%")),
                                   tabPanel("Expression Matrix (1)", plotOutput('matrix')),
                                   tabPanel("DE genes (2)", plotOutput('de_genes')),
                                   tabPanel("Marker genes (2)", plotOutput('mark_genes')),
                                   tabPanel("Cells outliers (2)", plotOutput('outliers')))
                }
                do.call(tabsetPanel, myTabs)
            })

            ## main reactive function for extraction of precalculated variables

            update_clustering <- observe({
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

            ## REACTIVE PANELS

            output$consensus <- renderPlot({
                withProgress(message = 'Plotting...', value = 0, {
                    pheatmap(values$consensus,
                             color = colour.pallete,
                             cluster_rows = values$hc,
                             cutree_rows = input$clusters,
                             cutree_cols = input$clusters)
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
                    pheatmap(values$dataset,
                             color = colour.pallete,
                             kmeans_k = 100,
                             cluster_cols = F,
                             show_rownames = F,
                             show_colnames = F,
                             gaps_col = values$col.gaps,
                             main = "Expression matrix is log2 scaled and clustered in 100 clusters by kmeans.
                                     Only the values of the cluster centers are shown.")
                })
            }, height = plot.height, width = plot.width)

            ## REACTIVE BUTTONS

            get_svm <- eventReactive(input$svm, {
                withProgress(message = 'Running SVM...', value = 0, {
                    prediction <- support_vector_machines1(values$dataset, study.dataset, "linear")
                    return(prediction)
                })
            })

            get_de_genes <- eventReactive(input$get_de_genes, {
                validate(
                    need(try(!is.null(rownames(dataset))), "\nNo gene names provided in the input expression matrix!")
                )
                withProgress(message = 'Calculating DE genes...', value = 0, {
                    # prepare dataset for plotting
                    d <- values$dataset
                    if(with_svm) {
                        d <- prepare_dataset(d, study.dataset)
                    }
                    # define de genes
                    values$de.res <- kruskal_statistics(d, colnames(d))
                    # check the results of de_genes_main:
                    # global variable de.res
                    validate(
                        need(try(length(values$de.res) != 0), "\nUnable to find significantly (p-value < 0.05) differentially expressed genes from obtained clusters! Please try to change the number of clusters k and run DE analysis again.")
                    )

                    d.param <- de_gene_heatmap_param(head(values$de.res, 70))

                    pheatmap(d[names(head(values$de.res, 70)), ], color = colour.pallete,
                             show_colnames = F,
                             cluster_rows = F,
                             cluster_cols = F,
                             annotation_row = d.param$row.ann,
                             annotation_names_row = F,
                             gaps_col = values$col.gaps)
                })

            })

            get_mark_genes <- eventReactive(input$get_mark_genes, {
                validate(
                    need(try(!is.null(rownames(dataset))), "\nNo gene names provided in the input expression matrix!")
                )
                withProgress(message = 'Calculating Marker genes...', value = 0, {
                    # prepare dataset for plotting
                    d <- values$dataset
                    if(with_svm) {
                        d <- prepare_dataset(d, study.dataset)
                    }
                    # define marker genes
                    values$mark.res <- get_marker_genes(d, as.numeric(colnames(d)))
                    # check the results of mark_genes_main:
                    # global variable mark.res
                    validate(
                        need(try(dim(values$mark.res)[1] != 0), "\nUnable to find significant marker genes from obtained clusters! Please try to change the number of clusters k and run marker analysis again.")
                    )
                    d.param <- mark_gene_heatmap_param(values$mark.res, unique(colnames(d)))
                    pheatmap(d[rownames(d.param$mark.res.plot), ], color = colour.pallete,
                             show_colnames = F,
                             cluster_rows = F,
                             cluster_cols = F,
                             annotation_row = d.param$row.ann,
                             annotation_names_row = F,
                             gaps_row = d.param$row.gaps,
                             gaps_col = values$col.gaps)
                })
            })

            get_outl <- eventReactive(input$get_outliers, {
                withProgress(message = 'Calculating cell outliers...', value = 0, {
                    # prepare dataset for plotting
                    d <- values$dataset
                    if(with_svm) {
                        d <- prepare_dataset(d, study.dataset)
                    }

                    # compute outlier cells
                    values$outl.res <- unlist(outl_cells_main(d))

                    print(values$outl.res)

                    plot(values$outl.res,
                         col = names(values$outl.res),
                         type = "p", ylab = "Outliers", xlab = "Cells",
                         pch = 16, cex = 1.1)
                })
            })

            ## PANELS REACTIVE ON BUTTON CLICK

            output$svm_panel <- renderText({
                svm.prediction <- get_svm()
                "SVM finished!"
            })

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
                    paste0("k-", input$clusters, "-labels-", filename, ".csv")
                },
                content = function(file) {
                    write.table(data.frame(new.labels = values$new.labels, original.labels = values$original.labels),
                                file = file, row.names = F, quote = F, sep = "\t")
                }
            )

            output$de <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-de-genes-", filename, ".csv")
                },
                content = function(file) {
                    validate(
                        need(try(!is.null(values$de.res)), "\nPlease first run differential expression analysis by using \"Get DE genes\" button!")
                    )
                    write.table(data.frame(gene = names(values$de.res), p.value = values$de.res),
                                file = file, row.names = F, quote = F, sep = "\t")
                }
            )

            output$markers <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-markers-", filename, ".csv")
                },
                content = function(file) {
                    validate(
                        need(try(!is.null(values$mark.res)), "\nPlease first run marker genes analysis by using \"Get Marker genes\" button!")
                    )
                    write.table(data.frame(gene = rownames(values$mark.res),
                                           AUC = values$mark.res[,1],
                                           cluster = values$mark.res[,2]),
                                file = file, row.names = F, quote = F, sep = "\t")
                }
            )

            output$outl <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-outliers-", filename, ".csv")
                },
                content = function(file) {
                    validate(
                        need(try(!is.null(values$outl.res)), "\nPlease first run marker genes analysis by using \"Get Marker genes\" button!")
                    )
                    write.table(data.frame(cluster = names(values$outl.res), MCD.dist = values$outl.res),
                                file = file, row.names = F, quote = F, sep = "\t")
                }
            )

            session$onSessionEnded(function() { stopApp() } )
        },
        options = list(launch.browser = T)
    )

}

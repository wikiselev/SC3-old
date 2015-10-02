run_shiny_app <- function(filename, distances, dimensionality.reductions, cons.table, dataset, study.dataset, svm.num.cells, working.sample, study.sample, cell.names, study.cell.names) {

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

    cell.names <<- cell.names
    study.cell.names <<- study.cell.names
    original.labels <<- NULL
    new.labels <<- NULL

    de.res <<- NULL
    mark.res <<- NULL
    outl.res <<- list()

    if(dim(study.dataset)[2] > 0) {
        with_svm <<- TRUE
    } else {
        with_svm <<- FALSE
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

            ## REACTIVE PANELS

            update_clustering <- reactive({
                res <- cons.table[unlist(lapply(dist.opts, function(x){setequal(x, input$distance)})) &
                                      unlist(lapply(dim.red.opts, function(x){setequal(x, input$dimRed)})) &
                                      as.numeric(cons.table[ , 3]) == input$clusters, 4][[1]]
                res1 <- cons.table[unlist(lapply(dist.opts, function(x){setequal(x, input$distance)})) &
                                      unlist(lapply(dim.red.opts, function(x){setequal(x, input$dimRed)})) &
                                      as.numeric(cons.table[ , 3]) == (input$clusters - 1), 4][[1]]
                input.consensus <<- res[[1]]
                input.labels <<- res[[2]]
                input.labels1 <<- res1[[2]]
                input.hc <<- res[[3]]
                input.clusts <<- cutree(input.hc, input$clusters)
                input.cell.order <<- order.dendrogram(as.dendrogram(input.hc))
                input.silh <<- res[[4]]
                cell.names <<- cell.names[input.cell.order]
                study.cell.names <<- study.cell.names
                colnames(dataset) <<- input.clusts
                dataset <<- dataset[ , input.cell.order]

                original.labels <<- cell.names
                new.labels <<- reindex_clusters(colnames(dataset))

                col.gaps <<- as.numeric(colnames(dataset))
                col.gaps <<- which(diff(col.gaps) != 0)
            })

            output$consensus <- renderPlot({
                update_clustering()
                withProgress(message = 'Plotting...', value = 0, {
                    pheatmap(input.consensus,
                             color = colour.pallete,
                             cluster_rows = input.hc,
                             cutree_rows = input$clusters,
                             cutree_cols = input$clusters)
                })
            }, height = plot.height, width = plot.width)

            output$silh <- renderPlot({
                update_clustering()
                withProgress(message = 'Plotting...', value = 0, {
                    plot(input.silh, col = "black")
                })
            }, height = plot.height, width = plot.width)

            output$labels <- renderUI({
                update_clustering()

                labs1 <- list()
                cols <- iwanthue(input$clusters - 1)
                for(i in 1:(input$clusters - 1)) {
                    col <- cols[i]
                    ind <- unlist(strsplit(as.character(input.labels1[i, ]), " "))
                    for(j in ind) {
                        labs1[[j]] <- paste0("<font color=\"", col, "\">", j, "</font>")
                    }
                }

                labs <- paste0("<br/><font size=\"3\">Colours correspond to clusters obtained by clustering the data by <b>",
                               input$clusters - 1, "</b> clusters</font><br/>")
                labs <- c(labs, "<br/>")
                for(i in 1:input$clusters) {
                    ind <- unlist(strsplit(as.character(input.labels[i, ]), " "))
                    for(j in ind) {
                        labs <- c(labs, labs1[[j]])
                    }
                    labs <- c(labs, c("<br/>", "<hr>"))
                }

                HTML(paste0(labs))
            })

            output$matrix <- renderPlot({
                update_clustering()
                withProgress(message = 'Plotting...', value = 0, {
                    pheatmap(dataset,
                             color = colour.pallete,
                             kmeans_k = 100,
                             cluster_cols = F,
                             show_rownames = F,
                             show_colnames = F,
                             gaps_col = col.gaps,
                             main = "Expression matrix is log2 scaled and clustered in 100 clusters by kmeans.
                                     Only the values of the cluster centers are shown.")
                })
            }, height = plot.height, width = plot.width)

            ## REACTIVE BUTTONS

            get_svm <- eventReactive(input$svm, {
                withProgress(message = 'Running SVM...', value = 0, {
                    update_clustering()
                    prediction <- support_vector_machines1(dataset, study.dataset, "linear")
                    colnames(study.dataset) <<- prediction
                    return(prediction)
                })
            })

            get_de_genes <- eventReactive(input$get_de_genes, {
                validate(
                    need(try(!is.null(rownames(dataset))), "\nNo gene names provided in the input expression matrix!")
                )
                withProgress(message = 'Calculating DE genes...', value = 0, {
                    update_clustering()
                    # prepare dataset for plotting
                    d <- prepare_dataset(dataset, study.dataset)
                    # define de genes
                    de.res <<- kruskal_statistics(d, colnames(d))
                    # check the results of de_genes_main:
                    # global variable de.res
                    validate(
                        need(try(length(de.res) != 0), "\nUnable to find significantly (p-value < 0.05) differentially expressed genes from obtained clusters! Please try to change the number of clusters k and run DE analysis again.")
                    )

                    d.param <- de_gene_heatmap_param(head(de.res, 70))

                    pheatmap(d[names(head(de.res, 70)), ], color = colour.pallete,
                             show_colnames = F,
                             cluster_rows = F,
                             cluster_cols = F,
                             annotation_row = d.param$row.ann,
                             annotation_names_row = F,
                             gaps_col = col.gaps)
                })

            })

            get_mark_genes <- eventReactive(input$get_mark_genes, {
                validate(
                    need(try(!is.null(rownames(dataset))), "\nNo gene names provided in the input expression matrix!")
                )
                withProgress(message = 'Calculating Marker genes...', value = 0, {
                    update_clustering()
                    # prepare dataset for plotting
                    d <- prepare_dataset(dataset, study.dataset)
                    # define marker genes
                    mark.res.plot <- mark_genes_main(d)
                    # check the results of mark_genes_main:
                    # global variable mark.res
                    validate(
                        need(try(dim(mark.res)[1] != 0), "\nUnable to find significant marker genes from obtained clusters! Please try to change the number of clusters k and run marker analysis again.")
                    )
                    d.param <- mark_gene_heatmap_param(mark.res.plot)
                    pheatmap(d[rownames(mark.res.plot), ], color = colour.pallete,
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
                withProgress(message = 'Calculating cell outliers...', value = 0, {
                    update_clustering()
                    # prepare dataset for plotting
                    d <- prepare_dataset(dataset, study.dataset)

                    # compute outlier cells
                    outl_cells_main(d)

                    plot(unlist(outl.res),
                         col = names(unlist(outl.res)),
                         type = "p", ylab = "Outliers", xlab = "Cells",
                         pch = 16, cex = 1.1)

                })

            })

            ## PANELS REACTIVE ON BUTTON CLICK

            output$svm_panel <- renderText({
                svm.prediction <<- get_svm()
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
                    out <- data.frame(new.labels = new.labels, original.labels = original.labels)
                    write.table(out, file = file, row.names = F, quote = F, sep = "\t")
                }
            )

            output$de <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-de-genes-", filename, ".csv")
                },
                content = function(file) {
                    validate(
                        need(try(!is.null(de.res)), "\nPlease first run differential expression analysis by using \"Get DE genes\" button!")
                    )
                    nams <- names(de.res)
                    names(de.res) <- NULL
                    out <- data.frame(gene = nams, p.value = de.res)
                    write.table(out, file = file, row.names = F, quote = F, sep = "\t")
                }
            )

            output$markers <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-markers-", filename, ".csv")
                },
                content = function(file) {
                    validate(
                        need(try(!is.null(mark.res)), "\nPlease first run marker genes analysis by using \"Get Marker genes\" button!")
                    )
                    out <- data.frame(gene = rownames(mark.res), AUC = mark.res[,1],
                                      cluster = mark.res[,2])
                    write.table(out, file = file, row.names = F, quote = F, sep = "\t")
                }
            )

            output$outl <- downloadHandler(
                filename = function() {
                    paste0("k-", input$clusters, "-outliers-", filename, ".csv")
                },
                content = function(file) {
                    validate(
                        need(try(!is.null(outl.res)), "\nPlease first run marker genes analysis by using \"Get Marker genes\" button!")
                    )
                    out <- data.frame(gene = rownames(mark.res), AUC = mark.res[,1],
                                      cluster = mark.res[,2])
                    write.table(out, file = file, row.names = F, quote = F, sep = "\t")
                }
            )

            session$onSessionEnded(function() { stopApp() } )
        },
        options = list(launch.browser = T)
    )

}

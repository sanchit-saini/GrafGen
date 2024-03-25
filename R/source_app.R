#                                                                              X
createApp <- function(obj, metadata=NULL, id=NULL) {

    # Check obj and metadata
    check_grafpop(obj) 
    check_metadata(metadata, id)

    test_results  <- obj$table
    test_metadata <- setup_metadata(metadata, id, test_results, 1)
    train_results <- getTrainResults()

    ui     <- app_ui()
    server <- app_server
    app    <- shinyApp(ui = ui, server = server)
    list(app=app, reference_results=train_results, 
         user_results=test_results, user_metadata=test_metadata)
}

app_ui <- function() {

    ui <- basicPage(
        tags$head(tags$style(HTML(".selectize-input {width:150px; height: 30px; font-size: 12px}"))),
        actionButton("quit", label = "Close App"),
        h4(strong("Data Options")),
        uiOutput("selected_type_ui"),
        uiOutput("selected_type_filter_ui"),
        uiOutput("selected_group_ui"),
        uiOutput("selected_group_filter_ui"),
        plotlyOutput("myplot",height="800px",width="1200px")
    )
}

app_server <- function(input, output, session) {

    test_metadata <- get("user_metadata")
    test_results  <- get("user_results")
    train_results <- get("reference_results")

    observeEvent(input$quit,{ stopApp() }) # stop app

    A_percent <- E_percent <- F_percent <- GD1_x <- GD2_y <- NULL
    Group_Var <- Nearest_neighbor <- Refpop <- Refpop_n <- NULL
    Sample <- Separation_percent <- Type_Var <- test <- NULL

    # from metadata file, selects all variables that isn't ID (first column)
    output$selected_type_ui <- renderUI({
        list_choices <- colnames(test_metadata)[-1]
        selectizeInput("selected_type",
        "Variable to filter / annotate by shape:", 
        choices=list_choices,  
        options = list(placeholder = 'Click to select a variable', 
            onInitialize = I('function() { this.setValue(""); }')))
    })

    # filter for variable to shape by
    output$selected_type_filter_ui <- renderUI({
        req(input$selected_type)
        selectizeInput("selected_type_filter","Variable values to include:", 
        choices=unique(sort(test_metadata[ ,input$selected_type])),  
        multiple=TRUE,selected=NULL,
        options = list( placeholder = 'Click to select one or more values'))
    })

    # from metadata file, selects all variables that isn't ID (first column)
    output$selected_group_ui <- renderUI({
        list_choices <- colnames(test_metadata)[-1]
        selectizeInput("selected_group", 
        "Additional variable to filter / display:", 
        choices=list_choices,  
        options = list(placeholder = 'Click to select a variable', 
        onInitialize = I('function() { this.setValue(""); }')))
    })

    # filter for variable to color by
    output$selected_group_filter_ui <- renderUI({
        req(input$selected_group)
        selectizeInput("selected_group_filter","Variable values to include:",
        choices=unique(sort(test_metadata[ ,input$selected_group])),  
        multiple=TRUE, selected=NULL, 
        options = list(placeholder = 'Click to select one or more values'))
    })

    output$type_legend_ui <- renderUI({
        req(input$selected_type!="")
        test_results_clean <- results_clean()$test_results_clean

        # how many rows will the legend print
        num_legend_rows <- ceiling(nlevels(test_results_clean$Type_Var)/5)
        plotOutput("type_legend",height=paste0(num_legend_rows*27.5,"px"),
        width="1100px")
    })

    # manipulates both test and train data
    results_clean <- reactive({
        req(length(input$selected_type)>0)
        selected_type  <- input$selected_type 
        selected_group <- input$selected_group

        # label type and group if both selected
        if (selected_type!="" & selected_group!=""){
            test_results_clean <- test_results %>% 
            inner_join(test_metadata,by=c("Sample"="ID")) %>% 
            as.data.frame()
            test_results_clean[,"Type_Var"] <- 
                test_results_clean[,selected_type]
            test_results_clean[,"Group_Var"] <- 
                test_results_clean[,selected_group]
            test_results_clean <- test_results_clean %>%  
            mutate(Type_Var = factor(Type_Var,
            levels=unique(sort(test_metadata[ ,selected_type]))),
            Group_Var = factor(Group_Var,
            levels=unique(sort(test_metadata[ ,selected_group]))),
            text = paste0(selected_type,", ",selected_group,", ",
            " Sample: ",Type_Var,", ",Group_Var,", ",Sample, "<br />",
            "Refpop, Neighbor, Separation: ", Refpop, ", ",
            Nearest_neighbor, ",", round(Separation_percent,0),
            "% ", "<br />",
            "African, European, Asian Ancestry: ", round(F_percent,0), 
            "%, ",round(E_percent,0), "%, ",round(A_percent,0),"%"))
        } else if (selected_type!="" & selected_group==""){
            # label type if selected
            test_results_clean <- test_results %>% 
            inner_join(test_metadata,by=c("Sample"="ID")) %>% 
            as.data.frame()
            test_results_clean[,"Type_Var"] <- 
                test_results_clean[,selected_type]
            test_results_clean <- test_results_clean %>%  
            mutate(Type_Var = factor(Type_Var,
            levels=unique(sort(test_metadata[ ,selected_type]))),
            text = paste0(selected_type,", "," Sample: ",Type_Var,", ",
            Sample, "<br />", "Refpop, Neighbor, Separation: ", Refpop, 
            ", ",Nearest_neighbor, ", ",round(Separation_percent,0),"% ", 
            "<br />", "African, European, Asian Ancestry: ", 
            paste0(round(F_percent,0), "%, ",round(E_percent,0), "%, ",
            round(A_percent,0),"%")))
        } else if (selected_type=="" & selected_group!=""){
            # label group if selected
            test_results_clean <- test_results %>% 
            inner_join(test_metadata,by=c("Sample"="ID")) %>% 
            as.data.frame()
            test_results_clean[,"Group_Var"] <- 
                test_results_clean[,selected_group]
            test_results_clean <- test_results_clean %>%  
            mutate(Group_Var = factor(Group_Var,
            levels=unique(sort(test_metadata[ ,selected_group]))),
            text = paste0(selected_group, ", " ," Sample: ", 
            Group_Var,", ",Sample, 
            "<br />", "Refpop, Neighbor, Separation: ", Refpop, ", ",
            Nearest_neighbor, ", ", round(Separation_percent,0),"% ", 
            "<br />", "African, European, Asian Ancestry: ", 
            paste0(round(F_percent,0), "%, ",round(E_percent,0), "%, ",
            round(A_percent,0),"%")))
        } else{
            test_results_clean <- test_results %>%
            mutate(text = paste0("Sample: ",Sample, "<br />",
            "Refpop, Neighbor, Separation: ", Refpop, ", ",
            Nearest_neighbor, ", ", round(Separation_percent,0),
            "% ", "<br />","African, European, Asian Ancestry: ", 
            paste0(round(F_percent,0), 
            "%, ",round(E_percent,0), "%, ",round(A_percent,0),"%")))
        }

        # calculate frequency
        refpop_n <- test_results_clean %>%
            group_by(Refpop) %>%
            reframe(n=n()) %>%
            mutate(Refpop_n = paste0(gsub("hpgp","",Refpop),"\n","n=",n)) %>%
            select(-n)

        refpop_order <- test_results_clean %>%
            left_join(refpop_n,by="Refpop") %>%
            arrange(GD1_x) %>%
            mutate(Refpop_n=factor(Refpop_n,levels=unique(Refpop_n)))

        test_results_clean <- test_results_clean %>%
            left_join(refpop_n,by="Refpop") %>%
            mutate(Refpop_n=factor(Refpop_n,
            levels=levels(refpop_order$Refpop_n)),
            Refpop=factor(Refpop,levels=unique(refpop_order$Refpop))) %>%
            as.data.frame()

        # adds sample size of test set annotated to refpop for trainset
        train_results_clean <- train_results %>%
            left_join(refpop_n,by="Refpop") %>%
            mutate(Refpop_n=factor(Refpop_n,
            levels=levels(refpop_order$Refpop_n)),
            Refpop=factor(Refpop,levels=unique(refpop_order$Refpop)))

        return(list(test_results_clean=test_results_clean,
            train_results_clean=train_results_clean))
    })

    # filter test data
    test_results_filtered <- reactive({
        req(nrow(results_clean()$test_results_clean)>0)
        selected_type <- input$selected_type
        selected_group <- input$selected_group
        selected_type_filter <- input$selected_type_filter
        selected_group_filter <- input$selected_group_filter
        test_results_clean <- results_clean()$test_results_clean
        test_results_filtered <- test_results_clean

        # filters for shape variable if selected
        if (length(selected_type_filter)>0){ 
            test_results_filtered <- test_results_filtered %>% 
            filter(Type_Var %in% selected_type_filter) 
        }

        # filters for second annotation variable if selected
        if (length(selected_group_filter)>0){ 
            test_results_filtered <- test_results_filtered %>% 
            filter(Group_Var %in% selected_group_filter) 
        }

        return(test_results_filtered)
    })

    # plot data
    output$myplot <- renderPlotly({
        test_results_clean <- results_clean()$test_results_clean
        test_results_filtered <- test_results_filtered()
        train_results_clean <- results_clean()$train_results_clean
        selected_type <- input$selected_type 
        selected_group <- input$selected_group

        # if no subjects available in filter, display message
        if (nrow(test_results_filtered)==0){
            output <- ggplot(test_results_filtered, aes(GD1_x,GD2_y))+
            theme_void()+
            geom_text(aes(x=1,y=1,
                label="No subjects available given current filter(s)"), 
                size=9)
            output <- ggplotly(output)
        } else {
            # --------------------------------------------
            # SCATTERPLOT + COLOR LEGEND
            # ------------------------------------------------
            mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(9) 
            p1 <- ggplot(train_results_clean, aes(GD1_x, GD2_y))+
            # theme stuff
            theme_classic()+
            theme(panel.background = element_blank(), 
            panel.border = element_rect(fill = NA, color = "black"),
            legend.title = element_blank(), legend.position="bottom",
            axis.text=element_text(face="bold",color="black",size=15),
            axis.title = element_text(face="bold", color="black",size=15))+
            guides(fill = guide_legend(nrow = 1), 
            color=guide_legend(nrow=1))+
            scale_fill_manual(values = mycolors)+
            scale_color_manual(values = mycolors)+
            scale_x_continuous("GD1",breaks=pretty_breaks(n=10),
            limits=c(min(c(1,min(test_results_filtered$GD1_x))),
            max(c(1.825,max(test_results_filtered$GD1_x)))))+
            scale_y_continuous("GD2",breaks=pretty_breaks(n=10),
            limits=c(min(c(1,min(test_results_filtered$GD2_y))),
            max(c(1.325,max(test_results_filtered$GD2_y)))))+
            stat_ellipse(aes(color=Refpop_n,fill=Refpop_n),
            geom = "polygon",alpha=0.25)+  # ellipse of training data
            # make triangle
            geom_segment(x = 1.05, y = 1.1, xend = 1.7658, yend = 1.1)+
            geom_segment(x = 1.05, y = 1.1, xend = 1.4701, yend = 1.2897)+
            geom_segment(x = 1.4701, y = 1.2897, xend = 1.7658, yend = 1.1)

            # make test datapoints a shape if that variable is selected, 
            #   otherwise default to circle
            if (selected_type!=""){
                # if number of categories is <= 15, make each category a shape
                if (length(unique(test_results_clean$Type_Var))<=15){
                    tvec <- c(1,3,4,0,2,5,6,7,9:13,18,35)
                    tmp  <- seq_len(nlevels(test_results_clean$Type_Var))
                    mapped_shapes <- 
                    data.frame(Type_Var = levels(test_results_clean$Type_Var)) %>%
                    mutate(Shape = tvec[tmp]) %>%
                    filter(Type_Var %in% unique(test_results_filtered$Type_Var)) %>%
                    mutate(Type_Var=factor(Type_Var,
                    levels=levels(test_results_filtered$Type_Var)))

                    mapped_shapes_list <- as.list(mapped_shapes$Shape)
                    names(mapped_shapes_list) <- mapped_shapes$Type_Var
                    req(row(mapped_shapes)==length(mapped_shapes_list))

                    p1 <- p1 + 
                    # add data points
                    geom_point(data=test_results_filtered,
                    aes(color=Refpop_n, shape=Type_Var),size=2.5) + 
                    # add dummy data invisible points just to plot 
                    #   the shapes itself in a legend
                    geom_point(data=mapped_shapes%>%mutate(GD1_x=1,GD2_y=1),
                    aes(shape=Type_Var),color="black")+
                    # dummy data white point to hide black point
                    geom_point(data=data.frame(GD1_x=1,GD2_y=1),size=3,
                    color="white")+
                    scale_shape_manual(values=mapped_shapes_list)
                } else {
                    # if number of categories is >15, make just as a circle
                    p1 <- p1 + geom_point(data=test_results_filtered,
                    aes(color=Refpop_n, shape=Type_Var),size=2.5) +
                    scale_shape_manual(values=
                    rep(1,length(unique((test_results_filtered$Type_Var)))))
                }
            } else if (selected_type==""){
                # make just as a circle if no type variable is selected
                p1 <- p1 + geom_point(data=test_results_filtered,
                aes(color=Refpop_n), shape=1,size=2.5)
            }

            output <- ggplotly(p1) %>% 
            layout(legend = list(orientation = "h",y = -0.15))
            # pop hood of plotly and modify
            output$x$layout$legend$title$text <- "" # remove legend title
            output$x$layout$legend$font$size <- 13 # change legend font size

            # fix legend names for colors. ellipses (1-9), 
            # change hover info to remove sample size
            for(i in seq_len(9)){
                output$x$data[[i]]$text <- 
                gsub("\\(","",str_split(output$x$data[[i]]$name,"\n")[[1]][1])
                output$x$data[[i]]$name <- 
                gsub("\\(","",str_split(output$x$data[[i]]$name,",")[[1]][1])  
            }

            # count number of levels under hood for actual data points 
            #   if number of categories > 15
            n_levels <- length(output$x$data)
            if (selected_type!=""){
                if (length(unique(test_results_clean$Type_Var))<=15){
                    n_levels <- length(output$x$data) - 
                    length(unique(test_results_filtered$Type_Var)) - 1

                    # fix legend names for shapes and 
                    #    remove hovers for dummy data
                    ivec <- (n_levels+1):length(output$x$data)
                    for(i in ivec){
                        output$x$data[[i]]$name <- gsub("\\(","",
                        str_split(output$x$data[[i]]$name,",")[[1]][1]) 
                        output$x$data[[i]]$hoverinfo <- "none"
                    }
                }
            }

            # for points, modify hover info. 
            # Different sets of point levels within each country
            counter <- 13 # starting level for point data
            jlen <- length(unique(test_results_filtered$Refpop))
            while (counter <= n_levels){
                for (j in seq_len(jlen)){
                    current_refpop <- paste0("hpgp",gsub("\\(","",
                    str_split(output$x$data[[counter]]$name,"\n")[[1]][1]))
                    current_df <- test_results_filtered %>% 
                    filter(Refpop==current_refpop)

                    if (selected_type!=""){
                        klen <- length(unique(current_df$Type_Var))
                        # manually define hover info for points 
                        #    for each type category
                        for (k in seq_len(klen)){
                            current_type <- 
                            str_split(output$x$data[[counter]]$text,
                            pattern="Type_Var: ")[[1]][2]
                            output$x$data[[counter]]$text <- 
                        (current_df %>% filter(Type_Var==current_type))$text
                            output$x$data[[counter]]$showlegend <- FALSE 
                            counter <- counter+1
                        }
                    } else{
                        output$x$data[[counter]]$text <- current_df$text
                        output$x$data[[counter]]$showlegend <- FALSE 
                        counter <- counter+1
                    }
                }
            }
        }
        output
    })
}

setup_metadata <- function(x, id, res, which) {

    ids  <- res[, 1, drop=TRUE]
    ret  <- data.frame("ID"=ids, stringsAsFactors=FALSE)

    if (length(x)) {
        cx <- colnames(x)
        if (!length(id)) id <- cx[1]
        xids <- x[, id, drop=TRUE]
        rows <- match(ids, xids)
        tmp  <- is.na(rows)
        tmp0 <- !tmp
        if (all(tmp)) {
            stop("None of the ids match in meta data and results table")
        }
        if (any(tmp)) {
           warning("Not all ids match in meta data and results table")
        }

        rows <- rows[tmp0]
        x    <- x[rows, , drop=FALSE] 
        tmp  <- !(cx %in% id)
        add  <- cx[tmp]
        if (length(add)) {
            for (v in add) {
                ret[, v]     <- "NA"
                ret[tmp0, v] <- x[, v, drop=TRUE] 
            } 
        }
    }
    if (which == 2) return(ret)
    if (ncol(ret) < 3) {
        cx  <- colnames(ret)
        add <- c("Group", "Type")
        tmp <- !(add %in% cx)
        add <- add[tmp]
        for (v in add) ret[, v] <- "NA"
    } 
    ret
}

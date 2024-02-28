#                                                                              X
interactiveReferencePlot <- function() {

    train_results <- getTrainResults()

    Refpop <- Sample <- Nearest_neighbor <- Separation_percent <- NULL
    F_percent <- E_percent <- A_percent <- NULL
    label1 <- label2 <- label3 <- GD1_x <- GD2_y <- Refpop_n <- .data <- NULL
    Country <- Country.code <- NULL

    train_results_clean <- train_results %>%
    mutate(label1 = paste0(Country,", ",Sample),
        label2 = paste0(Refpop,", ",Nearest_neighbor, ", ",
        round(Separation_percent,0),"%"),
        label3 = paste0(round(F_percent,0), "%, ",round(E_percent,0),
            "%, ",round(A_percent,0),"%")) %>%
    rename("Country, Sample"=label1,
          "Refpop, Neighbor, Separation"=label2,
          "African, European, Asian Ancestry"=label3)

    # calculate frequency
    refpop_n <- train_results_clean %>%
        group_by(Refpop) %>%
        reframe(n=n()) %>%
        mutate(Refpop_n = paste0(gsub("hpgp","",Refpop),"\n","n=",n)) %>%
        select(-n)

    refpop_order <- train_results_clean %>%
    left_join(refpop_n,by="Refpop") %>%
    arrange(GD1_x) %>%
    mutate(Refpop_n=factor(Refpop_n,levels=unique(Refpop_n)))

    train_results_clean <- train_results_clean %>%
    left_join(refpop_n,by="Refpop") %>%
    mutate(Refpop_n=factor(Refpop_n,levels=levels(refpop_order$Refpop_n)),
        Refpop=factor(Refpop,levels=unique(refpop_order$Refpop)))

    nb.cols <- 9
    mycolors_fill <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

    # Define the number of colors you want. 9 groups, 9 colors
    nb.cols <- 60
    mycolors <- colorRampPalette(brewer.pal(60, "Set1"))(nb.cols)

    train_results_filtered <- train_results_clean # apply no filters

    p1 <- ggplot(train_results_clean, aes(GD1_x, GD2_y)) + theme_classic()+
        theme(panel.background = element_blank(), 
            panel.border = element_rect(fill =   NA, color = "black"),
            legend.title = element_blank(), 
            legend.position="bottom", 
            axis.text=element_text(face="bold", family="serif", 
                color="black",size=15),
            axis.title = element_text(face="bold", family="serif", 
                color="black", size=15))+
            guides(fill = guide_legend(nrow = 3), 
            color=guide_legend(nrow=3))+
            scale_fill_manual(values = mycolors_fill)+
            scale_color_manual(values = mycolors)+
            scale_x_continuous("GD1",breaks=pretty_breaks(n=10), 
                limits=c(1,1.8))+
            scale_y_continuous("GD2",breaks=pretty_breaks(n=10), 
                limits=c(1,1.4))+
            # ellipse of training data
            stat_ellipse(aes(fill=Refpop_n, label4=Refpop),geom = "polygon",
                alpha=0.25)+
            # make triangle
            geom_segment(x = 1.05, y = 1.1, xend = 1.7658, yend = 1.1)+
            geom_segment(x = 1.05, y = 1.1, xend = 1.4701, yend = 1.2897)+
            geom_segment(x = 1.4701, y = 1.2897, xend = 1.7658, yend = 1.1)+
            # data points
            geom_point(data=train_results_filtered, 
                aes(label1=.data[["Country, Sample"]],
                label2=.data[["Refpop, Neighbor, Separation"]],
                label3=.data[["African, European, Asian Ancestry"]],
                color=Country.code),size=0.7)

    output <- ggplotly(p1, tooltip=c("label1"="Country, Sample",
                               "label2"="Refpop, Neighbor, Separation",
                               "label3"="African, European, Asian Ancestry",
                               "label4"="Refpop")) %>%
    layout(legend = list(orientation = "h",y = -0.15)) %>% 
    config(showTips = FALSE) 

    output$x$layout$legend$title$text <- "" # remove legend title
    output$x$layout$legend$font$size <- 12 # change legend font size

    for(i in seq_len(length(output$x$data))){
        if (!is.null(output$x$data[[i]]$name)){ 
            output$x$data[[i]]$name <- 
            gsub("\\(","",strsplit(output$x$data[[i]]$name,",")[[1]][1]) 
        }

        # adds spaces for hover variables
        if (!is.null(output$x$data[[i]]$text)){ 
           output$x$data[[i]]$text <- gsub("_"," ",output$x$data[[i]]$text) 
        }
    }

    # for ellipses, change hover info to remove sample size
    ilen <- length(unique(train_results_clean$Refpop)) + 3
    for(i in seq_len(ilen)){
        output$x$data[[i]]$showlegend  <- FALSE # remove legends for points
        output$x$data[[i]]$legendgroup <- NULL # keeps ellipse legends on
        output$x$data[[i]]$visible     <- TRUE 
    }
    output
}

interactivePlot <- function(obj, metadata=NULL, id=NULL, type=NULL, 
    group=NULL) {

    check_grafpop(obj)
    check_metadata(metadata, id)
    check_variable(group, "group", metadata) 
    check_variable(type, "type", metadata) 

    Refpop <- Sample <- Nearest_neighbor <- Separation_percent <- NULL
    F_percent <- E_percent <- A_percent <- NULL
    label1 <- label2 <- label3 <- GD1_x <- GD2_y <- Refpop_n <- .data <- NULL
    Group <- Type <- NULL

    test_results  <- obj$table
    test_results  <- update_test_results(test_results, metadata, id, 
        group, type)
    train_results <- getTrainResults()

    test_results_clean <- test_results %>%
    mutate(label1 = paste0(Group,", ",Type,", ",Sample),
    label2 = paste0(Refpop,", ", Nearest_neighbor, ", ", 
    round(Separation_percent,0),"%"), label3 = paste0(round(F_percent,0),
    "%, ",round(E_percent,0), "%, ",round(A_percent,0),"%")) %>%
    rename("Group, Type, Sample"=label1,
    "Refpop, Neighbor, Separation"=label2,
    "African, European, Asian Ancestry"=label3)

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

    train_results_clean <- train_results %>%
    left_join(refpop_n,by="Refpop") %>%
    mutate(Refpop_n=factor(Refpop_n,levels=levels(refpop_order$Refpop_n)),
    Refpop=factor(Refpop,levels=unique(refpop_order$Refpop)))

    nb.cols <- 9
    mycolors_fill <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

    # Define the number of colors you want. 9 groups, 9 colors
    nb.cols <- length(unique(test_results[, "Group", drop=TRUE]))
    mycolors <- colorRampPalette(brewer.pal(17, "Set1"))(nb.cols)

    # train is plotted first to get ellipses
    p1 <- ggplot(train_results_clean, aes(GD1_x, GD2_y))+
        # theme stuff
        theme_classic()+
        theme(panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_blank(), legend.position="bottom", 
        axis.text=element_text(face="bold",family="serif", 
        color="black",size=15),
        axis.title = element_text(face="bold",family="serif", 
        color="black",size=15))+
        guides(fill = guide_legend(nrow = 3), color=guide_legend(nrow=3))+
        scale_fill_manual(values = mycolors_fill)+
        scale_color_manual(values = mycolors)+
        scale_x_continuous("GD1",breaks=pretty_breaks(n=10),limits=c(1,1.8))+
        scale_y_continuous("GD2",breaks=pretty_breaks(n=10),limits=c(1,1.4))+
        # ellipse of training data
        stat_ellipse(aes(fill=Refpop_n, label4=Refpop),
        geom = "polygon",alpha=0.25)+
        # make triangle
        geom_segment(x = 1.05, y = 1.1, xend = 1.7658, yend = 1.1)+
        geom_segment(x = 1.05, y = 1.1, xend = 1.4701, yend = 1.2897)+
        geom_segment(x = 1.4701, y = 1.2897, xend = 1.7658, yend = 1.1)+
        # data points
        geom_point(data=test_results_clean, 
        aes(label1=.data[["Group, Type, Sample"]],
        label2=.data[["Refpop, Neighbor, Separation"]],
        label3=.data[["African, European, Asian Ancestry"]],
        color=Group),size=1.25)

    output <- ggplotly(p1, tooltip=c("label1"="Group, Type, Sample",
    "label2"="Refpop, Neighbor, Separation",
    "label3"="African, European, Asian Ancestry",
    "label4"="Refpop")) %>%
    layout(legend = list(orientation = "h",y = -0.15)) %>% 
    config(showTips = FALSE) 

    output$x$layout$legend$title$text <- "" # remove legend title
    output$x$layout$legend$font$size  <- 12 # change legend font size

    for(i in seq_len(length(output$x$data))){
        #fix legend names
        if (!is.null(output$x$data[[i]]$name)){ 
            output$x$data[[i]]$name <- gsub("\\(","",
            str_split(output$x$data[[i]]$name,",")[[1]][1]) 
        }

        # adds spaces for hover variables
        if (!is.null(output$x$data[[i]]$text)){ 
            output$x$data[[i]]$text <- gsub("_"," ",output$x$data[[i]]$text) 
        }
    }

    # for ellipses, change hover info to remove sample size
    ilen <- length(unique(train_results_clean$Refpop)) + 3
    for(i in seq_len(ilen)){
        output$x$data[[i]]$showlegend <- FALSE 
        output$x$data[[i]]$legendgroup <- NULL 
        output$x$data[[i]]$visible  <- TRUE
    }
    output
}

update_test_results <- function(test_results, metadata, id, group, type) {

    # Add type, group
    x    <- setup_metadata(metadata, id, test_results, 2)
    rows <- match(test_results[, 1], x[, 1]) 
    tmp  <- !is.na(rows)
    rows <- rows[tmp]
    n    <- length(rows)
    v    <- "Type"
    test_results[, v] <- "NA"
    if (length(type) && n) {
        test_results[tmp, v] <- x[rows, type, drop=TRUE]
    }
    v2 <- "Group"
    test_results[, v2] <- "NA"
    if (length(group) && n) {
        test_results[tmp, v2] <- x[rows, group, drop=TRUE]
    } else {
        test_results[, v2] <- test_results[, v, drop=TRUE]
    }
    test_results  
}


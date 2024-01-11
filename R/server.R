server <- function(input, output, session) {
  
  observe({
    req(input$navigate_to)
    if (input$navigate_to == "requirements") {
      updateTabItems(session, "tabs", selected = "requirements")
    }
  })
  
  metabolomics_data <- reactive({
    req(input$Data)
    read.csv(input$Data$datapath, check.names = F, row.names = NULL)
  })
  
  metabolomics_meta <- reactive({
    req(input$MetaAnn)
    read.csv(input$MetaAnn$datapath)
  })
  
  samples_meta <- reactive({
    req(input$SampleAnn)
    samples_raw <- read.csv(input$SampleAnn$datapath, na.strings=c("NA","NaN", ""))
    
    # A helper function to determine if a value should be treated as NA
    is_special_na <- function(val) {
      return(val %in% c("", "-", "."))
    }
    
    # Iterate through each column of the matrix
    for(i in 1:ncol(samples_raw)) {
      # Convert special cases to NA
      samples_raw[, i][is_special_na(samples_raw[, i])] <- NA
      
      # If the column has less than 5 unique values and is not already a factor
      if(length(unique(samples_raw[, i])) < 5 && !is.factor(samples_raw[, i])) {
        samples_raw[, i] <- as.factor(samples_raw[, i])
        
        # Only modify invalid names
        levels(samples_raw[, i]) <- sapply(levels(samples_raw[, i]), function(val) {
          if (!identical(make.names(val, unique=TRUE), val) && !is_special_na(val)) {
            return(paste(names(samples_raw)[i], "_", val, sep=""))
          }
          return(val)
        })
      }
    }
    
    return(samples_raw)
  })

  output$raw_meta_data <- renderDataTable({
    req(input$Data)
    DT::datatable(head(metabolomics_data()),
                  options = list(scrollX=TRUE))
  })
  
  output$meta_ann <- renderDataTable({
    req(input$MetaAnn)
    DT::datatable(head(metabolomics_meta()),
                  options = list(scrollX=TRUE))
    
  })
  
  output$sample_ann <- renderDataTable({
    req(input$SampleAnn)
    
    DT::datatable(head(samples_meta()),
                  options = list(scrollX=TRUE))
    
  })
  
  imputed_data <- reactive({
    req(input$Data)
    req(input$threshold)
    
    #metabolites_data <- read.csv(input$Data$datapath, check.names = F, row.names = 1)
    
    # clean data
    int.mat <- metabolomics_data()
    rowNms <- rownames(int.mat)
    colNms <- colnames(int.mat)
    naNms <- sum(is.na(int.mat))
    num.mat <- suppressWarnings(apply(int.mat, 2, as.numeric))
    num.mat <- apply(int.mat, 2, function(x) as.numeric(gsub(",", "", x)))
    int.mat <- num.mat
    
    # remove columns with low variance
    rownames(int.mat) <- rowNms
    colnames(int.mat) <- colNms
    varCol <- suppressWarnings(apply(int.mat, 2, var, na.rm = T))
    constCol <- (varCol == 0 | is.na(varCol))
    constNum <- sum(constCol, na.rm = T)
    if (constNum > 0) {
      int.mat <- int.mat[, !constCol]
    }
    
    # Cleaning stats
    totalCount <- nrow(int.mat) * ncol(int.mat)
    naCount <- sum(is.na(int.mat))
    naPercent <- round(100 * naCount/totalCount, 1)
    
    # Removes columns with more missing the values than the threshold
    limit <- (input$threshold / 100)
    good.inx <- apply(is.na(int.mat), 2, sum) / nrow(int.mat) < limit
    int.mat1 <- as.data.frame(int.mat[,good.inx])
    
    # Imputes the data using k-nearest neighbor imputation
    new.mat2 <- t(impute::impute.knn(t(int.mat1))$data)
    
    return(new.mat2) 
  })
  
  cleaned_data <- reactive({
    req(imputed_data)
    new.mat2 <- imputed_data()
    # Log Transform
    min.val <- min(abs(new.mat2[new.mat2 != 0]))/2
    norm.data <- log2((new.mat2 + sqrt(new.mat2^2 + min.val))/2)
    
    return(norm.data)
  })
  
  # Reactive expression to hold removed columns
  removed_columns <- reactive({
    # Apply threshold logic to the dataset based on input$threshold
    # This is just a placeholder; you'll need to implement the actual logic
    # Example: cleaned_data <- data[data$percentage >= input$threshold, ]
    
    # Find the column names in metabolomics_data that are not in cleaned_data
    cols <- setdiff(colnames(metabolomics_data()), colnames(cleaned_data()))
    
    # Filter out any empty character strings
    cols[cols != ""]
  })
  
  # Info box for removed columns count
  output$removedBox <- renderInfoBox({
    infoBox(
      "Removed Metabolites", length(removed_columns()), icon = icon("minus-circle"),
      color = "red"
    )
  })
  
  # Info box for remaining columns count
  output$remainingBox <- renderInfoBox({
    infoBox(
      "Remaining Metabolites", ncol(cleaned_data()), icon = icon("check-circle"),
      color = "green"
    )
  })
  
  # Plotting the distribution before normalization (using imputed_data)
  output$beforeCleaningPlot <- renderPlot({
    # Flatten the data to a vector
    values_before <- as.numeric(as.matrix(imputed_data()))
    values_before <- values_before[!is.na(values_before)]  # Remove NAs
    
    # Suppress the specific warning and plot
    suppressWarnings({
      ggplot(data.frame(value = values_before), aes(x = value)) +
        geom_density(fill = "skyblue", alpha = 0.7) +
        labs(title = "Distribution Before Normalization", x = "Value", y = "Density") +
        theme_minimal()
    })
  })
  
  # Plotting the distribution after cleaning (using cleaned_data)
  output$afterCleaningPlot <- renderPlot({
    # Flatten the data to a vector
    values_after <- as.numeric(as.matrix(cleaned_data()))
    values_after <- values_after[!is.na(values_after)]  # Remove NAs
    
    # Suppress the specific warning and plot
    suppressWarnings({
      ggplot(data.frame(value = values_after), aes(x = value)) +
        geom_density(fill = "lightgreen", alpha = 0.7) +
        labs(title = "Distribution After Cleaning", x = "Value", y = "Density") +
        theme_minimal()
    })
  })
  
  cleaned_meta <- reactive({
    removed_metabolites <- (setdiff(metabolomics_meta()$BIOCHEMICAL, colnames(cleaned_data())))
    cleaned_meta <- (metabolomics_meta()[!metabolomics_meta()$BIOCHEMICAL %in% (removed_metabolites),] )
    rownames(cleaned_meta)=NULL
    return(cleaned_meta)
  })
  
  observe({
    
    options = colnames(samples_meta())
    
    # Filter the columns based on the number of unique values
    filtered_columns <- samples_meta() %>%
      select(where(function(x) length(unique(x)) < 6)) %>%
      colnames()
    
    updateSelectInput(session, "PCA_shape", choices = filtered_columns, selected = 0)
    updateSelectInput(session, "PCA_color", choices = options, selected = 0)
    updateSelectInput(session, "sample_name_var", choices = options, selected = 0)
    updateSelectInput(session, "factors_var", choices = options, selected = 0)
    updateSelectInput(session, "group_column", choices = options, selected = 0)
  })
  
  observe({
    updateSelectInput(session, "pathway", choices = colnames(metabolomics_meta()), selected = 0)
    #updateSelectInput(session, "super_pathway", choices = colnames(metabolomics_meta()), selected = 0)
  })
  
  PCAplot <- reactive({
    req(input$PCA_shape)
    req(input$PCA_color)
    
    df_pca <- prcomp(cleaned_data())
    df_out <- as.data.frame(df_pca$x)
    PCA_color <- input$PCA_color
    PCA_shape <- input$PCA_shape
    ggplot(df_out,aes(x=PC1,y=PC2,color=samples_meta()[, c(PCA_color)], shape=samples_meta()[, c(PCA_shape)]))+
      geom_point()+ggtitle("")+labs(color='')+
      geom_point(size=8,alpha=0.5)+ #Size and alpha just for fun
      theme(  plot.title = element_text(hjust = 0.5,size=15,face = "bold"),
              axis.text.x = element_text( size = 15, angle = 45, hjust = .5, vjust = 0.5, face = "plain"),
              axis.text.y = element_text( size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
              axis.title.x = element_text( size = 15, angle = 0, hjust = .5, vjust = 0, face = "bold"),
              axis.title.y = element_text( size = 15, angle = 90, hjust = .5, vjust = .5, face = "bold"),
              #legend.title=element_text(size=20),
              legend.title=element_blank(), # remove legend title name
              legend.text = element_text(size=15,face="plain"),
              strip.text = element_text(size = 15,face="plain") ,
              legend.position="right",
              
              # for transparent background
              panel.background = element_rect(fill = "transparent"), # bg of the panel
              plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
              panel.grid.major = element_blank(), # get rid of major grid
              panel.grid.minor = element_blank(), # get rid of minor grid
              legend.background = element_rect(fill = "transparent"), # get rid of legend bg
              legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
              axis.line = element_line(colour = "black") # adding a black line for x and y axis
      ) +xlab(paste0("PC 1 (", round(df_pca$sdev[1],1),"%)")) +
      ylab(paste0("PC 2 (", round(df_pca$sdev[2],1),"%)"))
  })
  
  output$PCA_plot <- renderPlot({
    req(input$Data)
    req(input$SampleAnn)
    req(input$PCA_shape)
    req(input$PCA_color)
    
    print(PCAplot())
  })
  
  #observeEvent (input$submit_factors,{})
  
  int_factors <- reactive({
    req(input$factors_var)
    
    return(input$factors_var)
  })
  
  F_for_all <- reactive({
    req(int_factors)
    
    interested_factors <- int_factors()
    # Add "SAMPLE_NAME" to the beginning of the vector
    interested_factors <- c("SAMPLE_NAME", interested_factors)
    
    samples_meta_interested <- subset(samples_meta(), select=interested_factors)
    
    F_for_all1 <- matrix(0,nrow=ncol(samples_meta_interested),ncol = ncol(cleaned_data()))
    cobre_all <- cbind(samples_meta_interested,cleaned_data())
    
    factors_formula <- paste(interested_factors[!interested_factors == "SAMPLE_NAME"], collapse = "+")
    formula_str <- paste0("cobre_all[,(i+ncol(samples_meta_interested))] ~ ", factors_formula)
    lm_formula <- as.formula(formula_str)
    
    for (i in 1:ncol(cleaned_data())){
      
      lm.out2 <- lm(lm_formula, data=cobre_all, na.action=na.omit) 
      aa=summary(lm.out2)
      F_for_all1[nrow(F_for_all1),i] <- mean(aa$coefficients[,2])
      
      zz=Anova(lm.out2)
      zz=na.omit(zz)
      F_for_all1[1:(nrow(F_for_all1)-1),i] <- zz$`F value`
    }
    
    return(F_for_all1)
  })
  
  SOVplot <- reactive({
    #req(input$sample_name_var)
    req(input$factors_var)
    req(int_factors)
    req(input$submit_factors)
    
    interested_factors <- int_factors()
    # Add "SAMPLE_NAME" to the beginning of the vector
    interested_factors <- c("SAMPLE_NAME", interested_factors)
    
    samples_meta_interested <- subset(samples_meta(), select=interested_factors)
    #str(samples_meta_interested)
    
    F_for_all1 <- F_for_all()
    
    zz1 <- apply(F_for_all1,1,mean)
    zz1 <- as.data.frame(t(zz1))
    dim(zz1)
    colnames(zz1) <- c(c(interested_factors[-1],"error"))
    zz2 <- melt(zz1)
    
    ggplot(data=zz2,aes(x=variable,y=value)) +
      geom_bar(stat="identity", aes(fill=variable)) + # use variable to determine fill color
      labs(title="Sources of variation (Overall)", x="Factors", y = "Mean F Ratio") +
      theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"), # increase title size and make it bold
            axis.text.x = element_text(angle = 45, hjust = 1, size=14), # increase size and make x-axis text bold
            axis.title.x = element_text(size=16, face ="bold"),
            axis.text.y = element_text(size=14), # increase size and make y-axis text bold
            axis.title.y = element_text(size=16, face="bold"), # increase y-axis label size and make it bold
            panel.grid.major = element_blank(), # remove major grid lines
            panel.grid.minor = element_blank()) # remove minor grid lines
    
  })
  
  output$overall_variation <- renderPlot({
    req(SOVplot)
    
    print(SOVplot())
  })
  
  pathway_var_plot <- reactive({
    #req(input$sample_name_var)
    req(input$factors_var)
    req(int_factors)
    req(input$submit_factors)
    req(input$pathway)
    
    interested_factors <- int_factors()
    # Add "SAMPLE_NAME" to the beginning of the vector
    interested_factors <- c("SAMPLE_NAME", interested_factors)
    
    F_for_all1 <- F_for_all()
    
    rownames(F_for_all1) <- c(c(interested_factors[-1],"error"))
    colnames(F_for_all1) <- colnames(cleaned_data())
    
    ss <- t(F_for_all1)
    
    sub_input <- input$pathway
    
    ss_ALL1 <- data.frame(ss, pathway=cleaned_meta()[[sub_input]])
    ss_ALL1 %>% group_by(pathway) %>% mutate(N_metaboliesinP=(n())) %>%summarise_all(funs(mean)) %>% 
      mutate(pathway=paste(pathway,N_metaboliesinP,sep="  ")) %>% dplyr::select(-N_metaboliesinP,-error) %>%as.data.frame()%>% melt()%>% 
      arrange(pathway,desc(value)) %>%
      ggplot(aes(x=pathway,y=value,fill=variable))+
      geom_bar(stat="identity")+
      labs(title="Sources of variation (Pathway Specific) ",x="", y = "Mean F Ratio")+
      theme(plot.title = element_text(hjust = 0.5))+ theme(axis.text.x = element_text(size=7,angle = 45, hjust =1)) +
      scale_fill_manual(values=c( "#00ba38", "#00bfc4","#619cff", "#b79f00"))
  })
  
  output$pathway_var <- renderPlot({
    req(pathway_var_plot)
    
    print(pathway_var_plot())
  })
  
  observe({
    req(input$group_column)
    unique_values <- unique(samples_meta()[[input$group_column]])
    updateSelectInput(session, "first_outcome", choices = unique_values, selected = 0)
    updateSelectInput(session, "second_outcome", choices = unique_values, selected = 0)
  })
  
  output_limma_table <- reactive({
    req(input$p.adjust)
    req(input$group_column)
    req(input$first_outcome)
    req(input$second_outcome)
    
    # Get the user-selected column and convert it to a factor
    type <- factor(as.character(samples_meta()[[input$group_column]]))

    # Create a design matrix based on the selected outcome
    design <- model.matrix(~0 + type)

    colnames(design) = levels(type)

    A <<- input$first_outcome
    B <<- input$second_outcome

    contrast <- makeContrasts(paste(A, B, sep="-"), levels=unique(samples_meta()[[input$group_column]]))

    # Fit the linear model
    fit <- limma::lmFit(as.matrix(t(cleaned_data())), design)

    fit <- limma::contrasts.fit(fit, contrast)

    # Compute empirical Bayes statistics
    fit <- limma::eBayes(fit)

    # Get the top table of results
    result <- limma::topTable(fit, adjust.method=input$p.adjust, number=Inf, p.value = 1, coef = 1)
  })
  
  output$limma_table <- renderDataTable({
    req(input$p.adjust)
    #req(input$p.value)
    req(input$group_column)
    req(input$first_outcome)
    req(input$second_outcome)
    req(output_limma_table)
    
    DT::datatable(head(output_limma_table()),
                  options = list(scrollX=TRUE))
  })
  
  volcano_plot <- reactive({
    # Add a column that describes up- or down-regulation
    output_limma_table <- output_limma_table() %>% 
      mutate(Regulation = case_when(
        logFC > 0 & -log10(adj.P.Val) > 1.3 ~ "Up",
        logFC < 0 & -log10(adj.P.Val) > 1.3 ~ "Down",
        TRUE ~ "Not significant"
      ))
    
    # Define color palette
    color_palette <- c("Up" = "red", "Down" = "blue", "Not significant" = "grey")
    
    # Plot
    p <- ggplot(output_limma_table, aes(x = logFC, y = -log10(adj.P.Val), color = Regulation)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = color_palette) +
      geom_hline(yintercept = 1.3, linetype = "dashed") +
      labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value", color = "Regulation") +
      theme_bw() +
      theme(
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    # Add labels for top 5 most significant points in each category
    top_up_genes <- output_limma_table %>% 
      filter(Regulation == "Up") %>% 
      arrange(adj.P.Val) %>% 
      head(5)
    
    top_down_genes <- output_limma_table %>% 
      filter(Regulation == "Down") %>% 
      arrange(adj.P.Val) %>% 
      head(5)
    
    p + 
      geom_text_repel(data = top_up_genes,
                      aes(label = rownames(top_up_genes)),
                      box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.5, "lines"),
                      segment.color = 'grey50',
                      segment.size = 0.2) +
      geom_text_repel(data = top_down_genes,
                      aes(label = rownames(top_down_genes)),
                      box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.5, "lines"),
                      segment.color = 'grey50',
                      segment.size = 0.2)
    
  })
  
  output$limma_volcano_plot <- renderPlot({
    req(volcano_plot)
    
    print(volcano_plot())
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      "PCA.png"
    },
    
    content = function(file) {
      ggsave(file, plot = PCAplot(), bg = "white", device = "png")
    }
  )
  
  output$downloadSOV <- downloadHandler(
    filename = function() {
      "SOV.png"
    },
    
    content = function(file) {
      ggsave(file, plot = SOVplot(), device ="png")
    }
  )
  
  output$downloadPathwayVar <- downloadHandler(
    filename = function() {
      "PathwayVar.png"
    },
    
    content = function(file) {
      req(input$plot_width, input$plot_height, input$plot_dpi)
      
      # Convert width, height, and dpi to numeric
      plot_width <- as.numeric(input$plot_width)
      plot_height <- as.numeric(input$plot_height)
      plot_dpi <- as.numeric(input$plot_dpi)
      
      ggsave(file, plot = pathway_var_plot(), width = plot_width, height = plot_height, dpi = plot_dpi, device ="png")
    }
  )
  
  output$download_limma <- downloadHandler(
    filename = function() {
      paste("limma_table-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv((output_limma_table()), file)
    }
  )
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("volcano_plot-", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      ggsave(file, plot = volcano_plot())
    }
  )
  
  # Download handler for removed columns
  output$download_removed <- downloadHandler(
    filename = function() {
      paste("removed_metabolites_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data.frame(`Removed Metabolites` = removed_columns()), file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  # Download handler for normalized (cleaned) data
  output$download_normalized <- downloadHandler(
    filename = function() {
      paste("normalized_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Using 'cleaned_data' as the source for normalized data
      write.csv(cleaned_data(), file, row.names = TRUE)
    },
    contentType = "text/csv"
  )
  
  output$download_meta_data <- downloadHandler(
    filename = function() {
      "example_metabolomics_data.csv"
    },
    content = function(file) {
      file.copy("www/example_metabolomics_data.csv", file)
    }
  )
  
  output$download_meta_ann <- downloadHandler(
    filename = function() {
      "example_metabolomics_annotation.csv"
    },
    content = function(file) {
      file.copy("www/example_metabolomics_annotation.csv", file)
    }
  )
  
  output$download_sample_ann <- downloadHandler(
    filename = function() {
      "example_sample_annotation.csv"
    },
    content = function(file) {
      file.copy("www/example_sample_annotation.csv", file)
    }
  )
}

library(shiny)
library(shinydashboard)
library(shinyalert)
library(shinycssloaders)
library(shinyWidgets)

if (!require("DT")) install.packages('DT')
library(DT)

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(tidyverse)

library(car) 
library(dplyr)
library(reshape)

library(limma)
library(tibble)
library(EnhancedVolcano)

library(rsconnect)

options(shiny.maxRequestSize = 150 * 1024 * 1024) # 150 MB in bytes

ui <- dashboardPage(
  dashboardHeader(title = "MetabolonR"),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("Getting Started", tabName = "dashboard", icon = icon("dashboard")),
      #menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
      menuItem("Preview Data", tabName = "preview", icon = icon("search")),
      menuItem("Clean Data", tabName = "filter", icon = icon("filter")),
      #menuItem("Clean Data", tabName = "clean", icon = icon("filter")),
      menuItem("PCA", tabName = "PCA", icon = icon("chart-bar")),
      menuItem("Souces of Variation", tabName = "Variation", icon = icon("greater-than-equal")),
      menuItem("Biomarker Discovery", tabName = "Limma", icon = icon("bolt")),
      menuItem("Data Formatting Requirements", tabName = "requirements", icon = icon("file-csv"))
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content - home page
      tabItem(tabName = "dashboard",
              fluidPage(
                column(12,
                       align="center",
                       div(img(src="logo.png", align = "center"),
                           style="display: inline-block;")),
                div(src="logo.png", style="text-align: center;"),
                column(12,
                       p(strong("MetabolonR"), "is an interactive shiny App that provides an analysis pipeline for metabolomics data"
                       ), style = "font-family: 'times'; font-si20pt",

                       br(),
                       tags$ul(
                         tags$li("Users can upload metabolomic data, samples annotation, and metabolomics annotation.", 
                                 a(href="#", onclick="Shiny.setInputValue('navigate_to', 'requirements', {priority: 'event'});", "Click here"), 
                                 " for data formatting requirements."),
                         tags$li("Users can apply a variance filter, remove high missing metabolites,
                                 impute metabolites with small missingness rate, and log transformation."),
                         tags$li("Using the information from the sample annotation file, users can
                                 identify how much of variance is explained by age, sex or other covariant factors."),
                         tags$li("Users can select two groups to identify highly significant
                                 metabolites between these two groups.")
                       )
                ),
                headerPanel(""),
                fluidRow(
                  column(5,
                         wellPanel(
                           h3("Upload Data Here"),
                           hr(),
                           div(
                             fileInput("Data","Upload Metabolomics Data", multiple = TRUE,
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv"))
                           ),
                           div(
                             fileInput("MetaAnn","Upload Metabolomics Annotations", multiple = TRUE,
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv"))
                           ),
                           div(
                             fileInput("SampleAnn","Upload Sample Annotations", multiple = TRUE,
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv"))
                           )
                         )
                  )
                )
              )
      ),
      # Second tab content - upload data
      tabItem(tabName = "upload",
              
      ),
      # Third tab content - select gene
      tabItem(tabName = "preview",
              fluidRow(
                box(width = "700", title = "Metabolomics Data", withSpinner(DT::dataTableOutput("raw_meta_data")))
              ),
              fluidRow(
                box(width = "700", title = "Metabolomics Annotations", withSpinner(DT::dataTableOutput("meta_ann")))
              ),
              fluidRow(
                box(width = "700", title = "Samples Annotations", withSpinner(DT::dataTableOutput("sample_ann")))
              )
      ),
      tabItem(tabName = "filter",
              fluidPage(
                
                # Sidebar layout with input and output definitions
                sidebarLayout(
                  sidebarPanel(
                    # Slider input for percent threshold with default value
                    sliderInput("threshold", 
                                "Select Missingness Rate:", 
                                min = 0, 
                                max = 100, 
                                value = 50, 
                                step = 1, 
                                post = "%"),
                    
                    # Add a line break for better separation
                    tags$br(),
                    
                    # Adding instructional text
                    wellPanel(
                      h4("Instructions:"),
                      tags$ul(
                        tags$li("Use the slider above to set a missingness rate"),
                        tags$li("Metabolites with a missingness rate above this threshold will be removed from the analysis."),
                        tags$li("Ex: The default is set to 50%, meaning metabolites with >50% missing values are deleted.")
                      )
                    ),
                    
                    # Details on data cleaning
                    wellPanel(
                      h4("Data Cleaning Details:"),
                      tags$ul(
                        tags$li("Metabolites with variance = 0 will be removed."),
                        tags$li("Metabolites with missingness rate higher than the selected value will be removed from the downstream analysis."),
                        tags$li("The remaining missing metabolite values will be imputed using the KNN algorithm implemented in the impute R package."),
                        tags$li("Data is then log transformed. Log transformation stabilizes variance, makes data more normal, and reduces the impact of outliers.")
                      )
                    ),
                    
                    # Button to download the list of removed columns with some styling and spacing
                    tags$div(
                      downloadButton("download_removed", "Download Removed Metabolites"),
                      style = "margin-bottom: 10px; display: block; width: 100%;"
                    ),
                    
                    # Button to download the normalized data with full width
                    tags$div(
                      downloadButton("download_normalized", "Download Normalized Data"),
                      style = "display: block; width: 100%;"
                    )
                  ),
                  
                  mainPanel(
                    # Display info boxes
                    fluidRow(
                      infoBoxOutput("removedBox"),
                      infoBoxOutput("remainingBox")
                    ),
                    
                    # Density plots for data distribution before and after cleaning
                    fluidRow(
                      column(6, plotOutput("beforeCleaningPlot", height = "400px")),
                      column(6, plotOutput("afterCleaningPlot", height = "400px"))
                    )
                  )
                )
              )
            ),
      tabItem(tabName = "clean",
              fluidPage(
                verticalLayout(
                  column(12, align = "center", style = "font:bold",
                         tags$h2("Screening and Filtration Report:"),
                         tags$h5("(We removed the metabolites with low variability and high missingness, then used KNN for imputation and log2 for transformation.)"),
                  ),
                  column(12,
                         box(width = "700", title = "Ensured All Metabolite Values Are Numeric", status = "info",
                             span(htmlOutput("naMessage"), style = "font-size: 25px"),
                             tags$h3("2. Calculated the percentage of missinginess:"),
                             infoBoxOutput("percentage")
                         )),
                  column(12,
                         box(width = "700", title = "Metabolite Data Filtration", status = "info", 
                             tags$h3("3. Removed metabolites with a constant or single value across metabolite samples."),
                             span(htmlOutput("deleted"), style = "font-size: 20px"),
                             tags$h3("4. Removed metabolites if they were missed across more than selected threshold of the sample."),
                             sliderInput("threshold", "Percent Threshold:",
                                         min = 0, max = 100,
                                         value = 50),
                             span(htmlOutput("fifdeleted"), style = "font-size: 20px"),
                             # tags$h3("-----------------------"),
                             # infoBoxOutput("Totalrem")
                             # uiOutput("RemoveFilter"),
                             uiOutput("RemoveFilter"),
                             tags$h3('For more detail on what was removed, based on pathway, navigate to'),
                             #actionBttn("Clicktab4", label = "What Was Removed?",style="pill", color = "success", size="sm")
                         )),
                  
                  column(12,
                         box(width = "700", title = "Imputation and Log Transformion", status = "info",
                             tags$h3("5. Imputed the remaining missed metabolites using K-Nearest Neighbours."),
                             tags$h3("6. Log transformed using log2 function."),
                             tags$h3("----------------------------------------------------------------------------------------------------------------------------------------"),
                             tags$h4("Imputed and Log Transformed Data:", style = "font:bold"),
                             withSpinner(DT::dataTableOutput("IMPUTED"))
                         )),
                  
                  column(12,
                         box(width = "700", title = "Filtration Complete", status = "info",
                             tags$h3("7. You may download modified data or filtration summary. Preview data distrubution on next page."),
                             downloadButton("downloadData", "Download Data"),
                             downloadButton("report", "Download Report"))
                  )
                ),
                # column(2,offset = 10,
                #        actionBttn("Click6", label = "Data Distribution",style="minimal", color = "success", size="lg"))
              )
      ),
      tabItem(tabName = "PCA",
              fluidRow(
                column(
                  3,
                  wellPanel(
                    h3("Settings"),
                    hr(),
                    div(
                      selectInput(inputId = "PCA_shape",
                                  label = "Select PCA Shape",
                                  choices = NULL)
                    ),
                    div(
                      selectInput(inputId = "PCA_color",
                                  label = "Select PCA Color",
                                  choices = NULL)
                    ),
                  ),
                ),
                column(
                  7,
                  wellPanel(
                    h3("PCA"),
                    hr(),
                    plotOutput("PCA_plot"), downloadButton("downloadPlot", "Download Plot")
                  )
                )
                
              )
      ),
      tabItem(
        tabName = "Variation",
        fluidRow(
          column(
            3,
            wellPanel(
              h3("Settings"),
              hr(),
              # div(
              #   selectInput(inputId = "sample_name_var",
              #               label = "Sample_Name column:",
              #               choices = NULL)
              # ),
              div(
                selectInput(inputId = "factors_var",
                            label = "Interested Factors:",
                            multiple = TRUE,
                            choices = NULL)
              ),
              div(
                actionBttn("submit_factors", label = "Submit", style="pill", color = "success", size="sm")
              )
            )
          ),
          column(
            7,
            wellPanel(
              h3("Overall Variation"),
              hr(),
              plotOutput("overall_variation"),  downloadButton("downloadSOV", "Download Plot")
            )
          )
        ),
        fluidRow(
          column(3,
                 wellPanel(
                   tags$h3("Sub/Super Pathway Selection"),
                   hr(),
                   div(
                     selectInput(inputId = "pathway",
                                 label = "Select sub/super pathway column:",
                                 choices = NULL)
                   )
                 ),
                 wellPanel(
                   tags$h4("Download Plot Settings"),
                   hr(),
                   selectInput("plot_width", "Select plot width:", choices = seq(5, 50, by = 1), selected = 10),
                   selectInput("plot_height", "Select plot height:", choices = seq(5, 50, by = 1), selected = 10),
                   selectInput("plot_dpi", "Select plot resolution (dpi):", choices = c(72, 96, 150, 300), selected = 300)
                 )
          ),
          column(9,
                 wellPanel(
                   tags$h3("Source of Variation per sub/super pathway"),
                   hr(),
                   splitLayout(
                     plotOutput("pathway_var", height = "300px", width = 1500)
                     # plotOutput("super_pathway_var", height = "300px", width = "100%")
                   ),
                   downloadButton("downloadPathwayVar", "Download Plot")
                 )
          )
        )
      ),
      tabItem(tabName = "Limma",
              fluidRow(
                column(2,
                       wellPanel(
                         h3("Settings"),
                         hr(),
                         div(
                           tags$label("Select p.adjust method:"),
                           selectInput("p.adjust", NULL, 
                                       choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
                         ),
                         # div(
                         #   tags$label("Enter p-value:"),
                         #   numericInput("p.value", NULL, value = 0.05)
                         # ),
                         div(
                           tags$label("Select group/outcome column:"),
                           selectInput("group_column", NULL, NULL)
                         ),
                         div(
                           tags$label("Select first group:"),
                           selectInput("first_outcome", NULL, NULL)
                         ),
                         div(
                           tags$label("Select second group:"),
                           selectInput("second_outcome", NULL, NULL)
                         )
                       )
                ),
                column(10,
                       wellPanel(
                         h3("Analysis Results"),
                         hr(),
                         withSpinner(DT::dataTableOutput("limma_table")), downloadButton("download_limma", "Download Limma Table")
                       ),
                       column(8,
                              wellPanel(
                                h3("Volcano Plot"),
                                hr(),
                                plotOutput("limma_volcano_plot"), downloadButton("download_plot", "Download Plot")
                              )
                       )
                )
              )
      ),
      tabItem(tabName = "requirements",
              fluidPage(
                box(title = "Metabolomics Data", status = "primary", solidHeader = TRUE, 
                    p("1. ", strong("File Format:"), " The data should be in CSV format."),
                    p("2. ", strong("Row Names:"), " Each row should represent a unique sample. The ", em("rownames"), " (the entry at the start of each row without a corresponding column header) should contain the sample identifiers."),
                    p("3. ", strong("Columns:"), " Each column should represent a distinct gene. Column headers should indicate the gene names or IDs."),
                    p("4. ", strong("Values:"), " The values within the table should represent the gene expression levels or counts for each sample."),
                ),
                
                box(title = "Metabolomics Annotations", status = "primary", solidHeader = TRUE,
                    p("1. ", strong("File Format:"), " The data should be in CSV format."),
                    p("2. ", strong("BIOCHEMICAL Column:"), " There must be a column named '", code("BIOCHEMICAL"), "' (with exact spelling and capitalization). The values in this column should correspond to the column names from the gene data."),
                ),
                
                box(title = "Samples Annotations", status = "primary", solidHeader = TRUE,
                    p("1. ", strong("File Format:"), " The data should be in CSV format."),
                    p("2. ", strong("SAMPLE_NAME Column:"), " There must be a column named '", code("SAMPLE_NAME"), "' (with exact spelling and capitalization). The values in this column should correspond to the row names from the gene data."),
                ),
                
                box(title = "Example Datasets", status = "info", solidHeader = TRUE,
                    p("To assist with formatting, you can download example datasets below:"),
                    div(style = "display: inline-block; padding-right: 10px;",
                        downloadButton(outputId = "download_meta_data", label = "Metabolomics Data", icon = icon("file-download"))
                    ),
                    div(style = "display: inline-block; padding-right: 10px;",
                        downloadButton(outputId = "download_meta_ann", label = "Metabolomics Annotations", icon = icon("file-download"))
                    ),
                    div(style = "display: inline-block;",
                        downloadButton(outputId = "download_sample_ann", label = "Sample Annotations", icon = icon("file-download"))
                    )
                )
              )
      )
    )
  )
)

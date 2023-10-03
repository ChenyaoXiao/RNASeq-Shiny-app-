# Load packages ----------------------------------------------------------------
library(shiny)
library(readr)
library(here)
library(tibble)
library(ggplot2)
library(DT)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(htmlTable)

# Load helper functions and dataset type meta-----------------------------------
source('helpers.R')
dataset_type_meta <- read_excel("data/dataset_type_meta.xlsx")

# Define UI --------------------------------------------------------------------
ui <- fluidPage(
  # App title ----
  titlePanel("Search Gene Symbols in a Selected Dataset"),
  # Make error message red
  tags$head(tags$style(".shiny-output-error{color: red;}")),
  # Sidebar display
  sidebarLayout(
    # Inputs: 
    sidebarPanel(
      # Select program
      radioButtons(
        inputId = "program",
        label = "Select program:",
        choices = c("Program1","Program2"),
        inline = TRUE
      ),
      # Select model type
      radioButtons(
        inputId = "model_type",
        label = "Select model type:",
        choices = c("Cell Line","CDX","Syngeneic and Normal mouse",
                    "PDX","Normal mouse")
      ),
      # Select model name
      selectInput(
        inputId = "model_name",
        label = "Select model name:",
        choices = NULL
      ),
      # Contains normal issue or not
      checkboxInput("normal","Datasets contain normal tissue",value= FALSE),
      # Initialize dataset select
      selectInput(inputId ="dataset",
                  label = "Choose a dataset:",
                  choices = NULL),
      # Enter variable for look up
      textInput("gene","Add gene symbol, 
                separate each by semicolon. (case insensitive)"),
      # Select plot output
      selectInput(
        inputId = "Plot",
        label = "Select plot type",
        choices = c("bar chart", "box plot", "violin plot"),
        selected = "bar chart"
      ),
      # Select plot scale
      radioButtons(
        inputId = "scale",
        label = "Select TPM scale:",
        choices = c("TPM", "log2(TPM+1)")
      ),
      # Select fold change or TPM 
      radioButtons(
        inputId="table_choice",
        label="Choose table display:",
        choices = c("TPM","logFC")
      ),
      # This one is linked by the id 'download'
      downloadButton('downloadData',"Download the data"),
    width = 3),
    # Output
    mainPanel(
      tabsetPanel(
        tabPanel("Plot and table",plotOutput(outputId = "selected_graph"),
                 dataTableOutput(outputId = "tpmtable"), htmlOutput(outputId = "Helping")),
        tabPanel("Browse datasets", dataTableOutput("help_table"))
       ),
      width= 9)
    
  )
)

# Define server ----------------------------------------------------------------
server <- function(input, output, session) {
  # Show meta table for browsing all dataset information 
  output$help_table <- renderDataTable({
      DT::datatable(data = dataset_type_meta,
                    options = list(pageLength = 20),
                    rownames = FALSE)
    })
  # Update different select ui based on program and model
  observeEvent(c(input$model_type,input$program),{
    # Update Model name by program and model type
    choices1 <- subset(dataset_type_meta, Program == input$program & 
                             `Model Type`== input$model_type)$`Model Name`
    # Deal with model name that has multiple values
    choices2 <- sort(unique(unlist(str_split(choices1,", "))))
    if(length(choices1) == length(choices2))
    {updateSelectInput(inputId = "model_name", choices = choices1)}
    else {updateSelectInput(inputId = "model_name", choices = choices2)}
  })
  # For cell line model type with multiple values
  # Need convert to origin when choosing dataset
  check_name <- eventReactive(input$model_name,{
    dataset_type_meta$`Model Name`[
      (grepl(input$model_name,dataset_type_meta$`Model Name`))]
  })
  
  observeEvent(c(input$model_name,input$model_type,input$program),{
    # Update dataset selection based on program, model type and model name
    choices <- sort(subset(dataset_type_meta, Program == input$program & 
                        `Model Type`== input$model_type &
                        `Model Name`%in% check_name() )$Dataset)
    updateSelectInput(inputId = "dataset", choices = choices) 
  })
  # Update select for contains normal issue
  observeEvent(input$normal,{
    if(input$normal==TRUE){
    choices <- sort(subset(dataset_type_meta, Program == input$program & 
                        `Model Type`== input$model_type &
                        `Model Name`%in% check_name()
                      & `Normal Tissue`==TRUE)$Dataset)
    updateSelectInput(inputId = "dataset", choices = choices)}
    else {
      choices <- sort(subset(dataset_type_meta, Program == input$program & 
                          `Model Type`== input$model_type & 
                          `Model Name`%in% check_name())$Dataset)
    updateSelectInput(inputId = "dataset", choices = choices)
    }
  })
  # Multiple symbol selection
  gene_choices <-  reactive ({
    req(input$gene)
    unlist(strsplit(str_replace_all(input$gene, fixed(" "), ""),split = ";"))
  
  })
  # Transfer rds file to TPM data frame
  rna_rds <- reactive ({
    req(input$dataset)
    filter_TPM(str_c("data/",input$dataset,".RDS",sep=""))
  })
  # Create a subset of data filtering for chosen symbol
  rna_subset1 <- reactive({
    req(input$gene)
    # Ignore case in input gene
    # Use exact match
    # dplyr::filter(rna_rds(), grepl(paste0("\\b",gene_choices(),"\\b", collapse="|"), 
                                   # Symbol, ignore.case = TRUE))
    rna_subset <- dplyr::filter(rna_rds(), tolower(Symbol) %in% tolower(gene_choices()))
    # Check tgroup message: If the group name is "GX" or "GrX",
    # then use paste(treat, time, dose, sep='_')
    dplyr::mutate(rna_subset,tgroup = case_when(
      str_detect(tgroup,"^(G|Gr)\\d")==TRUE ~ paste(treat, time, dose, sep='_'),
                  TRUE ~ tgroup))
  })
  # Use partial match
  rna_subset2 <- reactive({
    req(input$gene)
    # Ignore case in input gene
    dplyr::filter(rna_rds(), grepl(paste0(gene_choices(), collapse="|"), 
                                   Symbol, ignore.case = TRUE))
  })
  # Transfer rds file to logFC data frame 
  logFC <- reactive ({
    req(input$dataset)
    filter_FC(str_c("data/",input$dataset,".RDS",sep=""))
  })
  FC_filter <- reactive({
    req(input$gene)
    # Ignore case in input gene
    dplyr::filter(logFC(), tolower(Symbol) %in% tolower(gene_choices()))
  })
  # Make choice for y laxis in figure and then plot
  # Choice 1: box plot
  plot1 <- reactive({
    ggplot(rna_subset1(), aes(x = tgroup, y = get(input$scale), 
                             fill=Symbol)) +
      geom_boxplot(position=position_dodge()) +
      theme_bw()+
      theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
      labs(x = "Treatment Group",y=input$scale)
  })
  # Choice 2: violin plot
  plot2 <- reactive({
   ggplot(rna_subset1(), mapping = aes(
     x = tgroup, y = get(input$scale), fill=Symbol))+
     geom_violin(trim=FALSE,bw=10,position=position_dodge(1))+
     theme_bw()+
     theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
     labs(x = "Treatment Group",y=input$scale)
 })
  # Define standard error of mean function
  std.error <- function(x) sd(x)/sqrt(length(x))
  # Get data summary for bar charts
  data_summary <-reactive({rna_subset1() %>%
   group_by(tgroup,Symbol)%>%
   summarize(mean=mean(get(input$scale)),se=std.error(get(input$scale)))})
  # Choice 3: Bar chart with error bar showing standard error
  plot3 <- reactive({
  ggplot(data_summary(), aes(x=tgroup, y=mean,fill=Symbol)) + 
    geom_bar(stat="identity", color="black", 
            position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                 position=position_dodge(.9))+
    theme_bw()+
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
    labs(x = "Treatment Group",y=input$scale)
  })
  # Return the requested graph
  graphInput <- reactive({
   switch(input$Plot,
          "box plot" = plot1(),
          "violin plot"= plot2(),
          "bar chart" = plot3()
   )
 })
  # Output the plot
  # Error display when input nonexistent genes
  output$selected_graph <- renderPlot({ 
    validate(
      need(all(tolower(gene_choices()) %in% tolower(rna_subset1()$Symbol)),
           paste0("Warning: Can't find symbol:",
                  gene_choices()[!tolower(gene_choices()) %in% tolower(rna_subset1()$Symbol)],
                  ". Either this gene is not expressed or this symbol 
                  is not a HGNC official symbol. Please double check your input."))
    )
    graphInput()
  })
  # Print data table for different choice 
  output$tpmtable <- renderDataTable({
    # If that gene is not found, the table could show any gene name containing the query.
    if(input$table_choice=="TPM" & nrow(rna_subset1())== 0){
      DT::datatable(data = rna_subset2()%>% 
                      dplyr::select(Symbol,ENTREZID,sample,TPM,'log2(TPM+1)',tgroup,
                                    treat,time,dose)%>%
                      rename("Treatment Group"=tgroup,Sample=sample,Treatment=treat,
                             Time=time,Dose=dose)%>%
                      mutate(TPM= round(TPM,2)),
                    options = list(pageLength = 10),
                    rownames = FALSE)
    }
    # if a gene is found, the table should only show that gene
    else if(input$table_choice=="TPM"){
      DT::datatable(data = rna_subset1()%>% 
                      dplyr::select(Symbol,ENTREZID,sample,TPM,'log2(TPM+1)',
                                    tgroup,treat,time,dose)%>%
                      rename("Treatment Group"=tgroup,Sample=sample,Treatment=treat,
                             Time=time,Dose=dose)%>%
                      mutate(TPM= round(TPM,2)),
                    options = list(pageLength = 10),
                    rownames = FALSE)
    }
    else{
      DT::datatable(data = FC_filter()%>%
                      mutate_if(is.double, round, 4),
                    options = list(scrollX = TRUE))
    }
  })
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      gene_set <- paste(gene_choices(),collapse = "_")
      str_c(input$dataset,"_gene_",gene_set, ".csv", sep="")
    },
    content = function(file) {
      if(input$table_choice=="TPM" & nrow(rna_subset1())== 0){
        write.csv(rna_subset2(), file, row.names = FALSE) 
      }
      else if(input$table_choice=="TPM"){
      write.csv(rna_subset1(), file, row.names = FALSE)}
      else{
      write.csv(FC_filter(), file, row.names = FALSE)}
    }
  )
  output$Helping <- renderPrint({
    if(input$table_choice=="logFC"){
      HTML(paste0("<b>","Message: If you see a blank column in logFC table, 
      that indicates it has a non-significant fold change.","</b>"))}
  })
 }


# Create a Shiny app object ----------------------------------------------------
shinyApp(ui=ui,server=server)

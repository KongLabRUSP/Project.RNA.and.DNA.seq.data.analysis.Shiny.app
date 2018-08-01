# |----------------------------------------------------------------------------------|
# | Project: RNA/DNA-seq pipeline                                                    |
# | Script:   Shiny app to streamline RNA and DNA seq analysis                       |
# | Authors:  Marissa (Meinizi) Zheng, Davit Sargsyan                                |   
# | Created:  07/22/2018                                                             |
# | Modified: 7/27/2018: Added folder path (using **shinyFiles** package), and       |
# |                      ability to run batch file by clicking **Run Program** button|                                                             |
# |----------------------------------------------------------------------------------|
options(stringsAsFactors = FALSE,
        shiny.maxRequestSize = 2^100)

require(shiny)
require(data.table)
require(DT)
library(shinydashboard)
library(shinythemes)
require(shinyFiles)

# Determine the OS
sysOS <- Sys.info()[['sysname']]

# Specify root folder----
if (sysOS == "Linux") {
  volumes <- c("Home" = "/home/administrator",
               "Data Storage" = "/datastorage")
} else if (sysOS == "Windows") {
  volumes <- c("Local Disk (C:)" = "C:")
}

ui <- dashboardPage(dashboardHeader(title = "Upstreaming Work",
                                    dropdownMenu(type = "notifications",
                                                 notificationItem(text = "Task complete at 86%",      # show here the % completion of the whole task, or with a pump up window
                                                                  icon = icon("exclamation-triangle"),
                                                                  status = "warning"))),
                    dashboardSidebar(sidebarMenu(menuItem(text = "Dashboard", 
                                                          tabName = "dashboard", 
                                                          icon = icon("dashboard")),
                                                 menuItem(text = "RNA-seq", 
                                                          tabName = "rna", 
                                                          icon = icon("th")),
                                                 menuItem(text = "DNA-seq", 
                                                          tabName = "dna", 
                                                          icon = icon("th")))),
                    dashboardBody(tabItems(tabItem(tabName = "dashboard",
                                                   h2("Workflow Overview & Instructions")),
                                           tabItem(tabName = "rna",
                                                   sidebarPanel(shinyDirButton(id = "directory", 
                                                                               label = "Folder select", 
                                                                               title = "Please select a folder"),
                                                                # fileInput(inputId = "filesIn",
                                                                #           label = "Import FastQ File(s)",
                                                                #           multiple = TRUE),
                                                                radioButtons(inputId = "isSingleEnd",
                                                                             label = "Select Single or Pair End",
                                                                             choices = list("Single End" = "true", 
                                                                                            "Pair End" = "false"),
                                                                             selected = "true"),
                                                                actionButton(inputId = "go",
                                                                             label = "Run Program")),
                                                   mainPanel(tags$h3("Working Directory"),
                                                             verbatimTextOutput(outputId = "wd"),
                                                             # tags$h3("Uploaded File"),
                                                             # DT:: dataTableOutput(outputId = "fileNames"),
                                                             tags$h3("File Directory"),
                                                             verbatimTextOutput(outputId = "fileDir"),
                                                             tags$h3("Program"),
                                                             verbatimTextOutput("directorypath"))),
                                           tabItem(tabNamvolumes = "dna",
                                                   sidebarPanel( ),
                                                   mainPanel()))))

server <- function(input, output, session) {
  # # Method 1: select files
  # # NOTE: this will force all files to be copied to a temp folder
  # #       and can be very slow. DO NOT USE for now.
  # output$fileNames <- DT::renderDT({
  #   validate(need(input$filesIn != "", ""))
  #   ne <- new.env()
  #   fnames <- data.table(input$filesIn$datapath)
  #   DT::datatable(fnames,
  #                 options = list(pageLength = 10),
  #                 selection = list(mode = "multiple"),
  #                 rownames = FALSE)
  # })
  
  # Method 2: specify folder containing FastQ file.
  # NOTE: all the files in the folder will be passed to the pipeline
  shinyDirChoose(input = input, 
                 id = "directory", 
                 roots = volumes, 
                 session = session, 
                 restrictions = system.file(package = "base"))
  
  output$wd <- renderPrint({
    getwd()
  })
  
  output$fileDir <- renderPrint({
      parseDirPath(roots = volumes, 
                   selection = input$directory)
  })
  
  output$directorypath <- renderPrint({
    wd <- getwd()
    # myString <-  parseDirPath(roots = volumes, 
    #                           selection = input$directory)
    myString <- paste('\"nextflow run rna_seq.nf --SingleEnd=\"',
                      input$isSingleEnd,
                      '\" --reads=\"',
                      parseDirPath(roots = volumes, 
                                   selection = input$directory),
                      '/*.fastq.gz\" -w \"/scratch/${USER}/NFWorkDir/${PWD}/work\" -resume\"',
                      sep = "")
    pgm <- paste(wd,
                 "/test_run.bat ",
                 myString,
                 " > ",
                 wd,
                 "/log.txt",
                 sep = "")
    
    # Print to Shiny textbox
    cat(myString)
  })
  
  observeEvent(input$go, {
    wd <- getwd()
    myString <- paste('\"nextflow run rna_seq.nf --SingleEnd=\"',
                      input$isSingleEnd,
                      '\" --reads=\"',
                      parseDirPath(roots = volumes, 
                                   selection = input$directory),
                      '/*.fastq.gz\" -w \"/scratch/${USER}/NFWorkDir/${PWD}/work\" -resume\"',
                      sep = "")
    pgm <- paste(wd,
                 "/test_run.bat ",
                 myString,
                 " > ",
                 wd,
                 "/log.txt",
                 sep = "")
    # Run the program - batch file printing myString to file
    system(pgm)
  })
}

shinyApp(ui, server)
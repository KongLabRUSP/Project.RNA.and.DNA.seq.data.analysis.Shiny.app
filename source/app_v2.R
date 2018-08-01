# |----------------------------------------------------------------------------------|
# | Project: RNA/DNA-seq pipeline                                                    |
# | Script:   Shiny app to streamline RNA and DNA seq analysis                       |
# | Authors:  Marissa (Meinizi) Zheng, Davit Sargsyan                                |   
# | Created:  07/22/2018                                                             |
# | Modified: 07/27/2018: Added folder path (using **shinyFiles** package), and      |
# |                      ability to run batch file by clicking **Run Program** button|  
# |           08/01/2018: app_v2 can now run NextFlow scripts (NextFlow must be      |
# |                       installed on a Linux machine first). A dummy 'pipeline.nf' |
# |                       script is created and tested for passing parameter values. |
# |----------------------------------------------------------------------------------|
# Reference: https://www.nextflow.io/docs/latest/getstarted.html
options(stringsAsFactors = FALSE)

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
                                                                radioButtons(inputId = "isSingleEnd",
                                                                             label = "Select Single or Pair End",
                                                                             choices = list("Single End" = "true", 
                                                                                            "Pair End" = "false"),
                                                                             selected = "true"),
                                                                actionButton(inputId = "go",
                                                                             label = "Run Program")),
                                                   mainPanel(tags$h3("Working Directory"),
                                                             verbatimTextOutput(outputId = "wd"),
                                                             tags$h3("File Directory"),
                                                             verbatimTextOutput(outputId = "fileDir"),
                                                             tags$h3("Program"),
                                                             verbatimTextOutput("directorypath"))),
                                           tabItem(tabNamvolumes = "dna",
                                                   sidebarPanel( ),
                                                   mainPanel()))))

server <- function(input, output, session) {
  # Specify folder containing FastQ file.
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
    pgm <- paste(getwd(),
                 '/nextflow run -process.echo true ',
                 getwd(),
                 '/rna-seq-genecount.nf',
                 # '/pipeline.nf',
                 ' --SingleEnd=\"',
                 input$isSingleEnd,
                 '\" --reads=\"',
                 parseDirPath(roots = volumes, 
                              selection = input$directory),
                 '/*.gz\" > ',
                 getwd(),
                 '/log.txt',
                 sep = "")

    # Print to Shiny textbox
    cat(pgm)
  })
  
  observeEvent(input$go, {
    pgm <- paste(getwd(),
                 '/nextflow run -process.echo true ',
                 getwd(),
                 '/rna-seq-genecount.nf',
                 # '/pipeline.nf',
                 ' --SingleEnd=\"',
                 input$isSingleEnd,
                 '\" --reads=\"',
                 parseDirPath(roots = volumes, 
                              selection = input$directory),
                 '/*.gz\" > ',
                 getwd(),
                 '/log.txt',
                 sep = "")

    # Run the program
    system(pgm)
  })
}

shinyApp(ui, server)
# |----------------------------------------------------------------------------------|
# | Project: RNA/DNA-seq pipeline                                                    |
# | Script: application                                                              |
# | Authors: Marissa (Meinizi) Zheng, Davit Sargsyan                                 |   
# | Created: 07/22/2018                                                              |
# | Modified:  7/26/2018                                                                      |
# |----------------------------------------------------------------------------------|
options(stringsAsFactors = FALSE)

require(shiny)
require(data.table)
require(DT)
library(shinydashboard)
library(shinythemes)

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
                                                   sidebarPanel(textInput(inputId = "URL", label = "Enter FASTQ file directory below£º"),
                                                                radioButtons(inputId = "isSingleEnd",
                                                                             label = "Select Single or Pair End",
                                                                             choices = list("Single End" = "true", "Pair End" = "false"),
                                                                             selected = "Single End")),
                                                   mainPanel(actionButton(inputId = "go",
                                                                          label = "Go"))),
                                           tabItem(tabName = "dna",
                                                   sidebarPanel( ),
                                                   mainPanel()))))

server <- function(input, output, session) {
  observeEvent(input$button, {
               shell.exec( "C:/Users/Pro/Documents/RNAseq_Pipeline/rna-seq-genecount.nf --prams.reads input$URL  --params.SingleEnd input$isSingleEnd " )})
  
  

}

shinyApp(ui, server)

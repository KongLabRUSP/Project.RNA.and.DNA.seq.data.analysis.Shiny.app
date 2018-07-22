# |----------------------------------------------------------------------------------|
# | Project: RNA/DNA-seq pipeline                                                    |
# | Script: application                                                              |
# | Authors: Marissa (Meinizi) Zheng, Davit Sargsyan                                 |   
# | Created: 07/22/2018                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
options(stringsAsFactors = FALSE)

require(shiny)
require(data.table)
require(DT)
library(shinydashboard)
library(shinythemes)

ui <- dashboardPage(dashboardHeader(title = "Shiny ICD",
                                    dropdownMenu(type = "notifications",
                                                 notificationItem(text = "5 new users today",
                                                                  icon = icon("users")),
                                                 notificationItem(text = "12 items delivered",
                                                                  icon = icon("truck"),
                                                                  status = "success"),
                                                 notificationItem(text = "Server load at 86%",
                                                                  icon = icon("exclamation-triangle"),
                                                                  status = "warning"))),
                    dashboardSidebar(sidebarMenu(menuItem(text = "Dashboard", 
                                                          tabName = "dashboard", 
                                                          icon = icon("dashboard")),
                                                 menuItem(text = "Mapping", 
                                                          tabName = "mapping", 
                                                          icon = icon("th")),
                                                 menuItem(text = "Convert", 
                                                          tabName = "convert", 
                                                          icon = icon("th")))),
                    dashboardBody(tabItems(tabItem(tabName = "dashboard",
                                                   h2("Hello1")),
                                           tabItem(tabName = "mapping",
                                                   sidebarPanel(radioButtons(inputId = "dataset",
                                                                             label = "Select List",
                                                                             choices = c("Diagnoses",
                                                                                         "Procedures"),
                                                                             selected = "Diagnoses"),
                                                                selectInput(inputId = "icd9_version",
                                                                            label = "ICD-9 Version",
                                                                            choices = 1:10),
                                                                uiOutput(outputId = "sub1In"),
                                                                uiOutput(outputId = "sub2In"),
                                                                uiOutput(outputId = "sub3In")),
                                                   mainPanel(DT:: dataTableOutput("tbl"),
                                                             br(),
                                                             actionButton(inputId = "do",
                                                                          label = "Save Selection"),
                                                             br(),
                                                             DT:: dataTableOutput("tbl2"),
                                                             br(),
                                                             downloadLink(outputId = "downloadData",
                                                                          label = "Download Selected Rows"),
                                                             br(),
                                                             downloadLink(outputId = "downloadMap",
                                                                          label = "Download Map of Selected Rows"))),
                                           tabItem(tabName = "convert",
                                                   sidebarPanel(fileInput(inputId = "browseMap",
                                                                          label = "Select Mapping File",
                                                                          multiple = FALSE),
                                                                fileInput(inputId = "browseData",
                                                                          label = "Select ICD Data File",
                                                                          multiple = FALSE),
                                                                uiOutput(outputId = "idColIn"),
                                                                uiOutput(outputId = "icdColsIn")),
                                                   mainPanel(DT:: dataTableOutput("tblMap"),
                                                             br(),
                                                             DT:: dataTableOutput("tblICD"))))))

server <- function(input, output, session) {
 output$tblMap<- DT::renderDT({
    validate(need(input$browseMap != "", ""))
    ne <- new.env()
    dt3 <- fread(input$browseMap$datapath)
    DT::datatable(head(dt3, 3),
                  options = list(pageLength = 10),
                  selection = list(mode = "multiple"),
                  rownames = FALSE)
  })
  
  output$tblICD <- DT::renderDT({
    validate(need(input$browseData != "", ""))
    ne <- new.env()
    fname <- load(file = input$browseData$datapath,
                  envir = ne)
    dt2 <- ne[[fname]]
    DT::datatable(head(dt2, 20),
                  options = list(pageLength = 10),
                  selection = list(mode = "multiple"),
                  rownames = FALSE)
  })
}

shinyApp(ui, server)
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# ======================================= RNA-seq ===============================
library(shiny)

ui <- fluidPage(
   
   titlePanel("RNA-seq Analysis Pipeline"),
   
   tabsetPanel(
     tabPanel("FASTQ to Count", 
              wellPanel(textInput("URL", "Put FASTQ File Directory£º"),
                        radioButtons("isSingleEnd", label = "Choose Single or Pair End:",
                        choices = list("Single End" = TRUE, "Pair End" = FALSE)),
                        actionButton("button", label = "Go"))
       
     ),
     tabPanel("tab2", "contents")
   )
)


server <- function(input, output) {
  
   observeEvent(input$button, {
     shell.exec( "C:/Users/Pro/Documents/RNAseq_Pipeline/rna-seq-genecount.nf --prams.reads input$URL  --params.SingleEnd input$isSingleEnd " )
   })
  
}

# Run the application 
shinyApp(ui = ui, server = server)






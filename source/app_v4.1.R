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
# |           08/13/2018: app_v3 added DNA program.                                  |
# |           08/20/2018: app_v4 RNA DEGseq analysis                                 |
# |           08/27/2018: app_v4 add VennDiagram and heatmap to DEGseq analysis, add |
# |                       exploratory analysis with more general experiment design   |
# |           09/03/2018: can run Deseq2 DE analysis, generate/ export result table  |
# |                       and MA plot                                                |
# |           09/10/2018: add dispersion estimation plot and count plot in Deseq2    |
# |                       analysis, add multiple data normalization method in        |
# |                       Exploratory Analysis Pannel and did some layout polish     |
# |----------------------------------------------------------------------------------|
# Reference: https://www.nextflow.io/docs/latest/getstarted.html
options(stringsAsFactors = FALSE)
options(shiny.maxRequestSize=30*1024^2) 
# install.packages("pheatmap")

require(shiny)
library(shinydashboard)
library(shinythemes)
require(shinyFiles)
library(shinyWidgets)

require(stringr)
require(readxl)
require(DT)
require(data.table)
require(MASS)
require(ggplot2)
require(DESeq2)
require(BiocParallel)
require(DEGseq)
require(knitr)
require(ggdendro)
require(VennDiagram)
library(RColorBrewer)
library(pheatmap)





# Determine the OS
sysOS <- Sys.info()[['sysname']]

# Specify root folder----
if (sysOS == "Linux") {
  volumes <- c("Home" = "/home/administrator/",
               "Data Storage" = "/datastorage/")
} else if (sysOS == "Windows") {
  volumes <- c("Local Disk (C:)" = "C:/")
}

ui <- dashboardPage(dashboardHeader(title = "NGS Pipeline",
                                    dropdownMenu(type = "notifications",
                                                 notificationItem(text = "",      
                                                                  icon = icon("exclamation-triangle"),
                                                                  status = "warning"))),
                    dashboardSidebar(sidebarMenu(menuItem(text = "Introduction", 
                                                          tabName = "introduction", 
                                                          icon = icon("dashboard")),
                                                 menuItem(text = "Generate Count Table", 
                                                          tabName = "count", 
                                                          icon = icon("th")),
                                                 menuItem(text = "RNA-seq Analysis", 
                                                          tabName = "rna-seq_analysis", 
                                                          icon = icon("th")),
                                                 menuItem(text = "DNA MethylSeq Analysis",
                                                          tabName = "dna-seq_analysis",
                                                          icon = icon("th")),
                                                 menuItem(text = "DNA vs RNA",
                                                          tabName = "dna_vs_rna",
                                                          icon = icon("th")))),
                    dashboardBody(tabItems(
                                           ## ---------- Introduction ------------
                                           tabItem(tabName = "dashboard",
                                                   h2("Workflow Overview & Instructions")),
                                           
                                           
                                           
                                           ## -------- generate count table -------------------
                                           tabItem(tabName = "count",
                                                   sidebarPanel(
                                                                wellPanel(tags$h5("Basic Setting:"),
                                                                          shinyDirButton(id = "directory", 
                                                                                         label = "Folder select", 
                                                                                         title = "Please select a folder"),
                                                                          radioButtons(inputId = "isSingleEnd",
                                                                                       label = "Select Single or Pair End",
                                                                                       choices = list("Single End" = "true", 
                                                                                                      "Pair End" = "false"),
                                                                                       selected = "true")),
                                                                wellPanel(tags$h5("Advanced Setting:"),
                                                                          radioButtons(inputId = "advanced_para",
                                                                                       label = NULL,
                                                                                       inline = TRUE,
                                                                                       choices = list("Hide" = "basic", 
                                                                                                      "Expand" = "advanced"),
                                                                                       selected = "basic"),
                                                                          conditionalPanel(condition = "input.advanced_para == 'advanced'",
                                                                                           textInput(inputId = "job_name",
                                                                                                     label = "job name (mmddy + x job)",
                                                                                                     value = paste0(substr(format(Sys.Date(),"%m%d%y"),1,4),
                                                                                                                    substr(format(Sys.Date(),"%m%d%y"),6,6),
                                                                                                                    "1")),
                                                                                           selectInput(inputId = "ntasks",
                                                                                                       label = "ntasks",
                                                                                                       choices = c(1),
                                                                                                       multiple = FALSE),
                                                                                           selectInput(inputId = "cpus_per_task",
                                                                                                       label = "cpus per task",
                                                                                                       choices = c(1:32),
                                                                                                       multiple = FALSE,
                                                                                                       selected = 2),
                                                                                           selectInput(inputId = "mem",
                                                                                                       label = "RAM",
                                                                                                       choices = c(paste0(c(16, 32, 64, 120, 192), "GB")),
                                                                                                       multiple = FALSE,
                                                                                                       selected = "64GB"),
                                                                                           selectInput(inputId = "time",
                                                                                                       label = "total run time limit (HH:MM:SS)",
                                                                                                       choices = c(paste0(c(2, 4, 6, 8, 12, 24),":00:00")),
                                                                                                       multiple = FALSE,
                                                                                                       selected = "2:00:00"),
                                                                                           textInput(inputId = "output",
                                                                                                     label = "STDOUT output file",
                                                                                                     value = "slurm.%N.%j.out"),
                                                                                           textInput(inputId = "error",
                                                                                                     label = "STDERR output file (optional)",
                                                                                                     value = "slurm.%N.%j.err"),
                                                                                           selectInput(inputId = "partition",
                                                                                                       label = "partition",
                                                                                                       choices = c("p_kongt_1", "main"),
                                                                                                       multiple = FALSE,
                                                                                                       selected = "p_kongt_1"))),
                                                                actionButton(inputId = "gorna",
                                                                             label = "Run RNA Program"),
                                                                actionButton(inputId = "godna",
                                                                             label = "Run DNA Program")),
                                                   mainPanel(tags$h3("Working Directory"),
                                                             verbatimTextOutput(outputId = "wd"),
                                                             tags$h3("File Directory"),
                                                             verbatimTextOutput(outputId = "fileDir"),
                                                             tags$h3("Program"),
                                                             verbatimTextOutput("directorypath"))),
                                                                            
                                           ## ---------------- rna seq analysis -----------------
                                           tabItem(tabName = "rna-seq_analysis",
                                                   tabsetPanel(
                                                     # ------- read in data -------------
                                                     tabPanel("Read In Data", 
                                                              wellPanel(
                                                                fluidRow(
                                                                  column(4,
                                                                         textInput(inputId = "project_name",
                                                                                   label = "Enter project name (no space) below:",
                                                                                   value = "Enter project name"))),
                                                                fluidRow(
                                                                  column(6,
                                                                         fileInput(inputId = "rna_count",
                                                                                   label = "Count Table File Select",
                                                                                   accept = c(".cvs")))),
                                                                fluidRow(
                                                                  column(6,
                                                                         fileInput(inputId = "rna_info",
                                                                                   label = "Sample Information File Select (for Deseq2)",
                                                                                   accept = c(".cvs"))))),
                                                              
                                                               
                                                              wellPanel(
                                                                tags$h4("View Count Table:"),
                                                                DT::dataTableOutput(outputId = 'display_rna_count')),
                                                              wellPanel(
                                                                tags$h4("View Sample Information Table:"),
                                                                DT::dataTableOutput(outputId = 'display_rna_info'))),
                                                     
                                                     
                                                     # ------- Exploratory Analysis ---------
                                                     tabPanel("Exploratory Analysis", 
                                                              wellPanel(
                                                                tags$h3("Data Preparation"),
                                                                fluidRow(
                                                                  column(4,
                                                                         pickerInput(inputId = "rna_expl_sample_choices",
                                                                                     label = "Select samples (same order as you want it to show in the plots):",
                                                                                     choices = NULL,
                                                                                     multiple = TRUE,
                                                                                     selected = NULL,
                                                                                     options = list(`actions-box` = TRUE))
                                                                         
                                                                         # selectInput(inputId = "rna_expl_sample_choices",
                                                                         #             label = "Select samples (same order as you want it to show in the plots):",
                                                                         #             choices = NULL,
                                                                         #             multiple = TRUE,
                                                                         #             selected = NULL)
                                                                         ),
                                                                  column(4,
                                                                         selectizeInput(inputId = "rna_expl_groupby_choices",
                                                                                        label = "Group by (select one, optional):",
                                                                                        choices = NULL,
                                                                                        multiple = TRUE,
                                                                                        selected = NULL))),
                                                                fluidRow(
                                                                  column(2,
                                                                         actionButton(inputId = "generate_result_expl",
                                                                                      label = "Generate Result")),
                                                                  column(2,
                                                                         actionButton(inputId = "export_plots_expl",
                                                                                      label = "Export Generated Plots"))),
                                                                fluidRow(
                                                                  column(10,
                                                                         verbatimTextOutput(outputId = "export_dir_expl"))
                                                                )),
                                                              wellPanel(
                                                                tags$h3("Count Histgram"),
                                                                plotOutput(outputId = "display_count_hist")),
                                                              
                                                              wellPanel(
                                                                tags$h3("Between-sample Distribution: Boxplot"),
                                                                fluidRow(
                                                                  column(6,
                                                                         plotOutput(outputId = "display_boxplot")))),
                                                              wellPanel(
                                                                fluidRow(
                                                                  column(6,
                                                                         wellPanel(
                                                                           tags$h3("Heatmap Sampel-to-sample Distance"),
                                                                           fluidRow(imageOutput(outputId = "display_heatmap_expl")))),
                                                                  column(6,
                                                                         wellPanel(
                                                                           tags$h3("PCA Plot"),
                                                                           fluidRow(plotOutput(outputId = "display_pca_expl"))))))
                                                              
                                                              
                                                              
                                                     ),
                                                     
                                                     
                                                     # ------------- DEGseq --------------
                                                     tabPanel("DE Analysis: DEGseq",
                                                              # Data preparation
                                                              wellPanel(
                                                                tags$h3("Data Preparation"),
                                                                fluidRow(
                                                                  column(4,
                                                                         selectInput(inputId = "dge_choices",
                                                                                     label = "Select 3 treatments in order: trt1 trt2 trt3",
                                                                                     choices = NULL,
                                                                                     multiple = TRUE,
                                                                                     selected = NULL)),
                                                                  column(3,
                                                                         actionButton(inputId = "generate_count_table",
                                                                                      label = "Generate Count Table"))),
                                                                fluidRow(
                                                                  column(4,
                                                                         textInput(inputId = "row_sum",
                                                                                   label = "Remove if total across 3 samples is <",
                                                                                   value = "10")),
                                                                  column(3,
                                                                         actionButton(inputId = "trim",
                                                                                      label = "Trim Count Table"))),
                                                                fluidRow(
                                                                  column(6,
                                                                         tags$h4("Count Table Summary:"),
                                                                         verbatimTextOutput("count_table_txtbox")),
                                                                  column(6,
                                                                         tags$h4("Trimmed Count Table Summary:"),
                                                                         verbatimTextOutput("trim_txtbox"),
                                                                         verbatimTextOutput("trimmed_count_table_txtbox")))),
                                                              
                                                              # DE analysis
                                                              wellPanel(
                                                                tags$h3("DE Analysis"),
                                                                fluidRow(
                                                                  column(3,
                                                                         selectInput(inputId = "q_value",
                                                                                     label = "q-value(Storey et al. 2003) <",
                                                                                     choices = c(0.01, 0.05, 0.1, 0.25),
                                                                                     selected = 0.01)),
                                                                  column(3,
                                                                         selectInput(inputId = "fold_change",
                                                                                     label = "Obs(log2(fold change)) >=",
                                                                                     choices = c(0.3, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                                                                     selected = 1)),
                                                                  column(5,
                                                                         actionButton(inputId = "run_DEGexp",
                                                                                      label = "RUN DEGexp"),
                                                                         p("Please wait for the result tables show up..."))),
                                                                
                                                                fluidRow(
                                                                  column(3,
                                                                         actionButton(inputId = "export_table_trt2_trt1",
                                                                                      label = "Export trt2_trt1 result table")),
                                                                  column(9,
                                                                         # show export path
                                                                         verbatimTextOutput(outputId = "export_dir1"))),
                                                                
                                                                fluidRow(
                                                                  column(3,
                                                                         actionButton(inputId = "export_table_trt3_trt2",
                                                                                      label = "Export trt3_trt2 result table")),
                                                                  column(9,
                                                                         # show export path
                                                                         verbatimTextOutput(outputId = "export_dir2"))),
                                                                wellPanel(
                                                                  tags$h4("trt2-trt1 Result Table:"),
                                                                  DT::dataTableOutput(outputId = "result_table1")),
                                                                
                                                                wellPanel(
                                                                  tags$h4("trt3-trt2 Result Table:"),
                                                                  DT::dataTableOutput(outputId = "result_table2")),
                                                                
                                                                wellPanel(
                                                                  tags$h4("MA Plot:"),
                                                                  fluidRow(
                                                                    column(3,
                                                                           textInput(inputId = "ma_title1",
                                                                                     label = "Enter trt2-trt1 MA plot title below:",
                                                                                     value = "")),
                                                                    column(3,
                                                                           textInput(inputId = "ma_title2",
                                                                                     label = "Enter trt3-trt2 MA plot title below:",
                                                                                     value = "")),
                                                                    column(2,
                                                                           actionButton(inputId = "generate_maplot",
                                                                                        label = "Generate MA Plot")),
                                                                    column(2,
                                                                           actionButton(inputId = "export_maplot",
                                                                                        label = "Export MA Plot"))),
                                                                  
                                                                  fluidRow(
                                                                    verbatimTextOutput(outputId = "export_maplot_dir")),
                                                                  
                                                                  fluidRow(
                                                                    column(6,
                                                                           wellPanel(
                                                                             plotOutput(outputId = "result_plot1"))),
                                                                    column(6,
                                                                           wellPanel(
                                                                             plotOutput(outputId = "result_plot2")))),
                                                                  fluidRow(
                                                                    column(6,
                                                                           verbatimTextOutput(outputId = "result_kable1")),
                                                                    column(6,
                                                                           verbatimTextOutput(outputId = "result_kable2"))))),

                                                              
                                                              wellPanel(
                                                                tags$h3("Changes in Gene Expression"),
                                                                fluidRow(
                                                                  column(3,
                                                                         actionButton(inputId = "generate_up.dn_dn.up_table_and_heatmap",
                                                                                      label = "Generate Results")),
                                                                  column(3,
                                                                         actionButton(inputId ="export_change_in_gene_exp",
                                                                                      label = "Export Table and Plots"))),
                                                                fluidRow(
                                                                  verbatimTextOutput(outputId ="up.dn_number_txtbox")),
                                                                fluidRow( 
                                                                  verbatimTextOutput(outputId ="display_change_in_gene_exp_dir")),
                                                                fluidRow(
                                                                  column(4,
                                                                         plotOutput(outputId = "display_venn_diagram1")),
                                                                  column(4,
                                                                         plotOutput(outputId = "display_venn_diagram2"))),
                                                                fluidRow(
                                                                  column(6,
                                                                         wellPanel(
                                                                           plotOutput(outputId = "display_up.dn_dn.up_heatmap")
                                                                         )),
                                                                  column(6,
                                                                         wellPanel(DT::dataTableOutput(outputId = "display_up.dn_dn.up_table")))))),
                                                     
                                                     
                                                     # --------------- DEseq2 ----------------
                                                     tabPanel("DE Analysis: DEseq2",
                                                              
                                                              # Data Preparation
                                                              wellPanel(
                                                                tags$h3("Data Preparation"),
                                                                
                                                                fluidRow(
                                                                  column(3,
                                                                         actionButton(inputId = "generate_count_matrix",
                                                                                      label = "Generate Count Matrix"))),
                                                                fluidRow(tags$h3("")),
                                                                wellPanel(tags$h4("Count Matrix Summary:"),
                                                                          fluidRow(column(12,
                                                                                          verbatimTextOutput(outputId = "display_summary_of_count_matrix")))),
                                                                
                                                                
                                                                
                                                                fluidRow(
                                                                  column(6,
                                                                         actionButton(inputId = "generate_design_matrix",
                                                                                      label = "Generate Design Matrix"))),
                                                                fluidRow(tags$h3("")),
                                                                wellPanel(
                                                                  tags$h4("Design Matrix Summary:"),
                                                                  fluidRow(column(8,
                                                                                  verbatimTextOutput(outputId = "display_summary_of_design_matrix")))),
                                                                
                                                                fluidRow(tags$h3("")),
                                                                fluidRow(
                                                                  column(4,
                                                                         selectInput(inputId = "deseq2_design_formula",
                                                                                     label = "Select single or multi-factor design formula:",
                                                                                     choices = NULL,
                                                                                     multiple = TRUE,
                                                                                     selected = NULL)),
                                                                  column(2,
                                                                         radioButtons(inputId = "multi_factor_option",
                                                                                      label = "Select one if multi-factor:",
                                                                                      choices = list("Select None" = "none",
                                                                                                     "As Group" = "as_group",
                                                                                                     "With Interaction" = "w_interaction", 
                                                                                                     "Without Interaction" = "o_interaction"),
                                                                                      selected = NULL)),
                                                                  column(3,
                                                                         actionButton(inputId = "generate_deseqdt_from_matrix",
                                                                                      label = "Generate DESeqDataSet From Matrix"))),
                                                                fluidRow(tags$h3("")),
                                                                fluidRow(
                                                                  column(4,
                                                                         textInput(inputId = "row_sum_deseq2",
                                                                                   label = "Remove if total across samples is <",
                                                                                   value = "10")),
                                                                  column(3,
                                                                         actionButton(inputId = "trim_deseqdt",
                                                                                      label = "Trim DESeqDataSet"))),
                                                                
                                                                wellPanel(
                                                                  fluidRow(
                                                                    column(6,
                                                                           tags$h4("DESeqDataSet Summary:"),
                                                                           verbatimTextOutput("display_formula"),
                                                                           verbatimTextOutput(outputId = "display_deseqdt")),
                                                                    column(6,
                                                                           tags$h4("Trimmed DESeqDataSet Summary:"),
                                                                           verbatimTextOutput(outputId = "display_trim_deseqdt_result"),
                                                                           verbatimTextOutput(outputId = "display_trimmed_deseqdt"))))),
                                                              
                                                              # DE analysis
                                                              wellPanel(
                                                                tags$h3("DE Analysis"),
                                                                fluidRow(column(4,
                                                                                selectInput(inputId = "deseq2_p_value",
                                                                                            label = "FDR adjusted p-value <",
                                                                                            choices = c(0.01, 0.05, 0.1, 0.25),
                                                                                            selected = 0.01)),
                                                                         column(4,
                                                                                selectInput(inputId = "deseq2_fold_change",
                                                                                            label = "Obs(log2(fold change)) >=",
                                                                                            choices = c(0.3, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                                                                                            selected = 1)),
                                                                         column(3,
                                                                                actionButton(inputId = "run_deseq",
                                                                                             label = "Run DESeq"))),
                                                                fluidRow(column(4,
                                                                                selectInput(inputId = "select_contrast_add",
                                                                                            label = "Select names of the log2 fold changes to add:",
                                                                                            choices = NULL,
                                                                                            multiple = TRUE,
                                                                                            selected = NULL)),
                                                                         column(4,
                                                                                selectInput(inputId = "select_contrast_substract",
                                                                                            label = "Select names of the log2 fold changes to substract (optional):",
                                                                                            choices = NULL,
                                                                                            multiple = TRUE,
                                                                                            selected = NULL)),
                                                                         column(2,
                                                                                actionButton(inputId = "generate_deseq_results",
                                                                                             label = "Generate Results")),
                                                                         column(2,
                                                                                actionButton(inputId = "export_deseq_results",
                                                                                             label = "Export Results"))),
                                                                fluidRow(verbatimTextOutput(outputId = "deseq2_results_dir")),
                                                                wellPanel(tags$h4("Result Table:"),
                                                                          DT::dataTableOutput(outputId = "display_dtf_res_contrast")),
                                                                wellPanel(tags$h4("MA Plot:"),
                                                                          fluidRow(column(3,
                                                                                          textInput(inputId = "enter_deseq2_ma_title",
                                                                                                    label = "Enter MA plot title below:",
                                                                                                    value = ""))),
                                                                          fluidRow(column(6,
                                                                                          plotOutput(outputId = "display_deseq2_ma")),
                                                                                   column(6,
                                                                                          tableOutput(outputId = "deseq2_sign_number")))),
                                                                wellPanel(tags$h4("Estimated Dispersion:"),
                                                                          fluidRow(column(6,
                                                                                          imageOutput(outputId = "deseq2_estimated_dispersion")))),
                                                                wellPanel(tags$h4("Plot Counts"),
                                                                          fluidRow(column(3,
                                                                                          textInput(inputId = "enter_gene_name",
                                                                                                    label = "Enter one gene name:")),
                                                                                   column(3,
                                                                                          selectInput(inputId = "intgroup",
                                                                                                      label = "Select the name of interested group(s):",
                                                                                                      choices = NULL,
                                                                                                      multiple = TRUE,
                                                                                                      selected = NULL)),
                                                                                   column(3,
                                                                                          actionButton(inputId = "generate_count_plot",
                                                                                                       label = "Generate Count Plot"))),
                                                                          fluidRow(column(6,
                                                                                          plotOutput(outputId = "deseq2_count_plot")))))))),
                                           
                                           ## ---------------- dna methylseq analysis ----------------
                                           tabItem(tabName = "dna-seq_analysis",
                                                   tabsetPanel(
                                                     tabPanel("Exploratory Analysis"),
                                                     tabPanel("DE Analysis: DDS ")
                                                   )),
                                           
                                           ## --------------- dna vs rna -----------------
                                           tabItem(tabName = "dna_vs_rna",
                                                   tabsetPanel(
                                                     tabPanel("DNA vs RNA")
                                                     # tabPanel()
                                                   )))))







## ------------- server -------------------
server <- function(input, output, session) {
  ##  ------------ upstream ---------------------
  # Specify folder containing FastQ file.
  # NOTE: all the files in the folder will be passed to the pipeline
  rv <- reactiveValues()
  
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
    rv$pgmrna <- paste0("sbatch --job-name=",
                        input$job_name,
                        " --ntasks=",
                        input$ntasks,
                        " --cpus-per-taks=",
                        input$cpus_per_task,
                        " --mem=",
                        input$mem,
                        " --time=",
                        input$time,
                        " --output=",
                        input$output,
                        " --error=",
                        input$error,
                        ifelse(input$partition == "main", 
                               "",
                               paste0(" --partition=",
                                      input$partition)),
                        " module use /projects/community/modulefiles", 
                        " module load nextflow",
                        " export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)",
                        " export NXF_OPTS=\'-Xms1g -Xmx4g\'", 
                        " export NF_Work_Dir=\"/scratch/${USER}/NFWorkDir/${PWD}/work\"",
                        " export IS_SINGLE=\"",
                        input$isSingleEnd,
                        "\"",
                        " export IN_FILES=\"",
                        ifelse(input$isSingleEnd == "true",
                               paste0(parseDirPath(roots = volumes, 
                                                   selection = input$directory),
                                      "/*.gz\""),
                               paste0(parseDirPath(roots = volumes, 
                                                   selection = input$directory),
                                       "/*_R{1,2}*.fastq.gz\"")),
                        " mkdir -p $NF_Work_Dir",
                        " date",
                        " SECONDS=0",
                        " srun nextflow run rna-seq.nf --SingleEnd=$IS_SINGLE --reads=$IN_FILES -w $NF_Work_Dir -with-trace -with-report ${SLURM_JOB_PARTITION}_${SLURM_JOB_NODELIST}_${SLURM_JOB_ID}_RNA-seq-nf.html \\-with-timeline ${SLURM_JOB_PARTITION}_${SLURM_JOB_NODELIST}_${SLURM_JOB_ID}_RNA-seq-nf-timeline.html -resume",
                        " date",
                        " duration=$SECONDS",
                        " echo \"$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed.\"",
                        " touch \"${SLURM_JOB_PARTITION}-${SLURM_JOB_NODELIST}-${SLURM_JOB_ID}-$(($duration / 3600 ))h_$(($(($duration % 3600)) / 60 ))m_$(($duration % 60))s.time\""
                        )
      
      
    
    rv$pgmdna <- paste0("sbatch --job-name=",
                        input$job_name,
                        " --ntasks=",
                        input$ntasks,
                        " --cpus-per-taks=",
                        input$cpus_per_task,
                        " --mem=",
                        input$mem,
                        " --time=",
                        input$time,
                        " --output=",
                        input$output,
                        " --error=",
                        input$error,
                        ifelse(input$partition == "main", 
                               "",
                               paste0(" --partition=",
                                      input$partition)),
                        " module use /projects/community/modulefiles", 
                        " module load nextflow",
                        " export NXF_OPTS=\'-Xms1g -Xmx4g\'",
                        " export NF_Work_Dir=\"/scratch/${USER}/NFWorkDir/${PWD}/work\"",
                        " mkdir -p $NF_Work_Dir",
                        " IS_SINGLE=\"",
                        input$isSingleEnd,
                        "\"",
                        " IN_FILES=\"",
                        ifelse(input$isSingleEnd == "true",
                               paste0(parseDirPath(roots = volumes, 
                                                   selection = input$directory),
                                      "/*.gz\""),
                               paste0(parseDirPath(roots = volumes, 
                                                   selection = input$directory),
                                      "/*_R{1,2}*.fastq.gz\"")),
                        # add the code inbetween?
                        " mkdir -p $NF_Work_Dir",
                        " date",
                        " SECONDS=0",
                        " srun nextflow run methyl-seq.nf -w $NF_Work_Dir --SingleEnd=$IS_SINGLE --reads=$IN_FILES -with-trace -with-report DNA-methyl.html  -with-timeline DNA-methyl-timeline.html -resume",
                        " EXIT_STATUS=$?",
                        " date",
                        " echo \"$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed.\"",
                        " touch \"${SLURM_JOB_PARTITION}-${SLURM_JOB_NODELIST}-${SLURM_JOB_ID}-$(($duration / 3600 ))h_$(($(($duration % 3600)) / 60 ))m_$(($duration % 60))s.time\"",
                        " if  [ $EXIT_STATUS  -eq  0 ]; then echo \" Pipeline completed. Removing WorkDir files\" rm -rf $NF_Work_Dir else echo \"Pipeline not completed.\" fi" 
                        )
    
    # Print to Shiny textbox
    cat(paste('RNA Program:',rv$pgmrna,'\n','DNA Program:', rv$pgmdna, sep = '\n'))
  })
  
  observeEvent(input$gorna, {
    # trigger nextflow to run rv$pgmrna or rv$pgmdna
    
  })
  
  observeEvent(input$godna, {
    pgmdna <- paste(getwd(),
                    '/nextflow run -process.echo true ',
                    getwd(),
                    # '/dna-methyl-nf-v2.nf',
                    '/pipeline.nf',
                    ' --SingleEnd=\"',
                    input$isSingleEnd,
                    '\" --reads=\"',
                    parseDirPath(roots = volumes, 
                                 selection = input$directory),
                    '/*.gz\" > ',
                    getwd(),
                    '/dnalog.txt',
                    sep = "")
    # Run the DNA program
    system(pgmdna)
  })
  
  
  
  
  
  ## ---------------------- downstream: RNA seq ---------------------
  ## ------------  Read in Data and update SelectInput -------------
  
  ## read in, display count table, update SelectInput choice
  observeEvent(input$rna_count,{
    
    
    infile_rna_count <- input$rna_count
    
    # process if count table is selected
    if(!is.null(infile_rna_count)){
      
      # read in data
      tb <- fread(infile_rna_count$datapath, skip = 1)
      
      if(colnames(tb)[1] != "Geneid"){
        showNotification("Error: count table is incorrect",
                         type = "error",
                         duration = 15)
      }else{
        trt.names <- str_extract(colnames(tb)[7:ncol(tb)], "[^.]+")
        colnames(tb)[7:ncol(tb)] <- trt.names
        colnames(tb)[1] <- "gene"
        tb <- as.data.frame(tb)
        rv$ct <- tb
        
        # display count table
        output$display_rna_count = DT::renderDataTable({
          DT::datatable(rv$ct,
                        options = list(scrollX = TRUE))})
        
        # update SelectInput
        updateSelectInput(session, 
                          inputId = "dge_choices", 
                          label = "Select 3 treatments in order: trt1 trt2 trt3", 
                          choices = colnames(rv$ct)[7: ncol(rv$ct)], 
                          selected = NULL)
        
        updatePickerInput(session, 
                          inputId = "rna_expl_sample_choices", 
                          label = "Select samples (same order as you want it to show in the plots):", 
                          choices = colnames(rv$ct)[7: ncol(rv$ct)], 
                          selected = NULL)
      
        # updateSelectInput(session, 
        #                   inputId = "rna_expl_sample_choices", 
        #                   label = "Select samples (same order as you want it to show in the plots):", 
        #                   choices = colnames(rv$ct)[7: ncol(rv$ct)], 
        #                   selected = NULL)
      }
      
    }
    
  })

  
  ## read in info table and update SelectInput choice
  observeEvent(input$rna_info,{
    infile_rna_count <- input$rna_count
    infile_rna_info <- input$rna_info
    # process if count and info are read in
    if(!is.null(infile_rna_count) & !is.null(infile_rna_info)){
      
      tb <- read.csv(infile_rna_info$datapath, header = TRUE)
      if(nrow(tb) == length(7:ncol(rv$ct))){
        tb[, colnames(tb)] <- lapply(tb, as.factor)
        tb <- cbind("Samples" = colnames(rv$ct)[7:ncol(rv$ct)], tb)
        rv$info <- tb
        
        output$display_rna_info = DT::renderDataTable({
          DT::datatable(rv$info,
                        options = list(scrollX = TRUE))})
        
        updateSelectizeInput(session, 
                             inputId = "rna_expl_groupby_choices", 
                             label = "Group by (select one, optional):", 
                             choices = colnames(rv$info)[2: ncol(rv$info)], 
                             selected = NULL)
        updateSelectInput(session, 
                          inputId = "deseq2_design_formula", 
                          label = "Select single or multi-factor design formula:", 
                          choices = colnames(rv$info)[2: ncol(rv$info)], 
                          selected = NULL)
        updateSelectInput(session, 
                          inputId = "intgroup", 
                          label = "Select the name of interested group(s):", 
                          choices = colnames(rv$info)[2: ncol(rv$info)], 
                          selected = NULL)
      }else{
        showNotification("Error: row number of your info table should be the same as your sample number",
                         type = "error",
                         duration = 15)
      }
      
    }
  })
  
  
  
  ## ------------ Exploratory Analysis -------------
  
  ## generate design matrix log2(base+1)
  observeEvent(input$rna_expl_sample_choices,{
    tmp <- rv$ct[ ,c("gene",input$rna_expl_sample_choices)]
    rv$expl_matrix <- cbind(tmp[, 1],log2(tmp[, -1] +1))
    rv$expl_matrix_long <- melt(data = rv$expl_matrix,
                                id.vars = 1,
                                measure.vars = 2:ncol(rv$expl_matrix),
                                variable.name = "Samples",
                                value.name = "Count")
    
  })
  
  ## generate selected info matrix
  observeEvent(input$rna_expl_groupby_choices,{
    rv$info_expl <- rv$info[ , c("Samples",input$rna_expl_groupby_choices)]
  })
  
  
  ## generate result plots
  observeEvent(input$generate_result_expl, {
    # count hist
    rv$count_hist <- ggplot(data = rv$expl_matrix_long,
                            aes(x = Count)) + 
      xlab("") +
      ylab(expression(log[2](count + 1))) +
      geom_histogram(colour = "white", fill = "#525252") + 
      facet_wrap( ~ Samples, ncol = 3)
    output$display_count_hist <- renderPlot({
      print(rv$count_hist)
    })
    
    # box plot
    if(is.null(input$rna_expl_groupby_choices)){
      rv$boxplot <- ggplot(data = rv$expl_matrix_long,
                           aes(x = Samples,
                               y = Count)) +
        geom_boxplot() +
        xlab("") +
        ylab(expression(log[2](count + 1)))
    }else{
      df_boxplot <- merge(rv$expl_matrix_long, rv$info_expl, by = "Samples")    
      rv$boxplot <- ggplot(data = df_boxplot,
                           aes(x = Samples,
                               y = Count,
                               group = input$rna_expl_groupby_choices,
                               fill = input$rna_expl_groupby_choices)) +
        geom_boxplot(aes(fill = df_boxplot[ ,ncol(df_boxplot)])) +
        xlab("") +
        ylab(expression(log[2](count + 1)))
        
    }
    output$display_boxplot <- renderPlot({
      print(rv$boxplot)
    })
    
    # heatmap s2s
    sample_dist <- dist(t(rv$expl_matrix[,-1]))
    rv$sample_dist_matrix <- as.matrix(sample_dist)
    # colnames(rv$sample_dist_matrix) <- NULL
    rv$heatmap <- pheatmap(rv$sample_dist_matrix,
                           col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255))
    
    file_name <- paste(input$project_name,
                       "_expl_heatmap_s2s_dist",
                       ".tiff",
                       sep = "")
    file_path <- paste(getwd(),
                       "/tmp/",
                       file_name,
                       sep = "")
    tiff(filename = file_path,
         height = 6,
         width = 6,
         units = 'in',
         res = 300,
         compression = "lzw+p")
    print(rv$heatmap)
    graphics.off()
    
    output$display_heatmap_expl <- renderImage({
      list(src = file_path,
           width = 400,
           height = 400)
    }, deleteFile = FALSE)
    
    
    # pca
    tmp <- as.data.frame(t(rv$expl_matrix[, -1]))
    df_pca <- prcomp(tmp)
    df_out <- as.data.frame(df_pca$x)
    df_out$Samples <- input$rna_expl_sample_choices
    
    if(is.null(input$rna_expl_groupby_choices)){
      rv$pca <- ggplot(data = df_out,
                       aes(x = PC1,
                           y = PC2)) +
        geom_point()
    }else{
      
      df_out <- merge(df_out, rv$info_expl, by = "Samples")
      
      percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
      percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
      
      rv$pca <- ggplot(data = df_out,
                       aes(x = PC1,
                           y = PC2)) +
        geom_point(data = df_out,
                   # depend on "group-by"
                   aes(fill = df_out[, ncol(df_out)]),
                   shape = 21,
                   size = 3,
                   alpha = 0.5) +
        xlab(percentage[2]) + 
        ylab(percentage[3]) +
        geom_text(data = df_out,
                  aes(x = PC1 + 20,
                      y = PC2,
                      label = df_out[, ncol(df_out)]),
                  size = 2,
                  hjust = 0.5) +
        labs(fill = input$rna_expl_groupby_choices)
      
    }
    
    output$display_pca_expl <- renderPlot({
      print(rv$pca)
    })
    
  })
  
  ## generate and display count histgram
  # observeEvent(input$generate_count_hist,{
  #   rv$count_hist <- ggplot(data = rv$expl_matrix_long,
  #                           aes(x = Count)) + 
  #     xlab("") +
  #     ylab(expression(log[2](count + 1))) +
  #     geom_histogram(colour = "white", fill = "#525252") + 
  #     facet_wrap( ~ Samples, ncol = 3)
  #   output$display_count_hist <- renderPlot({
  #     print(rv$count_hist)
  #   })
  # })
  
  ## generate and display boxplot
  # observeEvent(input$generate_boxplot,{
  #   
  #   if(is.null(input$rna_expl_groupby_choices)){
  #     rv$boxplot <- ggplot(data = rv$expl_matrix_long,
  #                          aes(x = Samples,
  #                              y = Count)) +
  #       geom_boxplot() +
  #       xlab("") +
  #       ylab(expression(log[2](count + 1)))
  #   }else{
  #     df_boxplot <- merge(rv$expl_matrix_long, rv$info_expl, by = "Samples")    
  #     rv$boxplot <- ggplot(data = df_boxplot,
  #                          aes(x = Samples,
  #                              y = Count)) +
  #       # depend on group-by
  #       geom_boxplot(aes(fill = df_boxplot[ ,ncol(df_boxplot)])) +
  #       xlab("") +
  #       ylab(expression(log[2](count + 1))) +
  #       labs(fill = input$rna_expl_groupby_choices)
  #     
  #   }
  #   output$display_boxplot <- renderPlot({
  #     print(rv$boxplot)
  #   })
  #   
  # })
  
  ## generate and display heatmap
  # observeEvent(input$generate_heatmap_s2s,{
  #   sample_dist <- dist(t(rv$expl_matrix[,-1]))
  #   rv$sample_dist_matrix <- as.matrix(sample_dist)
  #   colnames(rv$sample_dist_matrix) <- NULL
  #   rv$heatmap <- pheatmap(rv$sample_dist_matrix,
  #                          clustering_distance_rows=sample_dist,
  #                          clustering_distance_cols=sample_dist,
  #                          col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255))
  #   
  #   file_name <- paste(input$project_name,
  #                      "_expl_heatmap_s2s_dist",
  #                      ".tiff",
  #                      sep = "")
  #   file_path <- paste(getwd(),
  #                      "/tmp/",
  #                      file_name,
  #                      sep = "")
  #   tiff(filename = file_path,
  #        height = 6,
  #        width = 6,
  #        units = 'in',
  #        res = 300,
  #        compression = "lzw+p")
  #   print(rv$heatmap)
  #   graphics.off()
  #   
  #   output$display_heatmap_expl <- renderImage({
  #     list(src = file_path,
  #          width = 400,
  #          height = 400)
  #   }, deleteFile = FALSE)
  #   
  #   
  # })
  
  ## PCA
  # observeEvent(input$generate_pca,{
  #   tmp <- as.data.frame(t(rv$expl_matrix[, -1]))
  #   df_pca <- prcomp(tmp)
  #   df_out <- as.data.frame(df_pca$x)
  #   df_out$Samples <- input$rna_expl_sample_choices
  #   
  #   if(is.null(input$rna_expl_groupby_choices)){
  #     rv$pca <- ggplot(data = df_out,
  #                      aes(x = PC1,
  #                          y = PC2)) +
  #       geom_point()
  #   }else{
  #     
  #     df_out <- merge(df_out, rv$info_expl, by = "Samples")
  #     
  #     percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
  #     percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
  #     
  #     rv$pca <- ggplot(data = df_out,
  #                      aes(x = PC1,
  #                          y = PC2)) +
  #       geom_point(data = df_out,
  #                  # depend on "group-by"
  #                  aes(fill = df_out[, ncol(df_out)]),
  #                  shape = 21,
  #                  size = 3,
  #                  alpha = 0.5) +
  #       xlab(percentage[2]) + 
  #       ylab(percentage[3]) +
  #       geom_text(data = df_out,
  #                 aes(x = PC1 + 20,
  #                     y = PC2,
  #                     label = df_out[, ncol(df_out)]),
  #                 size = 2,
  #                 hjust = 0.5) +
  #       labs(fill = input$rna_expl_groupby_choices)
  #     
  #   }
  #   
  #   output$display_pca_expl <- renderPlot({
  #     print(rv$pca)
  #   })
  # })
  
  ## export all generated plots
  observeEvent(input$export_plots_expl,{
    file_name <- paste(input$project_name,
                       "_expl_count_hist",
                       ".tiff",
                       sep = "")
    file_path <- paste(getwd(),
                       "/tmp/",
                       file_name,
                       sep = "")
    tiff(filename = file_path,
         height = 6,
         width = 6,
         units = 'in',
         res = 300,
         compression = "lzw+p")
    print(rv$count_hist)
    graphics.off()
    
    file_name1 <- paste(input$project_name,
                       "_expl_boxplot",
                       ".tiff",
                       sep = "")
    file_path1 <- paste(getwd(),
                       "/tmp/",
                       file_name1,
                       sep = "")
    tiff(filename = file_path1,
         height = 6,
         width = 6,
         units = 'in',
         res = 300,
         compression = "lzw+p")
    print(rv$boxplot)
    graphics.off()
    
    file_name2 <- paste(input$project_name,
                        "_expl_heatmap_s2s_dist",
                        ".tiff",
                        sep = "")
    file_path2 <- paste(getwd(),
                        "/tmp/",
                        file_name2,
                        sep = "")
    
    file_name3 <- paste(input$project_name,
                        "_expl_pca",
                        ".tiff",
                        sep = "")
    file_path3 <- paste(getwd(),
                        "/tmp/",
                        file_name3,
                        sep = "")
    tiff(filename = file_path3,
         height = 6,
         width = 6,
         units = 'in',
         res = 300,
         compression = "lzw+p")
    print(rv$pca)
    graphics.off()
    
    
    
    
    output$export_dir_expl <- renderPrint({
      cat(
        paste("plots export to:",
              file_path,
              file_path1,
              file_path2,
              file_path3,
              sep = "\n"))
    })

  })
  
  
  ## -------------  DEGseq  ----------------
  
  ## create and summary count table specified by user
  observeEvent(input$generate_count_table,{
    if(length(input$dge_choices) == 3){
      req(rv$ct)
      # create count table
      rv$ct1 <- droplevels(rv$ct[ , c("gene",
                                      input$dge_choices)])
      # summary count table
      output$count_table_txtbox <- renderPrint({
        summary(rv$ct1[-1])})
    }
  })
  
  
  ## trim count table and display result
  observeEvent(input$trim,{
    if(!is.null(rv$ct1)){
      # trim count table
      rv$ct2 <- droplevels(subset(rv$ct1,
                                  rowSums(rv$ct1[, -1]) > as.numeric(input$row_sum)))
      # display trimmed result
      output$trim_txtbox <- renderPrint({
        cat(
          paste(nrow(rv$ct2),
                "genes left, down from",
                nrow(rv$ct1),
                "genes",
                sep = " "))})
      output$trimmed_count_table_txtbox <- renderPrint({
        summary(rv$ct2)
      })
    }
  })
  
  
  ## run DEGexp, read in result table, add column mu base+1
  observeEvent(input$run_DEGexp,{
    
    ## clear old results
    output$result_table1 <- DT::renderDataTable({})
    output$result_table2 <- DT::renderDataTable({})
    # export dir of result table
    output$export_dir1 <- renderPrint({})
    output$export_dir2 <- renderPrint({})
    # export dir of ma plot
    output$export_maplot_dir <- renderPrint({})
    # ma plot
    output$result_plot1 <- renderPlot({})
    output$result_plot2 <- renderPlot({})
    output$result_kable1 <- renderPrint({})
    output$result_kable2 <- renderPrint({})
    # result of change in gene exp
    output$display_venn_diagram1 <- renderPlot({})
    output$display_venn_diagram2 <- renderPlot({})
    output$display_up.dn_dn.up_table = DT::renderDataTable({})
    output$display_up.dn_dn.up_heatmap <- renderPlot({})
    output$display_change_in_gene_exp_dir <- renderPrint({})
    output$up.dn_number_txtbox <- renderPrint({})
    
    
    withProgress(message = "Running DEGexp",
                 value = 0,
                 expr = {
                   incProgress(0.2, detail = "Processing trt2-trt1")
                   DEGexp(geneExpMatrix1 = rv$ct2,
                          geneCol1 = 1, 
                          expCol1 = 3, 
                          groupLabel1 = colnames(rv$ct2)[3],
                          
                          geneExpMatrix2 = rv$ct2,
                          geneCol2 = 1, 
                          expCol2 = 2,
                          groupLabel2 = colnames(rv$ct2)[2],
                          
                          # depend on fold_change
                          foldChange = as.numeric(input$fold_change),
                          # depend on q_value
                          qValue = as.numeric(input$q_value),
                          thresholdKind = 5, 
                          rawCount = TRUE,
                          normalMethod = "none",
                          method = "MARS",
                          outputDir = "tmp")
                   rv$table_trt2_trt1 <- fread(paste(getwd(),"/tmp/output_score.txt", sep = ""))

                   incProgress(0.3, detail = "Processing trt3-trt2")
                   DEGexp(geneExpMatrix1 = rv$ct2,
                          geneCol1 = 1,
                          expCol1 = 4,
                          groupLabel1 = colnames(rv$ct2)[4],
                          
                          geneExpMatrix2 = rv$ct2,
                          geneCol2 = 1,
                          expCol2 = 3,
                          groupLabel2 = colnames(rv$ct2)[3],
                          foldChange = as.numeric(input$fold_change),
                          qValue = as.numeric(input$q_value),
                          thresholdKind = 5,
                          rawCount = TRUE,
                          normalMethod = "none",
                          method = "MARS",
                          outputDir = "tmp")
                   rv$table_trt3_trt2 <- fread(paste(getwd(),"/tmp/output_score.txt", sep = ""))
                   
                   incProgress(0.2, detail = "Displaying Result Tables")
                   # table_trt2_trt1
                   rv$table_trt2_trt1[ , mu := ifelse(value1 != 0 & value2 != 0, 
                                                      (log2(value1) + log2(value2))/2,
                                                      ifelse(value1 != 0 & value2 == 0, 
                                                             (log2(value1) + log2(value2+1))/2,
                                                             ifelse(value1 == 0 & value2 != 0, 
                                                                    (log2(value1 + 1) + log2(value2))/2,
                                                                    (log2(value1 + 1) + log2(value2 + 1))/2)))
                                       ]
                   rv$table_trt2_trt1$clr <- "black"
                   rv$table_trt2_trt1$clr[rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$q_value)] <- "purple"
                   rv$table_trt2_trt1$clr[rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$q_value) 
                                          & abs(rv$table_trt2_trt1$`log2(Fold_change)`) >= as.numeric(input$fold_change)] <- "red"
                   rv$table_trt2_trt1$pch <- 46
                   rv$table_trt2_trt1$pch[rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$q_value)] <- 3
                   rv$table_trt2_trt1$pch[rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$q_value) 
                                          & abs(rv$table_trt2_trt1$`log2(Fold_change)`) >= as.numeric(input$fold_change)] <- 4
                   rv$table_trt2_trt1$pch <- factor(rv$table_trt2_trt1$pch,
                                                    levels = c(46, 3, 4))
                   
                   rv$kable1 <- kable(table(rv$table_trt2_trt1$clr,
                                            rv$table_trt2_trt1$pch))
                   
                   rv$table_trt2_trt1_round <- rv$table_trt2_trt1
                   rv$table_trt2_trt1_round[, c(4:9,11)] <- round(rv$table_trt2_trt1_round[, c(4:9,11)], 2)
                   
                   output$result_table1 = DT::renderDataTable({
                     DT::datatable(rv$table_trt2_trt1_round,
                                   filter = "top",
                                   options = list(scrollX = TRUE,
                                                  columnDefs = list(list(searchable = FALSE, targets = c(1:4,6:7,11:13))))) # %>%
                       # formatRound(c("log2(Fold_change)",
                       #               "log2(Fold_change) normalized", 
                       #               "z-score", 
                       #               "p-value", 
                       #               "q-value(Benjamini et al. 1995)", 
                       #               "q-value(Storey et al. 2003)",
                       #               "Signature(q-value(Storey et al. 2003) < 0.01)",
                       #               "mu"), 2)
                     })
                   
                   # table_trt3_trt2
                   rv$table_trt3_trt2[ , mu := ifelse(value1 != 0 & value2 != 0, 
                                                      (log2(value1) + log2(value2))/2,
                                                      ifelse(value1 != 0 & value2 == 0, 
                                                             (log2(value1) + log2(value2+1))/2,
                                                             ifelse(value1 == 0 & value2 != 0, 
                                                                    (log2(value1 + 1) + log2(value2))/2,
                                                                    (log2(value1 + 1) + log2(value2 + 1))/2)))
                                       ]   
                   rv$table_trt3_trt2$clr <- "black"
                   rv$table_trt3_trt2$clr[rv$table_trt3_trt2$`q-value(Storey et al. 2003)` < as.numeric(input$q_value)] <- "purple"
                   rv$table_trt3_trt2$clr[rv$table_trt3_trt2$`q-value(Storey et al. 2003)` < as.numeric(input$q_value) 
                                          & abs(rv$table_trt3_trt2$`log2(Fold_change)`) >= as.numeric(input$fold_change)] <- "red"
                   rv$table_trt3_trt2$pch <- 46
                   rv$table_trt3_trt2$pch[rv$table_trt3_trt2$`q-value(Storey et al. 2003)` < as.numeric(input$q_value)] <- 3
                   rv$table_trt3_trt2$pch[rv$table_trt3_trt2$`q-value(Storey et al. 2003)` < as.numeric(input$q_value) 
                                          & abs(rv$table_trt3_trt2$`log2(Fold_change)`) >= as.numeric(input$fold_change)] <- 4
                   rv$table_trt3_trt2$pch <- factor(rv$table_trt3_trt2$pch,
                                                    levels = c(46, 3, 4))
                   
                   rv$kable2 <- kable(table(rv$table_trt3_trt2$clr,
                                            rv$table_trt3_trt2$pch))
                   
                   rv$table_trt3_trt2_round <- rv$table_trt3_trt2
                   rv$table_trt3_trt2_round[, c(4:9,11)] <- round(rv$table_trt3_trt2_round[, c(4:9,11)], 2)
                   
                   output$result_table2 = DT::renderDataTable({
                     DT::datatable(rv$table_trt3_trt2_round,
                                   filter = "top",
                                   options = list(scrollX = TRUE,
                                                  columnDefs = list(list(searchable = FALSE, targets = c(1:4,6:7,11:13))))) # %>%
                       # formatRound(c("log2(Fold_change)",
                       #               "log2(Fold_change) normalized", 
                       #               "z-score", 
                       #               "p-value", 
                       #               "q-value(Benjamini et al. 1995)", 
                       #               "q-value(Storey et al. 2003)",
                       #               "Signature(q-value(Storey et al. 2003) < 0.01)",
                       #               "mu"), 2)
                     })
                   Sys.sleep(3)
                   incProgress(0.3, detail = "Processing 100%")
                   Sys.sleep(1)
                 })
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  })
  
  
  ## export result table trt2_trt1
  observeEvent(input$export_table_trt2_trt1,{
    file_name <- paste(input$project_name, 
                       "_RNAseq_DEGseq_", 
                       colnames(rv$ct1)[3], 
                       "-", 
                       colnames(rv$ct1)[2],
                       ".csv", 
                       sep="")
    
    file_path <- paste(getwd(),
                       "/tmp/",
                       file_name,
                       sep= "")
    # display export dir
    output$export_dir1 <- renderPrint({
      cat(paste("table export to:", file_path))
    })
    # wirte CSV file
    observe({
      write.csv(rv$table_trt2_trt1,
                file = file_path,
                row.names = FALSE)})
  })
  
  
  ## export result table trt3_trt2
  observeEvent(input$export_table_trt3_trt2,{
    file_name <- paste(input$project_name, 
                       "_RNAseq_DEGseq_", 
                       colnames(rv$ct1)[4], 
                       "-", 
                       colnames(rv$ct1)[3],
                       ".csv", 
                       sep="")
    
    file_path <- paste(getwd(),
                       "/tmp/",
                       file_name,
                       sep= "")
    # display export dir
    output$export_dir2 <- renderPrint({
      cat(paste("table export to:", file_path))
    })
    # wirte CSV file
    observe({
      write.csv(rv$table_trt3_trt2,
                file = file_path,
                row.names = FALSE)})
  })
  
  
  ## generate and display MA plot
  observeEvent(input$generate_maplot,{
    
    # create MA plot1
    rv$p1 <- ggplot(rv$table_trt2_trt1,
                    aes(x = mu,
                        y = `log2(Fold_change) normalized`,
                        colour = clr,
                        shape = pch)) +
      scale_shape_manual(name = "Legend:",
                         labels = c("No significance",
                                    paste("p-Value < ", 
                                          as.numeric(input$q_value)),
                                    paste("p-Value < ", 
                                          as.numeric(input$q_value), 
                                          " & abs(log2) >= ", 
                                          as.numeric(input$fold_change))),
                         values = c(46, 3, 4)) +
      scale_color_manual(name = "Legend:",
                         values = c("black",
                                    "purple",
                                    "red"),
                         labels = c("No significance",
                                    paste("p-Value < ", 
                                          as.numeric(input$q_value)),
                                    paste("p-Value < ", 
                                          as.numeric(input$q_value), 
                                          " & abs(log2) >= ", 
                                          as.numeric(input$fold_change)))) +
      geom_hline(yintercept = c(-as.numeric(input$fold_change), as.numeric(input$fold_change)),
                 lty = 2) +
      ggtitle(input$ma_title1)  +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            legend.position = "top") + 
      geom_point() +
      labs(x = "log2Mean", y = "log2FoldChange")
    
    # create MA plot2 
    rv$p2 <- ggplot(rv$table_trt3_trt2,
                    aes(x = mu,
                        y = `log2(Fold_change)`,
                        colour = clr,
                        shape = pch)) +
      scale_shape_manual(name = "Legend:",
                         labels = c("No significance",
                                    paste("p-Value < ", 
                                          as.numeric(input$q_value)),
                                    paste("p-Value < ", 
                                          as.numeric(input$q_value), 
                                          " & abs(log2) >= ", 
                                          as.numeric(input$fold_change))),
                         values = c(46, 3, 4)) +
      scale_color_manual(name = "Legend:",
                         values = c("black",
                                    "purple",
                                    "red"),
                         labels = c("No significance",
                                    # depend on q_value and fold_change
                                    paste("p-Value < ", 
                                          as.numeric(input$q_value)),
                                    paste("p-Value < ", 
                                          as.numeric(input$q_value), 
                                          " & abs(log2) >= ", 
                                          as.numeric(input$fold_change)))) +
      geom_hline(yintercept = c(-as.numeric(input$fold_change), as.numeric(input$fold_change)),
                 lty = 2) +
      # depend on ma_title2
      ggtitle(input$ma_title2)  +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            legend.position = "top") + 
      geom_point() +
      labs(x = "log2Mean", y = "log2FoldChange")
    
    # display MA plot1
    output$result_plot1 <- renderPlot({
      print(rv$p1)
    })
    
    # display kable1
    output$result_kable1 <- renderPrint({
      rv$kable1
    })
    
    # display MA plot2
    output$result_plot2 <- renderPlot({
      print(rv$p2)
    })
    
    # display kable2
    output$result_kable2 <- renderPrint({
      rv$kable2
    })
  })
  
  
  ## export MA plots
  observeEvent(input$export_maplot,{
    file_name1 <- paste(input$project_name,
                        "_MA_plot_",
                        colnames(rv$ct1)[3], 
                        "-", 
                        colnames(rv$ct1)[2],
                        ".tiff",
                        sep = "")
    file_path1 <- paste(getwd(),
                        "/tmp/",
                        file_name1,
                        sep = "")
    file_name2 <- paste(input$project_name,
                        "_MA_plot_",
                        colnames(rv$ct1)[4], 
                        "-", 
                        colnames(rv$ct1)[3],
                        ".tiff",
                        sep = "")
    file_path2 <- paste(getwd(),
                        "/tmp/",
                        file_name2,
                        sep = "")
    # export plot1
    tiff(filename = file_path1,
         height = 6,
         width = 6,
         units = 'in',
         res = 300,
         compression = "lzw+p")
    print(rv$p1)
    
    # export plot2
    tiff(filename = file_path2,
         height = 6,
         width = 6,
         units = 'in',
         res = 300,
         compression = "lzw+p")
    print(rv$p2)
    graphics.off()
    
    # display export dir
    output$export_maplot_dir <- renderPrint({
      cat(paste("plot export to:", file_path1, file_path2, sep = "\n"))
    })
  })
  
  
  
  
  ## generate and display venn diagram, sign_gene_table, and two circle heatmap
  observeEvent(input$generate_up.dn_dn.up_table_and_heatmap,{
    withProgress(message = "Generating results",
                 value = 0,
                 expr = {
                   incProgress(0.2, detail = "Generating venn diagram")
                   # display up.dn dn.up number
                   rv$trt2_trt1_up <- rv$table_trt2_trt1[rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$q_value)
                                                         & rv$table_trt2_trt1$`log2(Fold_change)` >= as.numeric(input$fold_change), ]$GeneNames
                   rv$trt2_trt1_dn <- rv$table_trt2_trt1[rv$table_trt2_trt1$`q-value(Storey et al. 2003)` < as.numeric(input$q_value) 
                                                         & rv$table_trt2_trt1$`log2(Fold_change)` <= -as.numeric(input$fold_change), ]$GeneNames
                   rv$trt3_trt2_up <- rv$table_trt3_trt2[rv$table_trt3_trt2$`q-value(Storey et al. 2003)` < as.numeric(input$q_value) 
                                                         & rv$table_trt3_trt2$`log2(Fold_change)` >= as.numeric(input$fold_change), ]$GeneNames
                   rv$trt3_trt2_dn <- rv$table_trt3_trt2[rv$table_trt3_trt2$`q-value(Storey et al. 2003)` < as.numeric(input$q_value) 
                                                         & rv$table_trt3_trt2$`log2(Fold_change)` <= -as.numeric(input$fold_change), ]$GeneNames
                   
                   rv$up.dn <- rv$trt2_trt1_up[rv$trt2_trt1_up %in% rv$trt3_trt2_dn]
                   rv$dn.up <- rv$trt2_trt1_dn[rv$trt2_trt1_dn %in% rv$trt3_trt2_up]
                   
                   rv$up.dn_dn.up_genes <- unique(c(rv$up.dn,
                                                    rv$dn.up))
                   
                   rv$merge1  <- rv$table_trt2_trt1[rv$table_trt2_trt1$GeneNames %in% rv$up.dn_dn.up_genes, c("GeneNames","log2(Fold_change)")]
                   rv$merge2  <- rv$table_trt3_trt2[rv$table_trt3_trt2$GeneNames %in% rv$up.dn_dn.up_genes, c("GeneNames","log2(Fold_change)")]
                   
                   rv$up.dn_dn.up_table <- merge(rv$merge1, rv$merge2, by = "GeneNames")
                   colnames(rv$up.dn_dn.up_table) <- c("gene",paste(colnames(rv$ct2)[3],"-",colnames(rv$ct2)[2]), paste(colnames(rv$ct2)[4],"-",colnames(rv$ct2)[3]))
                   rv$up.dn_dn.up_table <- rv$up.dn_dn.up_table[order(rv$up.dn_dn.up_table[ ,2], decreasing = TRUE), ]
                   output$up.dn_number_txtbox <- renderPrint({
                     cat(paste0("trt2_trt1_up:", 
                                length(rv$trt2_trt1_up), 
                                " trt2_trt1_dn:", 
                                length(rv$trt2_trt1_dn), 
                                " trt3_trt2_up:", 
                                length(rv$trt3_trt2_up), 
                                " trt3_trt2_dn:", 
                                length(rv$trt3_trt2_dn),
                                " up.dn:",
                                length(rv$up.dn),
                                " dn.up:",
                                length(rv$dn.up)))
                   })
                   
                   # display venn diagram
                   output$display_venn_diagram1 <- renderPlot({
                     rv$venn_diagram1 <- draw.pairwise.venn(area1 = length(rv$trt2_trt1_up),
                                                            area2 = length(rv$trt3_trt2_dn),
                                                            cross.area = length(rv$up.dn),
                                                            scaled = TRUE,
                                                            col = c("green3", "firebrick"))
                   })
                   
                   output$display_venn_diagram2 <- renderPlot({
                     rv$venn_diagram2 <- draw.pairwise.venn(area1 = length(rv$trt2_trt1_dn),
                                                            area2 = length(rv$trt3_trt2_up),
                                                            cross.area = length(rv$dn.up),
                                                            scaled = TRUE,
                                                            col = c("firebrick", "green3"))
                   })
                   
                   incProgress(0.2, detail = "Generating result table")
                   # display up.dn_dn.up_table
                   output$display_up.dn_dn.up_table = DT::renderDataTable({
                     DT::datatable(rv$up.dn_dn.up_table,
                                   options = list(scrollX = TRUE))})
                   
                   
                   incProgress(0.2, detail = "Generating two circle heatmap")
                   # display two circle heatmap
                   # create two circle heatmap non-normalized foldchange
                   rv$up.dn_dn.up_table_long <- melt(data = rv$up.dn_dn.up_table,
                                                     id.vars = 1,
                                                     measure.vars = 2:3,
                                                     variable.name = "Comparison",
                                                     value.name = "Gene Expression Diff")
                   rv$up.dn_dn.up_table_long$Comparison <- factor(rv$up.dn_dn.up_table_long$Comparison,
                                                                  levels = c(paste(colnames(rv$ct2)[4],"-",colnames(rv$ct2)[3]),
                                                                             paste(colnames(rv$ct2)[3],"-",colnames(rv$ct2)[2])))
                   rv$lvls <- rv$up.dn_dn.up_table_long[rv$up.dn_dn.up_table_long$Comparison == paste(colnames(rv$ct2)[3],"-",colnames(rv$ct2)[2]), ]
                   rv$up.dn_dn.up_table_long$gene <- factor(rv$up.dn_dn.up_table_long$gene,
                                                            levels = rv$lvls$gene[order(rv$lvls$`Gene Expression Diff`)])
                   
                   rv$p3 <- ggplot(data = rv$up.dn_dn.up_table_long) +                        
                     coord_polar("y",
                                 start = 0,
                                 direction = -1) +
                     geom_tile(aes(x = as.numeric(Comparison),     
                                   y = gene,                       
                                   fill = `Gene Expression Diff`), 
                               color = "white") +                   
                     geom_text(data = rv$up.dn_dn.up_table_long[Comparison == paste(colnames(rv$ct2)[3],"-",colnames(rv$ct2)[2]), ],  
                               aes(x = rep(1.75,
                                           nlevels(gene)),
                                   y = gene,
                                   label = unique(gene),
                                   angle = 90 + seq(from = 0,
                                                    to = 360,
                                                    length.out = nlevels(gene))[as.numeric(gene)]),
                               hjust = 0) +
                     scale_fill_gradient2(low = "red",             
                                          high = "green",
                                          mid = "grey",
                                          midpoint = 0,            
                                          name = "Gene Expr Diff") + 
                     scale_x_continuous(limits = c(0,
                                                   max(as.numeric(rv$up.dn_dn.up_table_long$Comparison)) + 0.5),
                                        expand = c(0, 0)) +
                     scale_y_discrete("",
                                      expand = c(0, 0)) +
                     ggtitle(paste("Changes in Gene Expression", 
                                   "q-value >", 
                                   input$q_value, 
                                   "log2foldchange >=", 
                                   input$fold_change)) +    
                     theme(plot.title = element_text(hjust = 0.5),
                           axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.y = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks.y = element_blank())
                   output$display_up.dn_dn.up_heatmap <- renderPlot({
                     print(rv$p3)
                   })
                   Sys.sleep(2)
                   incProgress(0.4, detail = "Processing 100%")
                   Sys.sleep(1)
                 })
    
  })
  
  
  ## export up.dn_dn.up_table two-circle heatmap venn diagrams
  observeEvent(input$export_change_in_gene_exp,{
    
    # export up.dn_dn.up_table
    file_name <- paste(input$project_name, 
                       "_RNAseq_DEGseq_", 
                       colnames(rv$ct1)[4], 
                       "-", 
                       colnames(rv$ct1)[3],
                       "-", 
                       colnames(rv$ct1)[2],
                       "_sign_gene.csv", 
                       sep="")
    
    file_path <- paste(getwd(),
                       "/tmp/",
                       file_name,
                       sep= "")
    observe({
      write.csv(rv$up.dn_dn.up_table,
                file = file_path,
                row.names = FALSE)})
    
    # export heatmap
    file_name1 <- paste(input$project_name,
                        "_change_in_gene_exp_",
                        colnames(rv$ct1)[4],
                        "-",
                        colnames(rv$ct1)[3], 
                        "-", 
                        colnames(rv$ct1)[2],
                        "_heatmap.tiff",
                        sep = "")
    
    file_path1 <- paste(getwd(),
                        "/tmp/",
                        file_name1,
                        sep = "")
    
    tiff(filename = file_path1,
         height = 8,
         width = 8,
         units = 'in',
         res = 300,
         compression = "lzw+p")
    print(rv$p3)
    graphics.off()
    
    # export venn diagrams
    file_name2 <- paste(input$project_name,
                        "_change_in_gene_exp_",
                        colnames(rv$ct1)[4],
                        "-",
                        colnames(rv$ct1)[3], 
                        "-", 
                        colnames(rv$ct1)[2],
                        "_venn1.tiff",
                        sep = "")
    file_path2 <- paste(getwd(),
                        "/tmp/",
                        file_name2,
                        sep = "")
    tiff(filename = file_path2,
         height = 4,
         width = 5,
         units = 'in',
         res = 300,
         compression = "lzw+p")
    grid.draw(rv$venn_diagram1)
    graphics.off()
    
    
    file_name3 <- paste(input$project_name,
                        "_change_in_gene_exp_",
                        colnames(rv$ct1)[4],
                        "-",
                        colnames(rv$ct1)[3], 
                        "-", 
                        colnames(rv$ct1)[2],
                        "_venn2.tiff",
                        sep = "")
    file_path3 <- paste(getwd(),
                        "/tmp/",
                        file_name3,
                        sep = "")
    tiff(filename = file_path3,
         height = 4,
         width = 5,
         units = 'in',
         res = 300,
         compression = "lzw+p")
    grid.draw(rv$venn_diagram2)
    graphics.off()
    
    
    # display export dir
    output$display_change_in_gene_exp_dir <- renderPrint({
      cat(paste("table export to:", 
                file_path, 
                "plots export to:", 
                file_path1, 
                file_path2,
                file_path3,
                sep = "\n"))
    })
    
  })
  
  
  ## ---------------- DEseq2 -----------------------
  observeEvent(input$generate_count_matrix,{
    dt <- rv$ct[, 7:ncol(rv$ct)]
    rv$countdata <- as.matrix(dt)
    rownames(rv$countdata) <- rv$ct$gene
    output$display_summary_of_count_matrix <- renderPrint({
      summary(rv$countdata)
    })
  })
  
  observeEvent(input$generate_design_matrix,{
    output$display_summary_of_design_matrix <- renderPrint({
      summary(rv$info)
    })
  })
  
  
  
  observeEvent(input$generate_deseqdt_from_matrix,{
    # single factor
    if(length(input$deseq2_design_formula) == 1 &
       input$multi_factor_option == "none"){
      design <- as.formula(paste("~",
                                 input$deseq2_design_formula))
      # multi factor as_group
    }else if(length(input$deseq2_design_formula) > 1 &
             input$multi_factor_option == "as_group"){
      rv$info$group <- factor(apply(rv$info[ , input$deseq2_design_formula],
                                    1 , 
                                    paste , 
                                    collapse = "_" ))
      design <- as.formula("~ group")
      # multi factor without interaction
    }else if(length(input$deseq2_design_formula) > 1 &
             input$multi_factor_option == "o_interaction"){
      design <- as.formula(paste("~",
                                 paste(input$deseq2_design_formula, 
                                       collapse = "+")))
      # multi factor with interaction
    } else if(length(input$deseq2_design_formula) > 1 &
              input$multi_factor_option == "w_interaction"){
      design <- as.formula(paste("~",
                                 paste(input$deseq2_design_formula, 
                                       collapse = "+"),
                                 
                                 paste("+", 
                                       paste(input$deseq2_design_formula,
                                             collapse = ":"))))
    }
    
    if(class(design) == "formula"){
      output$display_formula <- renderPrint({
        cat(
          paste("The design formula:",
                paste(as.character(design),
                      collapse = " "))
        )
      })
      rv$dds <- DESeqDataSetFromMatrix(countData = rv$countdata,
                                       colData = rv$info,
                                       design = design)
      output$display_deseqdt <- renderPrint({
        rv$dds
      })
    }else{
      output$display_formula <- renderPrint({
        "Formula is not correct"
      })
      output$display_deseqdt <- renderPrint({})
    }
  })
  
  
  observeEvent(input$trim_deseqdt,{
    rv$dds_trimmed <- rv$dds[rowSums(counts(rv$dds)) >= as.numeric(input$row_sum_deseq2), ]
    output$display_trim_deseqdt_result <- renderPrint({
      cat(
        paste(nrow(counts(rv$dds_trimmed)),
              "genes left, down from",
              nrow(counts(rv$dds)),
              "genes")
      )
    })
    output$display_trimmed_deseqdt <- renderPrint({
      rv$dds_trimmed
    })
    
  })
  
  ## run DESeq
  observeEvent(input$run_deseq,{
    snowparam <- SnowParam(workers = snowWorkers(), 
                           type = "SOCK")
    register(snowparam, 
             default = TRUE)
    
    withProgress(message = "Running DESeq", 
                 value = 0,
                 expr = {
                   incProgress(0.2, detail = "Processing 20%")
                   incProgress(0.3, detail = "Processing 50%")
                   incProgress(0.2, detail = "Processing 70%")
                   rv$dds_res <- DESeq(rv$dds_trimmed,
                                       fitType = "local",
                                       parallel = TRUE)
                   incProgress(0.3, detail = "Processing 100%")
                   Sys.sleep(1)
                 })
    # update contrast choice
    updateSelectInput(session,
                      inputId = "select_contrast_add",
                      label = "Select names of the log2 fold changes to add:",
                      choices = resultsNames(rv$dds_res),
                      selected = NULL)
    
    updateSelectInput(session,
                      inputId = "select_contrast_substract",
                      label = "Select names of the log2 fold changes to substract (optional):",
                      choices = resultsNames(rv$dds_res),
                      selected = NULL)  
  })
  
  
  
  ## generate result table and plot
  observeEvent(input$generate_deseq_results,{
    withProgress(message = "Extracting results with contrast", 
                 value = 0,
                 expr = {
                   incProgress(0.2, detail = "Preparing result table")
                   if(is.null(input$select_contrast_substract)){
                     contrast_list <- list(input$select_contrast_add)
                   }else{
                     contrast_list <- list(input$select_contrast_add,
                                           input$select_contrast_substract)
                   }
                   rv$dds_res_contrast <- results(rv$dds_res,
                                                  contrast = contrast_list)
                   rv$dtf_res_contrast <- data.frame(gene = rownames(rv$dds_res_contrast),
                                                     do.call("cbind", 
                                                             rv$dds_res_contrast@listData))
                   rv$dtf_res_contrast$clr <- "black"
                   rv$dtf_res_contrast$clr[rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value)] <- "purple"
                   rv$dtf_res_contrast$clr[rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) & abs(rv$dtf_res_contrast$log2FoldChange) >= as.numeric(input$deseq2_fold_change)] <- "red"
                   rv$dtf_res_contrast$pch <- 46
                   rv$dtf_res_contrast$pch[rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value)] <- 3
                   rv$dtf_res_contrast$pch[rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) & abs(rv$dtf_res_contrast$log2FoldChange) >= as.numeric(input$deseq2_fold_change)] <- 4
                   rv$dtf_res_contrast$pch <- factor(rv$dtf_res_contrast$pch,
                                                     levels = c(46, 3, 4))
                   rv$dtf_res_contrast_round <- rv$dtf_res_contrast
                   rv$dtf_res_contrast_round[, c(2:7)] <-round(rv$dtf_res_contrast_round[, c(2:7)], 2)
                   
                   output$display_dtf_res_contrast <- DT::renderDataTable({
                     DT::datatable(rv$dtf_res_contrast_round,
                                   filter = "top",
                                   options = list(scrollX = TRUE,
                                                  columnDefs = list(list(searchable = FALSE, targets = c(1:2, 4:6, 8:9))))) #%>%
                       # formatRound(c("baseMean",
                       #               "log2FoldChange", 
                       #               "lfcSE", 
                       #               "stat", 
                       #               "pvalue", 
                       #               "padj"), 2)
                   })
                   
                   incProgress(0.2, detail = "Preparing MA plot")
                   rv$deseq2_ma <- ggplot(rv$dtf_res_contrast,
                                          aes(x = log(baseMean + 1),
                                              y = `log2FoldChange`,
                                              colour = clr,
                                              shape = pch)) +
                     scale_shape_manual(name = "Legend:",
                                        values = c(46, 3, 4),
                                        labels = c("No significance",
                                                   paste("p-Value < ", 
                                                         as.numeric(input$deseq2_p_value)),
                                                   paste("p-Value < ", 
                                                         as.numeric(input$deseq2_p_value), 
                                                         " & abs(log2) >= ", 
                                                         as.numeric(input$deseq2_fold_change)))) +
                     scale_color_manual(name = "Legend:",
                                        values = c("black",
                                                   "purple",
                                                   "red"),
                                        labels = c("No significance",
                                                   paste("p-Value < ", 
                                                         as.numeric(input$deseq2_p_value)),
                                                   paste("p-Value < ", 
                                                         as.numeric(input$deseq2_p_value), 
                                                         " & abs(log2) >= ", 
                                                         as.numeric(input$deseq2_fold_change)))) +
                     geom_hline(yintercept = c(-as.numeric(input$deseq2_fold_change), as.numeric(input$deseq2_fold_change)),
                                lty = 2) +
                     theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.background = element_blank(),
                           axis.line = element_line(colour = "black"),
                           plot.title = element_text(hjust = 0.5),
                           legend.position = "top") + 
                     geom_point() +
                     labs(x = "log2(BaseMean + 1)", y = "log2FoldChange")
                   
                   output$display_deseq2_ma <- renderPlot({
                     print(rv$deseq2_ma)
                   })
                   
                   # sign gene number
                   result <- c("up-regulated DEG",
                               "non DEG",
                               "down-regulated DEG")
                   number <- c(
                     sum(rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) & rv$dtf_res_contrast$log2FoldChange >= as.numeric(input$deseq2_fold_change), na.rm = TRUE),
                     nrow(rv$dtf_res_contrast) - sum(rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) & abs(rv$dtf_res_contrast$log2FoldChange) >= as.numeric(input$deseq2_fold_change), na.rm = TRUE),
                     sum(rv$dtf_res_contrast$padj < as.numeric(input$deseq2_p_value) & rv$dtf_res_contrast$log2FoldChange <= - as.numeric(input$deseq2_fold_change), na.rm = TRUE))
                   deseq2_sign_number_tb <- data.frame(result, number)
                   
                   output$deseq2_sign_number <- renderTable(deseq2_sign_number_tb)
                   
                   # estiated dispersion plot
                   incProgress(0.2, detail = "Preparing estimated dispersion plot")
                   output$deseq2_estimated_dispersion <- renderImage({
                     file_name <- paste(input$project_name,
                                        "_DESeq2_disp_est",
                                        ".tiff",
                                        sep = "")
                     file_path <- paste(getwd(),
                                        "/tmp/",
                                        file_name,
                                        sep = "")
                     tiff(filename = file_path,
                          height = 6,
                          width = 6,
                          units = 'in',
                          res = 300,
                          compression = "lzw+p")
                     plotDispEsts(rv$dds_res)
                     graphics.off()
                     
                     list(src = file_path,
                          width = 500,
                          height = 450)
                   }, deleteFile = FALSE)
                   
                   incProgress(0.4, detail = "Processing 100%")
                   Sys.sleep(2)
                 })
    # renew MA plot
    observeEvent(input$enter_deseq2_ma_title,{
      rv$deseq2_ma_w_title <- rv$deseq2_ma + 
        ggtitle(input$enter_deseq2_ma_title)
      output$display_deseq2_ma <- renderPlot({
        print(rv$deseq2_ma_w_title)
      })
    })
    
  })
  
  
  
  ## generate count plot
  observeEvent(input$generate_count_plot,{
    if(is.null(input$enter_gene_name) | is.null(input$intgroup)){
      showNotification("Please enter one gene name and select interested group(s)",
                       type = "error",
                       duration = 15)
    }else if(input$enter_gene_name %in% rownames(rv$dds_trimmed)  == FALSE ){
      showNotification("Please enter the exact correct gene name",
                       type = "error",
                       duration = 15)
    }else{
      dt <- plotCounts(rv$dds_res, gene = input$enter_gene_name, intgroup = input$intgroup, returnData = TRUE)
      if(length(input$intgroup) > 1){
        dt$group <- apply(dt[, -1], 1, paste, collapse= ":")
      }
      
      rv$count_plot <- ggplot(dt, aes(x=dt[,ncol(dt)], y=count)) +
        geom_point(position=position_jitter(w=0.1,h=0)) +
        labs(y = "normalized count", x = NULL)
      output$deseq2_count_plot  <- renderPlot({print(rv$count_plot)})
    }
  })
  
  # export result
  observeEvent(input$export_deseq_results,{
    file_name <- paste(input$project_name,
                       "_DESeq2_MA",
                       input$enter_deseq2_ma_title,
                       ".tiff",
                       sep = "")
    file_path <- paste(getwd(),
                       "/tmp/",
                       file_name,
                       sep = "")
    tiff(filename = file_path,
         height = 6,
         width = 6,
         units = 'in',
         res = 300,
         compression = "lzw+p")
    print(rv$deseq2_ma_w_title)
    graphics.off()
    
    file_name1 <- paste(input$project_name,
                        "_DESeq2_result_table.csv",
                        sep = "")
    file_path1 <- paste(getwd(),
                        "/tmp/",
                        file_name1,
                        sep= "")
    write.csv(rv$dtf_res_contrast,
              file = file_path1,
              row.names = FALSE)
    
    file_name2 <- paste(input$project_name,
                       "_DESeq2_DispEst.tiff",
                       sep = "")
    file_path2 <- paste(getwd(),
                       "/tmp/",
                       file_name2,
                       sep = "")
    
    if(!is.null(rv$count_plot)){
      file_name3 <- paste(input$project_name,
                          "_DESeq2_CountPlot.tiff",
                          sep = "")
      file_path3 <- paste(getwd(),
                          "/tmp/",
                          file_name3,
                          sep = "")
      tiff(filename = file_path3,
           height = 6,
           width = 6,
           units = 'in',
           res = 300,
           compression = "lzw+p")
      print(rv$count_plot)
      graphics.off()
    }
    
    
    
    # display export dir
    
    output$deseq2_results_dir <- renderPrint({
      cat(paste("result table export to:",
                file_path1, 
                "plots export to:", 
                file_path,
                file_path2,
                ifelse(!is.null(rv$count_plot),file_path3, ""),
                sep = "\n"))
    })
  })
  
  
  
  
  
  
  
  # end  
}

shinyApp(ui, server)



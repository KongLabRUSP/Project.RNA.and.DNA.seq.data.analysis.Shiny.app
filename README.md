##  Project: RNA and DNA sequencing data analysis Shiny app
### Data Analysis: Marissa (Meinizi) Zheng, Davit Sargsyan 
### Created: 07/22/2018 

---    

## Table of Contents
[File Legend](#leg)
[Daily Logs](#logs)  

## File Legend<a name="files"></a>
### Scripts
**app_v2.R**: current (08/01/2018) version of the app taht can run NextFlow script on Linux.    
**app_v1.R**: first version of the app. Can run batch files on Windows.    
**nextflow**: NextFlow program file.    
**pipeline.nf**: example NextFlow scriptwith passing arguments.    
**rna-seq-genecount.nf**: Bala's current (08/01/2018)  RNA pipeline script.    
**run.nextflow.R**: test running NextFlow script from R.     
**run_pgm.R**: test running NextFlow script from R and shinyFiles package.    
**test_run.bat**: test batch file to run from R.    
**testApp.R**: an example  of a Shiny Dashboard app.  

## Daily Logs<a name="logs"></a>
### 08/09/2018
* **app_v3** added DNA program. 

### 08/27/2018
* **app_v4** app_v4 add VennDiagram and heatmap to DEGseq analysis, add exploratory analysis with more general experiment design.

### 08/20/2018
* **app_v4** app_v4 add RNA DEGseq analysis.

### 08/13/2018
* **app_v3** app_v3 added DNA upstream pannel.

### 08/01/2018
* **app_v2** can now run NextFlow scripts (NextFlow must be installed on a Linux machine first). A dummy 'pipeline.nf' script is created and tested for passing parameter values.

### 07/27/2018
* Added folder path (using **shinyFiles** package), and ability to run batch file by clicking **Run Program** button.

### 07/26/2018
* modified RNA-seq pipeline uploaded        
* shiny app is modified with dashboard template     

### 07/21/2018
* Repository created       
* Template apps are uploaded    

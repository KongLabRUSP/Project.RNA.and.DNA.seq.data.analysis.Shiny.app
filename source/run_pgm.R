wd <- getwd()
pgm <- paste(wd,
             "/source/test_run.bat ABC 123 > ",
             wd,
             "/source/log.txt",
             sep = "")
pgm
system(pgm)

require(shinyFiles)
library(help = "shinyFiles")
?shinyFilesExample

shinyFilesExample()
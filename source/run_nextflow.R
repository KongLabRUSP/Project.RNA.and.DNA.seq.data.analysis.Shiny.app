# /home/administrator/nextflow run -process.echo true /home/administrator/Documents/Project.RNA.and.DNA.seq.data.analysis.Shiny.app/source/hello.nf
# /home/administrator/nextflow run -process.echo true hello.nf

system("/home/administrator/nextflow run -process.echo true /home/administrator/hello.nf")

pgm <- paste(getwd(),
             "/source/nextflow run -process.echo true ",
             getwd(),
             "/source/hello.nf",
             sep = "")
pgm
system(pgm)


# ./nextflow run pipeline.nf --db=/path/to/blast/db --query=/path/to/query.fasta

pgm <- paste(getwd(),
             "/source/nextflow run -process.echo true ",
             getwd(),
             "/source/pipeline.nf ",
             "--db=ABC --query=DEF",
             sep = "")
pgm
system(pgm)

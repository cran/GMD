
library(GMD)

## data
cat("\n\nRunning example in the data `cage': \n")
readline("Press `Enter' to continue..[ENTER]")
example(cage)

## gmd
cat("\n\nRunning example in the function `gmd': \n")
readline("Press `Enter' to continue..[ENTER]")
example(gmd)

## plot.gmd
cat("\n\nRunning example in the function `plot.gmd': \n")
readline("Press `Enter' to continue..[ENTER]")
example(plot.gmd)

## gmdm
cat("\n\nRunning example in the function `gmdm': \n")
readline("Press `Enter' to continue..[ENTER]")
example(gmdm)

## plot.gmdm
cat("\n\nRunning example in the function `plot.gmdm': \n")
readline("Press `Enter' to continue..[ENTER]")
example(plot.gmdm)

## citation
choice <- readline("\n\nPrinting the citation?[y/n]:")
if(choice!="n"){
  print(citation("GMD"))
}

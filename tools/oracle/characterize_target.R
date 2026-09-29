# Run the pinned, unchanged FilterTarget.R on each small geometry probe.
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
cases <- jsonlite::read_json(file.path(root,"cases.json"),simplifyVector=FALSE)
results <- list()
for (name in names(cases)) {
  folder <- file.path(root,"cases",name)
  command <- c("/work/source/excavator2/lib/R/FilterTarget.R",
               file.path(folder,"target.bed"),file.path(folder,"output"),"panel","synthetic","100",
               file.path(folder,"chromosomes.tsv"),file.path(folder,"gaps.tsv"))
  status <- system2("Rscript",command,stdout=file.path(folder,"stdout.log"),stderr=file.path(folder,"stderr.log"))
  results[[name]] <- list(status=as.integer(status))
}
jsonlite::write_json(results,file.path(root,"results.json"),auto_unbox=TRUE,pretty=TRUE)

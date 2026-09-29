# Run the complete unchanged inference script with centromere boundary variants.
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]; program <- args[2]; output <- args[3]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
cases <- list(both_arms=c(7000,11000), no_short=c(0,50),
  no_short_inside=c(50,500), no_long=c(19000,20000),
  single_short=c(150,175), single_long=c(17925,17950),
  endpoints=c(6000,12100), inside=c(5000,13000))
rows <- lapply(names(cases), function(name) {
  folder <- file.path(output,name)
  for (sample in c("Test1","Test2"))
    dir.create(file.path(folder,"Results",sample),recursive=TRUE,showWarnings=FALSE)
  bounds <- cases[[name]]
  centromeres <- file.path(folder,"centromeres.txt")
  writeLines(c("chrom\tchromStart\tchromEnd",
    paste("chr1",bounds[1],bounds[2],sep="\t"),"chr2\t7000\t11000"),centromeres)
  command <- c(file.path(program,"lib/R/EXCAVATORInferenceExome.R"),folder,
    file.path(root,"target"),file.path(root,"paired.yaml"),"paired",program,
    "synthetic",file.path(root,"prepared"),file.path(program,"parameters.yaml"),centromeres)
  log <- file.path(folder,"legacy.log")
  status <- system2("Rscript",command,stdout=log,stderr=log)
  data.frame(case=name,start=bounds[1],end=bounds[2],status=status)
})
write.table(do.call(rbind,rows),file.path(output,"cases.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)

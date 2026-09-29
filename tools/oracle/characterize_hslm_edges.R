# Pure R wrapper boundaries; do not pass undefined degenerate inputs to Fortran.
args <- commandArgs(trailingOnly=TRUE)
source(args[1])
dir.create(args[2], recursive=TRUE, showWarnings=FALSE)
encode <- function(x) paste(format(x, digits=17, trim=TRUE), collapse=",")
cases <- list(
  sole_below=c(0,3), sole_equal=c(0,4), sole_above=c(0,5),
  first_short=c(0,4,10), middle_short=c(0,6,10,16),
  last_short=c(0,6,10), adjacent_short=c(0,3,6,12),
  all_short=c(0,3,6,9), no_short=c(0,5,10)
)
rows <- lapply(names(cases), function(name) {
  breaks <- cases[[name]]
  values <- seq_len(tail(breaks,1))/8
  filtered <- FilterSeg(breaks,4)
  data.frame(case=name, breaks=encode(breaks), min_windows=4,
    values=encode(values), filtered=encode(filtered),
    reconstructed=encode(as.vector(SegResults(rbind(values),filtered))))
})
write.table(do.call(rbind,rows), file.path(args[2],"filter.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
profiles <- list(constant=rep(0,10), single=0.1, two=c(-1,1),
                 three=c(-1,0,1), varying=c(-1,-0.5,0,0.5,1))
rows <- lapply(names(profiles), function(name) {
  values <- profiles[[name]]
  p <- ParamEstSeq(rbind(values),0.1)
  data.frame(case=name, values=encode(values), mi=encode(p$mi),
             smu=encode(p$smu), sepsilon=encode(p$sepsilon))
})
write.table(do.call(rbind,rows), file.path(args[2],"parameters.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
message <- tryCatch({ SortState(1L); "unexpected success" },
                    error=function(e) conditionMessage(e))
writeLines(message,file.path(args[2],"single-window-error.txt"))

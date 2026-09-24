# Small algorithm fixture, not biological truth or a replacement acceptance sample.
args <- commandArgs(trailingOnly = TRUE)
source(args[1])
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
set.seed(20260924)
seed_before <- .Random.seed
mdata <- c(-3.1, -3, -2.9, -1.1, -1, -0.9, -0.1, 0, 0.1,
           0.5, 0.58, 0.7, 0.95, 1, 1.2)
thru <- 0.35
thrd <- 0.5
fit <- EMFastCall(mdata, thru, thrd)
muvec <- fit$muvec
sdvec <- fit$sdvec
prior <- fit$prior
bound <- fit$bound
iterations <- fit$iter
posterior <- PosteriorP(mdata, muvec, sdvec, prior)
calls <- LabelAss(posterior, mdata)
seed_after <- .Random.seed
stopifnot(identical(sort(unique(as.numeric(calls[,1]))), c(-2, -1, 0, 1, 2)))
save(mdata, thru, thrd, muvec, sdvec, prior, bound, iterations, posterior,
     calls, seed_before, seed_after, file = file.path(args[2], "fastcall.RData"))
write.table(cbind(mdata, calls), file.path(args[2], "calls.tsv"), sep = "\t",
            row.names = FALSE, col.names = c("log2_ratio", "call", "probability"), quote = FALSE)
writeLines(capture.output(sessionInfo(), RNGkind()), file.path(args[2], "runtime.txt"))

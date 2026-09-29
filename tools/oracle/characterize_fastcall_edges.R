# Extra differential cases. Legacy functions themselves are not edited.
args <- commandArgs(trailingOnly = TRUE)
source(args[1])
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
cases <- list(
  uneven = c(-3.1,-3,-2.9,-1.1,-1,-0.9,-0.1,0,0.1,0.5,0.58,0.7,0.95,1,1.2,0.03,0.07),
  boundaries = c(-50,-20,-1.5,-1.3,-0.5,0.35,0.9,20,50),
  all_normal = rep(0,7),
  single = 0,
  extreme = c(-100,-10,0,10,100),
  nondefault = c(-3.1,-3,-1.1,-1,-0.2,0,0.2,0.4,0.58,0.8,1,1.2)
)
if (length(args) >= 3) {
  hslm <- read.table(args[3], header=TRUE, sep="\t", colClasses="character")
  x <- as.numeric(hslm$SegMean[!is.na(hslm$Log2R)])
  cases$paired_baseline <- x[c(TRUE,diff(x)!=0)]
}
# Record MStep returns; tracing does not alter calculations or consume RNG.
.steps <- list()
trace("MStep", exit=quote({
  .steps[[length(.steps)+1]] <<- returnValue()
}), print=FALSE)
for (name in names(cases)) {
  .steps <- list()
  mdata <- cases[[name]]
  thru <- if (name == "nondefault") 0.2 else 0.35
  thrd <- if (name == "nondefault") 0.4 else 0.5
  start <- StartCond(mdata,thru,thrd)
  initial_deviations <- start$sdvec
  initial_priors <- c(0.05,0.1,0.7,0.1,0.05)
  initial_statistic <- sum(PosteriorP(mdata,start$muvec,start$sdvec,initial_priors)*initial_priors)
  set.seed(20260924)
  seed_before <- .Random.seed
  fit <- EMFastCall(mdata,thru,thrd)
  muvec <- fit$muvec; sdvec <- fit$sdvec; prior <- fit$prior
  bound <- fit$bound; iterations <- fit$iter
  posterior <- PosteriorP(mdata,muvec,sdvec,prior)
  calls <- LabelAss(posterior,mdata)
  seed_after <- .Random.seed
  iteration_trace <- do.call(rbind,lapply(.steps,function(step) {
    statistic <- sum(PosteriorP(mdata,step$mu,step$sdev,step$prior)*step$prior)
    c(step$sdev,step$prior,statistic)
  }))
  save(mdata,thru,thrd,muvec,sdvec,prior,bound,iterations,posterior,calls,
       seed_before,seed_after,iteration_trace,initial_deviations,initial_statistic,
       file=file.path(args[2],paste0(name,".RData")))
}
untrace("MStep")
# Prove trace instrumentation did not change the original fit.
for (name in names(cases)) {
  saved <- new.env()
  load(file.path(args[2],paste0(name,".RData")), envir=saved)
  plain <- EMFastCall(saved$mdata,saved$thru,saved$thrd)
  stopifnot(identical(plain$muvec,saved$muvec), identical(plain$sdvec,saved$sdvec),
            identical(plain$prior,saved$prior), identical(plain$iter,saved$iterations),
            identical(plain$bound,saved$bound))
}
# Direct assignment probes include near ties and losing intermediate ties.
posterior <- rbind(c(0.5,0.5,0,0,0), c(0,0,1,0,0), rep(0.2,5),
                   c(0.499998,0.500002,0,0,0), c(0.49999,0.50001,0,0,0),
                   c(0,0,0,0.5,0.5))
set.seed(20260924); seed_before <- .Random.seed
calls <- LabelAss(posterior,rep(0,nrow(posterior))); seed_after <- .Random.seed
save(posterior,calls,seed_before,seed_after,file=file.path(args[2],"ties.RData"))
# Cellularity correction is applied in the inference script before fitting.
original <- c(-10,-3,-1,0,0.58,1,2)
cellularity <- 0.6
corrected <- 2^original/cellularity - (1-cellularity)/cellularity
corrected[corrected < 2^-5] <- 2^-5
corrected <- log2(corrected)
save(original,cellularity,corrected,file=file.path(args[2],"cellularity.RData"))

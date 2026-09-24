# Capture the unchanged legacy numerical routines and Python-owned wrapper policy.
args <- commandArgs(trailingOnly=TRUE)
source(file.path(args[1], "lib/R/LibraryJSLMIn.R"))
dyn.load(file.path(args[1], "lib/F77/FastJointSLMLibraryI.so"))
dir.create(args[2], recursive=TRUE, showWarnings=FALSE)
cases <- list(
  shifts=c(rep(-0.8,8),rep(0,9),rep(0.7,8))+rep(c(-0.02,0,0.02,0.01,-0.01),5),
  short=c(-0.8,-0.7,0,0.1,0.7,0.8),
  uneven=c(-0.9,-0.8,-0.7,-0.1,0.02,0.1,0.7,0.9,0.8,0.1,-0.1),
  single=0.1,
  ties=rep(0,4)
)
for (name in names(cases)) {
  values <- cases[[name]]; n <- length(values)
  positions <- as.integer(cumsum(rep(c(100,100,30000,700000,1),length.out=n)))
  if (name == "uneven") positions[3] <- positions[2]
  omega <- 0.1; theta <- 0.0001; distance <- 100000
  # Single/constant profiles exercise kernels with externally estimated parameters.
  estimation_values <- if (name %in% c("single","ties")) cases$shifts else values
  parameters <- ParamEstSeq(rbind(estimation_values),omega)
  means <- as.vector(MukEst(rbind(values),1)); mi <- parameters$mi
  smu <- parameters$smu; sepsilon <- parameters$sepsilon
  if (name == "ties") { means <- c(-0.1,0.1); mi <- 0 }
  k <- length(means); eta <- theta+(1-theta)*exp(log(theta)/(diff(positions)/distance))
  z <- .Fortran("transemisi",as.double(means),as.double(mi),as.double(eta),as.integer(n-1),
       as.double(values),as.integer(k),as.integer(1),as.double(smu),as.double(sepsilon),
       as.integer(n),G=matrix(0,k,k),P=matrix(0,k,k*(n-1)),emission=matrix(0,k,n))
  G <- z$G; transitions <- z$P; emissions <- z$emission
  initial <- log(rep(1/k,k))
  v <- .Fortran("bioviterbii",initial,transitions,emissions,as.integer(n),as.integer(k),
                path=integer(n),predecessors=matrix(0L,k,n))
  path <- v$path; predecessors <- v$predecessors
  # Independently replay scores from original matrices; verify traceback decisions.
  scores <- matrix(0,k,n); scores[,1] <- initial+emissions[,1]
  if (n>1) for (t in 2:n) for (j in 1:k) {
    candidates <- scores[,t-1]+transitions[,(t-2)*k+j]
    stopifnot(which.max(candidates)==predecessors[j,t])
    scores[j,t] <- max(candidates)+emissions[j,t]
  }
  stopifnot(which.max(scores[,n])==path[n])
  breaks <- c(0,which(diff(path)!=0),n)
  fw <- 2
  filtered <- FilterSeg(breaks,fw)
  segmented <- as.vector(SegResults(rbind(values),filtered))
  if (n>1 && name!="ties") stopifnot(identical(as.numeric(breaks),
      as.numeric(JointSegIn(rbind(values),rbind(means),mi,smu,sepsilon,positions,omega,theta,distance))))
  save(values,positions,estimation_values,omega,theta,distance,means,mi,smu,sepsilon,eta,
       G,transitions,emissions,initial,path,predecessors,scores,breaks,fw,filtered,segmented,
       file=file.path(args[2],paste0(name,".RData")))
}

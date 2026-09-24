# Small non-normal, two-chromosome prepared-data fixture for analysis and writers.
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
dir.create(root,recursive=TRUE,showWarnings=FALSE)
chrom <- rep(c("chr1","chr2"),each=120)
pos <- rep(c(seq(100,6000,by=100),seq(12100,18000,by=100)),2)
start <- pos-20; end <- pos+20
class <- rep(c("IN","OUT"),120)
profile <- rep(rep(c(-1,-0.8,0,0.58,0.9,0),each=20),2)+rep(c(-0.02,0,0.02,0.01,-0.01),48)
for (sample in c("Control1","Control2","Test1","Test2")) {
  counts <- if (sample=="Control1") rep(100,240) else if (sample=="Control2") rep(110,240) else
            (if (sample=="Test1") 100 else 110)*2^profile
  MatrixNorm <- cbind(chrom,pos,start,end,paste0("gene",seq_along(pos)),counts,class)
  path <- file.path(root,"prepared",sample,"RCNorm")
  dir.create(path,recursive=TRUE,showWarnings=FALSE)
  save(MatrixNorm,file=file.path(path,paste0(sample,".NRC.RData")))
}
target <- file.path(root,"target");dir.create(file.path(target,"FRB"),recursive=TRUE,showWarnings=FALSE)
writeLines("chr1 chr2",file.path(target,"synthetic_chromosome.txt"))
writeLines(c("chrom\tchromStart\tchromEnd","chr1\t7000\t11000","chr2\t7000\t11000"),file.path(root,"centromeres.txt"))
for (chr in c("chr1","chr2")) {
  FRBData <- cbind(start[chrom==chr],rep(c("A","C","G","T"),30))
  save(FRBData,file=file.path(target,"FRB",paste0("FRB.",chr,".RData")))
}
writeLines(c("T2: Test2","C2: Control2","C1: Control1","T1: Test1"),file.path(root,"paired.yaml"))
# Inference pooling receives only test labels after DataAnalysisParallel splitting.
writeLines(c("T1: Test1","T2: Test2"),file.path(root,"tests.yaml"))
for (mode in c("paired","pooling")) {
  output <- file.path(root,mode)
  for (sample in c("Test1","Test2")) dir.create(file.path(output,"Results",sample),recursive=TRUE,showWarnings=FALSE)
  if (mode=="pooling") {
    load(file.path(root,"prepared/Control1/RCNorm/Control1.NRC.RData"))
    MatrixNorm[,6] <- as.character((as.numeric(MatrixNorm[,6])+110)/2)
    dir.create(file.path(output,"Control/RCNorm"),recursive=TRUE,showWarnings=FALSE)
    save(MatrixNorm,file=file.path(output,"Control/RCNorm/Control.NRC.RData"))
  }
  command <- c(file.path(args[2],"lib/R/EXCAVATORInferenceExome.R"),output,target,
               file.path(root,if (mode=="paired") "paired.yaml" else "tests.yaml"),mode,args[2],
               "synthetic",file.path(root,"prepared"),file.path(args[2],"parameters.yaml"),
               file.path(root,"centromeres.txt"))
  stopifnot(system2("Rscript",command)==0)
}

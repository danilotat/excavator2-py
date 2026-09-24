# Run the unchanged preparation and analysis scripts against tiny BAM fixtures.
args<-commandArgs(trailingOnly=TRUE)
root<-args[1]; program<-args[2]
legacy<-file.path(root,"legacy");target<-file.path(legacy,"target")
load(file.path(legacy,"prepared/Test1/RCNorm/Test1.NRC.RData"))
MyTarget<-MatrixNorm[,c(1,3,4,5,7)]
save(MyTarget,file=file.path(target,"synthetic.RData"))
for (folder in c("GCC","MAP")) dir.create(file.path(target,folder),showWarnings=FALSE)
for (chr in c("chr1","chr2")) {
  n<-sum(MyTarget[,1]==chr);GCContent<-rep(0.5,n);MapMed<-rep(1,n)
  save(GCContent,file=file.path(target,"GCC",paste0("GCC.",chr,".RData")))
  save(MapMed,file=file.path(target,"MAP",paste0("Map.",chr,".RData")))
}
prepared<-file.path(root,"legacy-prepared")
for (sample in c("Control1","Control2","Test1","Test2")) {
  folder<-file.path(prepared,sample)
  for (sub in c("RC","RCNorm","Images")) dir.create(file.path(folder,sub),recursive=TRUE,showWarnings=FALSE)
  stopifnot(system2("bash",c(file.path(program,"lib/bash/FiltBam.sh"),
     file.path(root,"bams",paste0(sample,".bam")),20,folder,program,sample,target))==0)
  stopifnot(system2("Rscript",c(file.path(program,"lib/R/EXCAVATORNormalizationExome.R"),
     program,folder,sample,target))==0)
  # Independent replay of the correction checkpoints, followed by a final-output check.
  source(file.path(program,"lib/R/LibraryExomeRC.R"))
  counts<-loadRC(file.path(folder,"RC"),list.files(file.path(folder,"RC")),c("chr1","chr2"))
  features<-loadTarget(target,c("chr1","chr2"))
  lengths<-as.integer(MyTarget[,3])-as.integer(MyTarget[,2])
  weighted<-counts/lengths;size<-mapped<-corrected<-normalized<-weighted
  for (kind in c("IN","OUT")) {
    i<-which(MyTarget[,5]==kind)
    size[i]<-if(kind=="IN") CorrectSize(weighted[i],lengths[i],5)$RCNorm else weighted[i]
    mapped[i]<-CorrectMAP(size[i],features[[3]][i]*100,5)$RCNorm
    corrected[i]<-CorrectGC(mapped[i],features[[2]][i]*100,5)$RCNorm
    normalized[i]<-corrected[i]
    normalized[i[which(normalized[i]==0)]]<-min(normalized[i][normalized[i]!=0])
  }
  load(file.path(folder,"RCNorm",paste0(sample,".NRC.RData")))
  stopifnot(identical(as.numeric(MatrixNorm[,6]),as.numeric(as.character(normalized))))
  save(counts,weighted,size,mapped,corrected,normalized,file=file.path(folder,"trace.RData"))
}
for (mode in c("paired","pooling")) {
  result<-file.path(root,paste0("legacy-prepared-",mode))
  for (sample in c("Test1","Test2")) dir.create(file.path(result,"Results",sample),recursive=TRUE,showWarnings=FALSE)
  if (mode=="pooling") stopifnot(system2("Rscript",c(file.path(program,"lib/R/PoolingCreateControl.R"),
       program,result,target,file.path(legacy,"paired.yaml"),mode,prepared))==0)
  stopifnot(system2("Rscript",c(file.path(program,"lib/R/EXCAVATORInferenceExome.R"),result,target,
      file.path(legacy,if(mode=="paired") "paired.yaml" else "tests.yaml"),mode,program,"synthetic",prepared,
      file.path(program,"parameters.yaml"),file.path(legacy,"centromeres.txt")))==0)
}

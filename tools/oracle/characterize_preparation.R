# Finite normalization fixtures and literal MakeReadCount wrapper probes.
args <- commandArgs(trailingOnly=TRUE)
source(file.path(args[1],"lib/R/LibraryExomeRC.R"))
dyn.load(file.path(args[1],"lib/F77/F4R.so"))
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
count_case <- function(reads, starts, ends, name) {
  nwin<-length(starts); residual<-0; jex<-1; counts<-rep(0,nwin)
  failed<-FALSE
  for (offset in seq(1,length(reads)+1,by=500000)) {
    rv<-reads[seq.int(offset,min(offset+499999,length(reads)))]
    if (offset>length(reads)) rv<-integer(0)
    n<-length(rv)
    if (!n) { failed<-TRUE;break }
    if (rv[n]<=ends[jex]) residual<-n
    else {
      out<-.Fortran("EXOMECOUNT",as.integer(rv),as.integer(nwin),as.integer(starts),
                    as.integer(ends),as.integer(n),as.integer(residual),as.integer(counts),as.integer(jex))
      residual<-out[[6]];counts<-out[[7]];jex<-out[[8]]
    }
    if (n<500000 || jex>nwin) break
  }
  save(reads,starts,ends,counts,residual,jex,failed,file=file.path(args[2],paste0(name,".RData")))
}
# Keep a final sentinel beyond all reads to avoid dereferencing out-of-bounds Fortran storage.
count_case(c(1L,10L,20L,21L,30L,40L,50L,80L),c(10L,20L,35L,100L),c(20L,25L,50L,200L),"endpoints")
count_case(c(10L,20L,21L,22L,50L,60L),c(10L,15L,30L,100L),c(25L,20L,50L,200L),"overlap")
count_case(c(110L,120L),c(100L),c(200L),"final_only")
for (n in c(499999L,500000L,500001L,1000001L)) {
  reads<-rep(15L,n)
  if (n>500000) reads[n]<-150L
  count_case(reads,c(10L,100L,300L),c(20L,200L,400L),paste0("chunk_",n))
}
counts<-c(0,1,2,10,20,40,3,7,9,0,8,2,9,12,13,14)
lengths<-c(5,10,12,5,10,12,20,20,30,30,40,40,50,50,60,60)
classes<-rep(c("IN","OUT"),8)
gc<-c(0,5,5.001,10,95,100,101,-1,20,25,30,40,50,60,70,80)/100
mappability<-rev(gc)
weighted<-counts/lengths;size<-mapped<-corrected<-normalized<-weighted
for (kind in c("IN","OUT")) {
  i<-which(classes==kind)
  size[i]<-if (kind=="IN") CorrectSize(weighted[i],lengths[i],5)$RCNorm else weighted[i]
  mapped[i]<-CorrectMAP(size[i],mappability[i]*100,5)$RCNorm
  corrected[i]<-CorrectGC(mapped[i],gc[i]*100,5)$RCNorm
  normalized[i]<-corrected[i]
  normalized[i[which(normalized[i]==0)]]<-min(normalized[i][normalized[i]!=0])
}
save(counts,lengths,classes,gc,mappability,weighted,size,mapped,corrected,normalized,
     file=file.path(args[2],"normalization.RData"))

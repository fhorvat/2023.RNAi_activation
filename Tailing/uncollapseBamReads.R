library(GenomicAlignments)

cat("\n @ pepap-functions loaded : \"uncollapseBamReads\"\n\n",sep="")

uncollapseBamReads <-
 function(
  BAMfile,xPAR=ScanBamParam(what=c("mapq","mrnm","mpos","isize","seq","qual"),tag=c("NH","HI","nM")),addUniqTag=F,UniqTagName="uniname",UniqTagSep="."
 ) {

 cat("\n ++ Loading BAMfile ++\n",sep="")
 cat("     ",BAMfile,"\n",sep="")
 xBAM <- readGAlignments( file=BAMfile,param=xPAR,use.names=T )

 if ( addUniqTag ) {
#  uni.I <- as.character(rep.int( x=seq_along(xBAM),times=as.integer(gsub( pattern="^[0-9]*[-]",replacement="",x=names(xBAM) )) ))
  uni.I <- as.integer(gsub( pattern="^[0-9]*[-]",replacement="",x=names(xBAM) ))
  uni.I <- unlist(sapply( X=uni.I,FUN=seq ))
 }

 cat(" ++ Uncollapsing ++\n",sep="")
 xBAM <-
  xBAM[
   rep.int(
    x     = seq_along(xBAM),
    times = as.integer( gsub( pattern="^[0-9]*[-]",replacement="",x=names(xBAM) ) )
   ) ]

 if ( addUniqTag ) {
  cat(" ++ Adding UniqTag ++\n",sep="")
  mcols(xBAM)[[UniqTagName]] <- paste( mcols(xBAM)[["qname"]],as.character(mcols(xBAM)[["HI"]]),uni.I,sep=UniqTagSep )
  names(xBAM)                <- mcols(xBAM)[[UniqTagName]]
 }

 cat(" ++ Finished ++\n\n",sep="")
 return(xBAM)

}


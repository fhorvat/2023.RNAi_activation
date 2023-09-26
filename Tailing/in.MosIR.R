library(GenomicRanges)
library(GenomicAlignments)
source("uncollapseBamReads.R")
source("clipped_seq.R")

#!#load( "/storage/brno1-cerit/home/pepap/DicerX/00.stats/02.DE/mirAnnot.dt.rda",verbose=T )
load( "mirAnnot.pepType.dt.rda",verbose=T )
mirAnnot.dt <- mirAnnot.pepType.dt

SUFF  <- ".se.Aligned.sortedByCoord.out.bam"
BSTOR <- "/path/to/brain/bams/"
BBAMS <-       system( command=paste0("/bin/ls ",BSTOR,"*/*",SUFF),intern=T )
BCONS <- gsub(  SUFF,"",gsub( "^.*[/]","",BBAMS ) )
HSTOR <- "/path/to/heart/bams/"
HBAMS <-       system( command=paste0("/bin/ls ",HSTOR,"*/*",SUFF),intern=T )
HCONS <- gsub(  SUFF,"",gsub( "^.*[/]","",HBAMS ) )
HCONS <- paste0( "Heart_",HCONS )
SSTOR <- "/path/to/spleen/bams/"
SBAMS <-       system( command=paste0("/bin/ls ",SSTOR,"*/*",SUFF),intern=T )
SCONS <- gsub(  SUFF,"",gsub( "^.*[/]","",SBAMS ) )
SCONS <- paste0( "Spleen_",SCONS )
LSTOR <- "/path/to/liver/bams/"
LBAMS <- head( system( command=paste0("/bin/ls ",LSTOR,"*/*",SUFF),intern=T ),15 )
LCONS <- gsub(  SUFF,"",gsub( "^.*[/]","",LBAMS ) )

aBAMS <- c(HBAMS,SBAMS,LBAMS,BBAMS)
aCONS <- c(HCONS,SCONS,LCONS,BCONS)

MINOVRL=15
#> remove spliced reads
RMV.SPL=as.logical("TRUE")
#> remove soft-clipped reads
RMV.SCL=as.logical("FALSE")
#> read-length range
RLENRAN=c(21,23)
#> alllowed number of edit distances
NMRANGE=c(0,100)
#> HI max allowed number
HI.LIM=NULL
#>
OUTNAME=paste0("sofClippedSeqs-MosIR-",RLENRAN[1],"-",RLENRAN[2],"nt")

#= FUN
dt2gr <- function(inDT,dt.seqnames="chr",dt.start="start",dt.end="end",dt.strand="strand",keep.metaData=T) {

 out.gr <-
  GRanges(
   seqnames = as.character( inDT[[dt.seqnames]] ),
   ranges   = IRanges(
    start   = as.integer( inDT[[dt.start]] ),
    end     = as.integer( inDT[[dt.end]] )
   ),
   strand   = as.character( inDT[[dt.strand]] )
  )

 if ( keep.metaData ) {
  mcols( out.gr ) <- inDT[ , colnames(inDT)[ !( colnames(inDT) %in% c(dt.seqnames,dt.start,dt.end,dt.strand) ) ] , with = F ]
 }

 return( out.gr )

}

#=

mirAnnot.gr <- dt2gr( inDT=mirAnnot.dt )

cat( "\n",sep="" )
if (F) {
objs.ls <- c()
for ( i in seq_along(aBAMS) ) {

 cat( " -> Loading sample : ",aCONS[i],"\n",sep="" )

 tmp.tga <- readGAlignments( file=aBAMS[i],param=ScanBamParam( what=c("qname","seq","qual"),which=GRanges("pCag-EGFP_mosIR:1-6660"),tag=c("NH","HI","nM","AS") ),use.names=T ) 

 if ( RMV.SPL ) { tmp.tga <- tmp.tga[ njunc(tmp.tga)  == 0              ] }
 if ( RMV.SCL ) { tmp.tga <- tmp.tga[ qwidth(tmp.tga) == width(tmp.tga) ] }

 tmp.tga     <- tmp.tga[ ( width(tmp.tga)         >= RLENRAN[1] ) & ( width(tmp.tga)         <= RLENRAN[2] ) ]
 tmp.tga     <- tmp.tga[ ( mcols(tmp.tga)[["nM"]] >= NMRANGE[1] ) & ( mcols(tmp.tga)[["nM"]] <= NMRANGE[2] ) ]

 tmp.cli <-
  clipped_seq(
   xBAM      = tmp.tga,
   MAX.HI    = HI.LIM,
   collapsed = T
  )
 tmp.cli[["INDEX"]] <- tmp.cli[,paste0(chromosome,":",coord3p,":",strand)]

 assign( x=paste0(aCONS[i],".cli.dt"),value=tmp.cli )
 objs.ls                            <- c( objs.ls,paste0(aCONS[i],".cli.dt") )

 rm( list=c("tmp.tga","tmp.cli") )
 gc( verbose=T )

 cat( "\n",sep="" )

}

 save( list=objs.ls,file=paste0(OUTNAME,"-",pepDate(),".rda") )

} else {
 objs.ls <- load( "sofClippedSeqs-MosIR-21-23nt-20230831.rda",verbose=T )
}

if (T) {

cat( "\n",sep="" )
scSeq.ls <- c()
for ( iob in objs.ls ) {

 cat( " => ",iob,"\n",sep="" )

 tmp.dt   <- get(iob)
 tmp.dt   <- tmp.dt[,{list( c=sum(rCount) )},by=c("aln.width","clipped_3p")][ order(c,decreasing=T) ]
 scSeq.ls <- unique( c(scSeq.ls,tmp.dt[ tmp.dt[ , (cumsum(c)/sum(c))<=0.99 ] , clipped_3p ]) )
 
 rm( list=c("tmp.dt") )
 gc( verbose=T )

}

cat( "\n",sep="" )
tot.dt <- data.table()
for ( iob in objs.ls ) {

 cat( " => ",iob,"\n",sep="" )

 tmp.dt   <- get(iob)
 tmp.dt[["SAMPLE"]] <- sub( "[.]cli[.]dt$","",iob )
 tot.dt <- rbind( tot.dt,tmp.dt[ clipped_3p %in% scSeq.ls ] )
 
 rm( list=c("tmp.dt") )
 gc( verbose=T )

}

tot.dt[ clipped_3p=="" ][["clipped_3p"]] <- "NONE"
tot.dcast.dt                             <- dcast( tot.dt,SAMPLE+aln.width~clipped_3p,value.var="rCount",fun=sum,drop=F,fill=0 )
tot.dcast.dt                             <- tot.dcast.dt[,c("SAMPLE","aln.width","NONE",colnames(tot.dcast.dt)[!(colnames(tot.dcast.dt)%in%c("SAMPLE","aln.width","NONE"))] ),with=F]

save( list=c("scSeq.ls","tot.dt","tot.dcast.dt"),file="combRobjs.MosIR.rda" )

library(pheatmap)
dir.create( path=paste0("heatmap-overview-",pepDate()),showWarnings=F )
for ( iS in aCONS ) {

 cat( " => ",iS,"\n",sep="" )

 pheatmap(
  mat=log10(as.matrix( x=tot.dcast.dt[ order(aln.width,decreasing=F) ][ SAMPLE==iS,colnames(tot.dcast.dt)!="SAMPLE",with=F ],rownames="aln.width" )+1),
  cluster_rows=F,cluster_cols=F,main=paste0(iS," : log10+1"),
  filename=paste0("heatmap-overview-",pepDate(),"/",iS,".log10-heatmap.pdf")
 )

}

}



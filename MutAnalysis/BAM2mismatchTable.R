library(data.table)
library(GenomicAlignments)
library(Biostrings)

cat("\n @ pepap-functions loaded : \"getMismatchTable\",\"MismatchTable2HeatMap\",\"ExtractMMregs\",\"getMUTMATRIX\",\"sumMatches\"\n\n",sep="")

#= (1)
getMismatchTable <-
 function(
  ITGA,FREF,
  COLLAPSED   = T,
  RETURN.TGA  = F,
  RETURN.PA   = F
 ) {

if ( !exists("REF.dss") ) {
 cat( "\n HINT : Create a DNAStringSet object \"REF.dss\" instead of loading reference sequence for each BAM-file\n\n",sep="" )
 cat( "\n ++ READING REFERENCE FASTA ++\n", sep="" )
 REF.dss        <- readDNAStringSet(FREF)
 names(REF.dss) <- sub( " .*$","",names(REF.dss) )
}

cat( " ++ READING BAM ALIGNMENT ++\n",      sep="" )
tmp.tga      <- ITGA
tmp.qseq.HI  <- mcols(tmp.tga)[["HI"]]
tmp.qseq.NH  <- mcols(tmp.tga)[["NH"]]
if ( RETURN.TGA ) { return(tmp.tga) }
tmp.rlen     <- width( tmp.tga )

cat( " ++ EXTRACT QUERY SEQUENCES ++\n",    sep="" )
tmp.qseq.str <- as.character(strand(tmp.tga))
if ( sum( !( unique(tmp.qseq.str) %in% c("-","+") ) )!=0 ) { print(table(tmp.qseq.str)); stop(" !!! Unknown strand found !!!\n") }
tmp.qseq.dss <- mcols(tmp.tga)[["seq"]]
if ( COLLAPSED ) {
tmp.qseq.cnt <- as.numeric( sub( "^[0-9]*[-]","",mcols(tmp.tga)[["qname"]] ) )
} else {
tmp.qseq.cnt <- as.numeric("1")
}
tmp.qseq.beg <- start(tmp.tga)
tmp.qseq.beg[ tmp.qseq.str=="-" ] <- end(tmp.tga)[ tmp.qseq.str=="-" ]
tmp.qseq.beg <- tmp.qseq.beg + c(-1,+1)[ as.integer(tmp.qseq.str=="-")+1 ]
cat( " ++ EXTRACT REFERENCE SEQUENCES ++\n",sep="" )
tmp.qseq.chr <- as.character(seqnames(tmp.tga))
tmp.rseq.ir  <- extractList( x=IRanges( start=start(tmp.tga),end=end(tmp.tga) ),i=IRanges( start=seq(length(tmp.tga)),width=1 ) )
tmp.rseq.dss <- unlist(extractAt( x=REF.dss[ tmp.qseq.chr ],at=tmp.rseq.ir ))

cat( " ++ SEQUENCE LAYER ++\n",             sep="" )
tmp.rseq.sL                      <- sequenceLayer( x=tmp.rseq.dss,cigar=cigar(tmp.tga),from="reference",to="pairwise-dense" )
tmp.rseq.sL[ tmp.qseq.str=="-" ] <- reverseComplement( x=tmp.rseq.sL[ tmp.qseq.str=="-" ] )
tmp.qseq.sL                      <- sequenceLayer( x=tmp.qseq.dss,cigar=cigar(tmp.tga),from="query",    to="pairwise-dense" )
tmp.qseq.sL[ tmp.qseq.str=="-" ] <- reverseComplement( x=tmp.qseq.sL[ tmp.qseq.str=="-" ] )
cat( " ++ PAIRWISE ALIGNMENT ++\n",         sep="" )
tmp.comp.pA <- pairwiseAlignment( pattern=tmp.qseq.sL,subject=tmp.rseq.sL,type="global",gapOpening=1000,gapExtension=1000 )

if ( RETURN.PA ) { return(tmp.comp.pA) }

cat( " ++ LIST MISMATCHES ++\n",            sep="" )
tmp.mism.dt               <- as.data.table( mismatchTable( x=tmp.comp.pA ) )
if ( nrow(tmp.mism.dt)!=0 ) {
 if ( nrow(tmp.mism.dt[ PatternStart  !=SubjectStart     ])!=0 ) { print(tmp.mism.dt[ PatternStart    !=SubjectStart     ]); stop(" !!! PA   start differs !!!\n") }
 if ( nrow(tmp.mism.dt[ PatternEnd    !=SubjectEnd       ])!=0 ) { print(tmp.mism.dt[ PatternEnd      !=SubjectEnd       ]); stop(" !!! PA     end differs !!!\n") }
 if ( nrow(tmp.mism.dt[ PatternQuality!=SubjectQuality   ])!=0 ) { print(tmp.mism.dt[ PatternQuality  !=SubjectQuality   ]); stop(" !!! PA quality differs !!!\n") }
 tmp.mism.dt               <- tmp.mism.dt[ , c("PatternId","SubjectSubstring","PatternSubstring","PatternStart") , with=F ]
 colnames(tmp.mism.dt)     <- c("ID","GENOMIC","SEQUENCED","GLOC")
 tmp.mism.dt[["COUNT"]]    <- tmp.qseq.cnt[ tmp.mism.dt[["ID"]] ]
 tmp.mism.dt[["HI"]]       <- tmp.qseq.HI[ tmp.mism.dt[["ID"]] ]
 tmp.mism.dt[["NH"]]       <- tmp.qseq.NH[ tmp.mism.dt[["ID"]] ]
 tmp.qseq.add              <- tmp.mism.dt[["GLOC"]]
 tmp.qseq.add[ tmp.qseq.str[ tmp.mism.dt[["ID"]] ]=="-" ] <- (-1)*tmp.qseq.add[ tmp.qseq.str[ tmp.mism.dt[["ID"]] ]=="-" ]
 tmp.mism.dt[["GLOC"]]     <- paste0( tmp.qseq.chr[ tmp.mism.dt[["ID"]] ],":",tmp.qseq.beg[ tmp.mism.dt[["ID"]] ]+tmp.qseq.add,":",tmp.qseq.str[ tmp.mism.dt[["ID"]] ] )
 tmp.mism.dt[["RLEN"]]     <- tmp.rlen[ tmp.mism.dt[["ID"]] ]
 tmp.mism.dt[["ID"]]       <- (mcols(tmp.tga)[["qname"]])[ tmp.mism.dt[["ID"]] ]
 tmp.mism.dt[["MISMATCH"]] <- tmp.mism.dt[ , paste0( GENOMIC,">",SEQUENCED ) ]
}

cat( " ++ DONE ++\n",                       sep="" )
return( tmp.mism.dt )

}

#= (2)
MismatchTable2HeatMap <- function( IMMTAB,STR.BIAS=F,NH.BIAS=F,CHR.BIAS=F,RETURN.CNTTAB=F ) {

 if ( (STR.BIAS+NH.BIAS+CHR.BIAS)>1 ) { stop( "\n !!! Multiple biases requested, choose only one !!!\n\n" ) }

 inp.dt  <-
  data.table(
   MISM = IMMTAB[ , rep.int( x =                MISMATCH,times=COUNT ) ],
   RLEN = IMMTAB[ , rep.int( x =                    RLEN,times=COUNT ) ],
   STR  = IMMTAB[ , rep.int( x = sub( "^.*[:]","",GLOC ),times=COUNT ) ],
   NH   = IMMTAB[ , rep.int( x =                      NH,times=COUNT ) ],
   CHR  = IMMTAB[ , rep.int( x = sub( "[:].*$","",GLOC ),times=COUNT ) ]
  )

 if ( RETURN.CNTTAB ) { return( inp.dt ) }

 out.dt <- table( inp.dt[          , c("MISM","RLEN") , with=F ] )

# if ( STR.BIAS ) {
#  fwd.dt <- table( inp.dt[ STR=="+" , c("MISM","RLEN") , with=F ] )
#  rev.dt <- table( inp.dt[ STR=="-" , c("MISM","RLEN") , with=F ] )
#  out.dt <- list( FWD=fwd.dt,REV=rev.dt )
# }
 if ( STR.BIAS ) {
  out.dt <- table( inp.dt[          , c("MISM", "STR") , with=F ] )
 }

 if (  NH.BIAS ) {
  out.dt <- table( inp.dt[          , c("MISM",  "NH") , with=F ] )
 }

 if ( CHR.BIAS ) {
  out.dt <- table( inp.dt[          , c("MISM", "CHR") , with=F ] )
 }

 return( out.dt )

}

#= (3)
ExtractMMregs <- function( IMMTAB,MIN.DIST=100,IGNORE.STRAND=F ) {

 igr <- GRanges( rep.int( x=IMMTAB[["GLOC"]],times=IMMTAB[["COUNT"]] ) )
 ogr <- reduce( x=igr,min.gapwidth=MIN.DIST,ignore.strand=IGNORE.STRAND )
 mcols(ogr)[["mmreg"]] <- sprintf( fmt=paste0("mmreg%0",nchar(length(ogr)),"d"),seq(length(ogr)) )

 ovl <- findOverlaps( query=igr,subject=ogr,ignore.strand=IGNORE.STRAND )
 odt <- merge( data.table( mmreg=mcols(ogr)[["mmreg"]] ),data.table( mmreg=(mcols(ogr)[["mmreg"]])[ subjectHits(ovl) ],count=as.integer("1") ),by="mmreg",sort=F,all.x=T )
 if ( nrow(odt[ is.na(count) ])!=0 ) {
  odt[ is.na(count) ][["count"]] <- as.integer("0")
 }
 odt <- odt[ , { list( count=sum(count) ) } , by="mmreg" ]
 odt <- merge( data.table( mmreg=mcols(ogr)[["mmreg"]],gloc=as.character(ogr),width=width(ogr) ),odt,by="mmreg",sort=F,all.x=T )
 odt <- odt[ order(count,width,decreasing=T) ]

 return( odt )

}

#= (4)
getMUTMATRIX <- function( IMMTAB,RLENRANGE=c(0,100),NUCS=c("A","C","G","T") ) {

 if ( nrow(IMMTAB)==0 ) {

  return( matrix( data=0,     nrow=length(NUCS),ncol=length(NUCS),dimnames=list( NUCS,NUCS ) ) )

 } else                 {

  MUTS <- c()
  for( i in NUCS ){ for( j in NUCS ){ MUTS <- c( MUTS,paste0( i,">",j ) ) } }

  IMMTAB <- IMMTAB[ ( RLEN >= RLENRANGE[1] ) & ( RLEN <= RLENRANGE[2] ) ]

  OUTMAT                  <- IMMTAB[,{ list(mm.sum=sum(COUNT)) },by="MISMATCH"]
  OUTMAT                  <- setNames( object=OUTMAT[["mm.sum"]],nm=OUTMAT[["MISMATCH"]] )
  OUTMAT                  <- OUTMAT[ MUTS ]
  OUTMAT[ is.na(OUTMAT) ] <- 0

  return( matrix( data=OUTMAT,nrow=length(NUCS),ncol=length(NUCS),dimnames=list( NUCS,NUCS ),byrow=T ) )

 }

}

#= (1)
sumMatches <- function( ITGA,FREF,COLLAPSED=T,RETURN.PA   = T ) {

if ( !exists("REF.dss") ) {
 cat( "\n HINT : Create a DNAStringSet object \"REF.dss\" instead of loading reference sequence for each BAM-file\n\n",sep="" )
 cat( "\n ++ READING REFERENCE FASTA ++\n", sep="" )
 REF.dss        <- readDNAStringSet(FREF)
 names(REF.dss) <- sub( " .*$","",names(REF.dss) )
}

cat( " ++ READING BAM ALIGNMENT ++\n",      sep="" )
tmp.tga      <- ITGA
tmp.qseq.HI  <- mcols(tmp.tga)[["HI"]]
tmp.qseq.NH  <- mcols(tmp.tga)[["NH"]]
tmp.rlen     <- width( tmp.tga )

cat( " ++ EXTRACT QUERY SEQUENCES ++\n",    sep="" )
tmp.qseq.str <- as.character(strand(tmp.tga))
if ( sum( !( unique(tmp.qseq.str) %in% c("-","+") ) )!=0 ) { print(table(tmp.qseq.str)); stop(" !!! Unknown strand found !!!\n") }
tmp.qseq.dss <- mcols(tmp.tga)[["seq"]]
if ( COLLAPSED ) {
tmp.qseq.cnt <- as.numeric( sub( "^[0-9]*[-]","",mcols(tmp.tga)[["qname"]] ) )
} else {
tmp.qseq.cnt <- as.numeric("1")
}
tmp.qseq.beg <- start(tmp.tga)
tmp.qseq.beg[ tmp.qseq.str=="-" ] <- end(tmp.tga)[ tmp.qseq.str=="-" ]
tmp.qseq.beg <- tmp.qseq.beg + c(-1,+1)[ as.integer(tmp.qseq.str=="-")+1 ]
cat( " ++ EXTRACT REFERENCE SEQUENCES ++\n",sep="" )
tmp.qseq.chr <- as.character(seqnames(tmp.tga))
tmp.rseq.ir  <- extractList( x=IRanges( start=start(tmp.tga),end=end(tmp.tga) ),i=IRanges( start=seq(length(tmp.tga)),width=1 ) )
tmp.rseq.dss <- unlist(extractAt( x=REF.dss[ tmp.qseq.chr ],at=tmp.rseq.ir ))

cat( " ++ SEQUENCE LAYER ++\n",             sep="" )
tmp.rseq.sL                      <- sequenceLayer( x=tmp.rseq.dss,cigar=cigar(tmp.tga),from="reference",to="pairwise-dense" )
tmp.rseq.sL[ tmp.qseq.str=="-" ] <- reverseComplement( x=tmp.rseq.sL[ tmp.qseq.str=="-" ] )
tmp.qseq.sL                      <- sequenceLayer( x=tmp.qseq.dss,cigar=cigar(tmp.tga),from="query",    to="pairwise-dense" )
tmp.qseq.sL[ tmp.qseq.str=="-" ] <- reverseComplement( x=tmp.qseq.sL[ tmp.qseq.str=="-" ] )
cat( " ++ PAIRWISE ALIGNMENT ++\n",         sep="" )
tmp.comp.pA <- pairwiseAlignment( pattern=tmp.qseq.sL,subject=tmp.rseq.sL,type="global",gapOpening=1000,gapExtension=1000 )

tmp.match.dt               <-
 data.table(
  xA = sum(
   ( unlist(strsplit( x=rep.int( x=as.character(tmp.comp.pA@pattern),times=tmp.qseq.cnt ),split="",fixed=T )) == "A" ) &
   ( unlist(strsplit( x=rep.int( x=as.character(tmp.comp.pA@pattern),times=tmp.qseq.cnt ),split="",fixed=T )) == "A" )
       ),
  xC = sum(
   ( unlist(strsplit( x=rep.int( x=as.character(tmp.comp.pA@pattern),times=tmp.qseq.cnt ),split="",fixed=T )) == "C" ) &
   ( unlist(strsplit( x=rep.int( x=as.character(tmp.comp.pA@pattern),times=tmp.qseq.cnt ),split="",fixed=T )) == "C" )
       ),
  xG = sum(
   ( unlist(strsplit( x=rep.int( x=as.character(tmp.comp.pA@pattern),times=tmp.qseq.cnt ),split="",fixed=T )) == "G" ) &
   ( unlist(strsplit( x=rep.int( x=as.character(tmp.comp.pA@pattern),times=tmp.qseq.cnt ),split="",fixed=T )) == "G" )
       ),
  xT = sum(
   ( unlist(strsplit( x=rep.int( x=as.character(tmp.comp.pA@pattern),times=tmp.qseq.cnt ),split="",fixed=T )) == "T" ) &
   ( unlist(strsplit( x=rep.int( x=as.character(tmp.comp.pA@pattern),times=tmp.qseq.cnt ),split="",fixed=T )) == "T" )
       )
 )

cat( " ++ DONE ++\n",                       sep="" )
return( tmp.match.dt )

}


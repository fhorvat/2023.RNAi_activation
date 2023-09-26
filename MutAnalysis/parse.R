library(data.table)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)
source("BAM2mismatchTable.R")

ANNOTF="gc_vM15_RM_all.chr_patch_hapl_scaff.rda"

INFA  <- "GRCm38.p6.genome.pCag-EGFP_MosIR.fa.gz"
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

#> MosIR crds ... pCag-EGFP_mosIR:2508-3027 / pCag-EGFP_mosIR:3126-3645
#!#MosIR.gr <- GRanges( c("pCag-EGFP_mosIR:2508-3027:*","pCag-EGFP_mosIR:3126-3645:*") )
MosIR.gr <- GRanges( c("pCag-EGFP_mosIR:2508-3645:*") )

ordN <-
 c(

 "Heart_H13_WT","Heart_H14_WT","Heart_H15_WT",                   
 "Heart_H7_DcrHASom_EGFPmosIR_delPKR","Heart_H8_DcrHASom_EGFPmosIR_delPKR","Heart_H9_DcrHASom_EGFPmosIR_delPKR",
 "Heart_H10_DcrHASom_EGFPmosIR_PKR","Heart_H11_DcrHASom_EGFPmosIR_PKR","Heart_H12_DcrHASom_EGFPmosIR_PKR",
 "Heart_H1_DcrHAoo_EGFPmosIR_delPKR","Heart_H2_DcrHAoo_EGFPmosIR_delPKR","Heart_H3_DcrHAoo_EGFPmosIR_delPKR",
 "Heart_H4_DcrHAoo_EGFPmosIR_PKR","Heart_H5_DcrHAoo_EGFPmosIR_PKR","Heart_H6_DcrHAoo_EGFPmosIR_PKR",

 "Spleen_wtDicer_mosIR_PKR_01","Spleen_wtDicer_mosIR_PKR_02","Spleen_WT",
 "Spleen_DicerSom_mosIR_delPKR_01","Spleen_DicerSom_mosIR_delPKR_02","Spleen_DicerSom_mosIR_delPKR_03",
 "Spleen_DicerSom_mosIR_PKR_01","Spleen_DicerSom_mosIR_PKR_02","Spleen_DicerSom_mosIR_PKR_03",
 "Spleen_DicerX_mosIR_delPKR_01","Spleen_DicerX_mosIR_delPKR_02","Spleen_DicerX_mosIR_delPKR_03",
 "Spleen7_WT","Spleen8_WT","Spleen9_WT",
 "Spleen1_DcrHAoo_EGFPmosIR_delPKR","Spleen2_DcrHAoo_EGFPmosIR_delPKR","Spleen3_DcrHAoo_EGFPmosIR_delPKR",
 "Spleen4_DcrHAoo_EGFPmosIR_PKR","Spleen5_DcrHAoo_EGFPmosIR_PKR","Spleen6_DcrHAoo_EGFPmosIR_PKR",

 "Liver_WT01","Liver_WT02","Liver_WT03",
 "Liver_DcrHASom_EGFPmosIR_delPKR01","Liver_DcrHASom_EGFPmosIR_delPKR02","Liver_DcrHASom_EGFPmosIR_delPKR03",
 "Liver_DcrHASom_EGFPmosIR_PKR01","Liver_DcrHASom_EGFPmosIR_PKR02","Liver_DcrHASom_EGFPmosIR_PKR03",
 "Liver_DcrHAoo_EGFPmosIR_delPKR01","Liver_DcrHAoo_EGFPmosIR_delPKR02","Liver_DcrHAoo_EGFPmosIR_delPKR03",
 "Liver_DcrHAoo_EGFPmosIR_PKR01","Liver_DcrHAoo_EGFPmosIR_PKR02","Liver_DcrHAoo_EGFPmosIR_PKR03",

 "Brain07_WT","Brain08_WT","Brain09_WT",
 "Brain01_DcrHAoo_EGFPmosIR_delPKR","Brain02_DcrHAoo_EGFPmosIR_delPKR","Brain03_DcrHAoo_EGFPmosIR_delPKR",
 "Brain04_DcrHAoo_EGFPmosIR_PKR","Brain05_DcrHAoo_EGFPmosIR_PKR","Brain06_DcrHAoo_EGFPmosIR_PKR",

 "Thymus7_WT","Thymus8_WT","Thymus9_WT",
 "Thymus1_DcrHAoo_EGFPmosIR_delPKR","Thymus2_DcrHAoo_EGFPmosIR_delPKR","Thymus3_DcrHAoo_EGFPmosIR_delPKR",
 "Thymus4_DcrHAoo_EGFPmosIR_PKR","Thymus5_DcrHAoo_EGFPmosIR_PKR","Thymus6_DcrHAoo_EGFPmosIR_PKR",

 "B1_DcrSom_EGFPmosIR_delPKR","B2_DcrSom_EGFPmosIR_delPKR","B3_DcrSom_EGFPmosIR_delPKR",
 "Sp1_DcrSom_EGFPmosIR_delPKR","Sp2_DcrSom_EGFPmosIR_delPKR","Sp3_DcrSom_EGFPmosIR_delPKR",
 "Thy1_DcrSom_EGFPmosIR_delPKR","Thy2_DcrSom_EGFPmosIR_delPKR","Thy3_DcrSom_EGFPmosIR_delPKR"

 )

#> remove spliced reads
RMV.SPL=as.logical("TRUE")
#> remove soft-clipped reads
RMV.SCL=as.logical("FALSE")
#> read-length range
RLENRAN=c(21,23)
#> alllowed number of edit distances
NMRANGE=c(0,100)
#> alllowed number of multimapper hits
HIRANGE=c(1,100)
#>
OUTNAME=paste0("mutAnal.",RLENRAN[1],"to",RLENRAN[2],"nt")

#>># LOADING GENOMIC DNA SEQUENCE
if ( !exists("REF.dss") ) {
 cat( "\n ++ READING REFERENCE FASTA ++\n\n", sep="" )
 REF.dss        <- readDNAStringSet(INFA)
 names(REF.dss) <- sub( " .*$","",names(REF.dss) )
#>># correct found mismatch in reference
 REF.dss[["pCag-EGFP_mosIR"]][3352] <- "C"
}

#>># LOADING ANNOTATIONS
if ( !exists("mirAnnot.pepType.dt") ) {
 load( "mirAnnot.pepType.dt.rda",verbose=T )
}
if ( !exists("gc_annot.dt") ) {

 gc_annot.gr <- load( file=ANNOTF,verbose=T)
 gc_annot.gr <- get(gc_annot.gr)
 gc_annot.dt <- as.data.table(gc_annot.gr)
 gid2name.dt <- unique( gc_annot.dt[,c("gene_id","gene_name"),with=F] )

}
cat("\n",sep="")

#>># EXTRACTING MISMATCHES
if (T) {

 out.Nreads.dt  <- data.table()
 out.mmTabs.obj <- list()
 out.xMatch.obj <- list()
 out.mmMATs.obj <- list()
 for ( ibam in seq_along(aBAMS) ) {

  cat( " => ",aCONS[ibam],"\n",sep="" )
#!#  tmp.tga <-
#!#   uncollapseBamReads(
#!#    BAMfile     = aBAMS[ibam],
#!#    xPAR        = ScanBamParam( what=c("qname","seq","qual"),which=GRanges("pCag-EGFP_mosIR:1-6660"),tag=c("NH","HI","nM","AS") ),
#!#    addUniqTag  = T,
#!#    UniqTagName = "uniname",
#!#    UniqTagSep  = "."
#!#   )
#!#  tmp.tga <- readGAlignments( file=aBAMS[ibam],param=ScanBamParam( what=c("qname","seq","qual"),tag=c("NH","HI","nM","AS"),which=GRanges("pCag-EGFP_mosIR:1-6660") ),use.names=T ) 
  tmp.tga <- readGAlignments( file=aBAMS[ibam],param=ScanBamParam( what=c("qname","seq","qual"),tag=c("NH","HI","nM","AS"),which=MosIR.gr ),use.names=T ) 

  if ( RMV.SPL ) { tmp.tga <- tmp.tga[ njunc(tmp.tga)  == 0              ] }
  if ( RMV.SCL ) { tmp.tga <- tmp.tga[ qwidth(tmp.tga) == width(tmp.tga) ] }

  tmp.tga     <- tmp.tga[ ( width(tmp.tga)         >= RLENRAN[1] ) & ( width(tmp.tga)         <= RLENRAN[2] ) ]
  tmp.tga     <- tmp.tga[ as.character(strand(tmp.tga))=="+" ]
  print( length(tmp.tga) )
  tmp.tga     <- tmp.tga[ !duplicated(names(  tmp.tga))      ]
  print( length(tmp.tga) )

  if ( length(tmp.tga)!=0 ) {
   imatch.dt   <- sumMatches( ITGA=tmp.tga,FREF=INFA,COLLAPSED=T )
   NREADS      <- unique(mcols(tmp.tga)[["qname"]])
   NREADS      <- sum(as.numeric( sub( "^[0-9]*[-]","",NREADS ) ))
  } else                    {
   imatch.dt   <- data.table( xA=0,xC=0,xG=0,xT=0 )
   NREADS      <- 0
  }

  if ( length(tmp.tga[ mcols(tmp.tga)[["nM"]]!=0 ])!=0 ) {

   MMREADS     <- unique(mcols(tmp.tga[ mcols(tmp.tga)[["nM"]]!=0 ])[["qname"]])
   MMREADS     <- sum(as.numeric( sub( "^[0-9]*[-]","",MMREADS ) ))
   tmp.tga     <- tmp.tga[ ( mcols(tmp.tga)[["nM"]] >= NMRANGE[1] ) & ( mcols(tmp.tga)[["nM"]] <= NMRANGE[2] ) ]
#!#   tmp.tga     <- tmp.tga[ ( mcols(tmp.tga)[["HI"]] >= HIRANGE[1] ) & ( mcols(tmp.tga)[["HI"]] <= HIRANGE[2] ) ]

   imtb.dt                                               <-
    getMismatchTable(
     ITGA      = tmp.tga,
     FREF      = INFA,
     COLLAPSED = T
    )

  } else                    {

   MMREADS     <- 0
   imtb.dt     <- data.table()

  }

  out.mmTabs.obj                                        <- c( out.mmTabs.obj,list( NNN=imtb.dt ) )
  names(out.mmTabs.obj)[ names(out.mmTabs.obj)=="NNN" ] <- aCONS[ibam]

  out.Nreads.dt                                         <- rbind( out.Nreads.dt,data.table( LIB=aCONS[ibam],TOT.READS=NREADS,MISMATCHED.READS=MMREADS ) )

  if ( nrow(imtb.dt)!=0 ) {
  tmp.mat                                               <- getMUTMATRIX( IMMTAB=imtb.dt[ GRanges( GLOC ) %within% MosIR.gr ] )
  } else                  {
  tmp.mat                                               <- getMUTMATRIX( IMMTAB=imtb.dt                                      )
  }
  diag(tmp.mat)                                         <- unlist(imatch.dt)[c("xA","xC","xG","xT")]
  out.mmMATs.obj                                        <- c( out.mmMATs.obj,list( NNN=tmp.mat ) )
  names(out.mmMATs.obj)[ names(out.mmMATs.obj)=="NNN" ] <- aCONS[ibam]

  out.xMatch.obj                                        <- c( out.xMatch.obj,list( NNN=imatch.dt ) )
  names(out.xMatch.obj)[ names(out.xMatch.obj)=="NNN" ] <- aCONS[ibam]

  cat( "\n",sep="" )

 }

 assign(  x=paste0(OUTNAME,".mmTabs.obj"),value=out.mmTabs.obj )
 save( list=paste0(OUTNAME,".mmTabs.obj"),file=paste0(OUTNAME,"-",pepDate(),".mmTabs.obj.rda") )
 assign(  x=paste0(OUTNAME,".Nreads.dt"), value=out.Nreads.dt  )
 save( list=paste0(OUTNAME,".Nreads.dt"), file=paste0(OUTNAME,"-",pepDate(),".Nreads.dt.rda")  )
 assign(  x=paste0(OUTNAME,".mmMATS.obj"),value=out.mmMATs.obj )
 save( list=paste0(OUTNAME,".mmMATS.obj"),file=paste0(OUTNAME,"-",pepDate(),".mmMATS.obj.rda") )
 assign(  x=paste0(OUTNAME,".xMatch.obj"),value=out.xMatch.obj )
 save( list=paste0(OUTNAME,".xMatch.obj"),file=paste0(OUTNAME,"-",pepDate(),".xMatch.obj.rda") )

}

library(openxlsx)
TISSUE="none"
out.wb <- createWorkbook()
addWorksheet(   wb=out.wb,sheetName="TOTAL.READ.NUMBER" )
writeDataTable( wb=out.wb,sheet=    "TOTAL.READ.NUMBER", x=out.Nreads.dt,   colNames=T,rowNames=F,headerStyle=createStyle(textDecoration="bold") )
for ( NS in ordN ) {

 cat( " >> ",NS," <<\n",sep="" )

 out.dt <- out.mmMATs.obj[[NS]]
 cTIS   <- substr( x=NS,start=1,stop=5 )
 NSHT   <- cTIS
 if ( NSHT=="Splee" ) { NSHT <- "Spleen" }
 if ( NSHT=="Thymu" ) { NSHT <- "Thymus" }
 if ( cTIS=="Sm6_w" ) { REP  <- 1; LINE <- 15 }
 if ( grepl(pattern="^Sm",x=NSHT) | grepl(pattern="^B[1-3]",x=NSHT) | grepl(pattern="^Sp[1-3]",x=NSHT) | grepl(pattern="^Thy[1-3]",x=NSHT) ) { NSHT <- "LastKohort"; cTIS <- "LastKohort" }
 if ( TISSUE!=cTIS ) {
  addWorksheet(  wb=out.wb,sheetName=NSHT )
  TISSUE <- cTIS
  REP    <- 1
  LINE   <- 1
 }
 writeData(      wb=out.wb,sheet=    NSHT, x=NS,                          startCol=1+(7*(REP-1)),startRow=LINE                    )
 writeData(      wb=out.wb,sheet=    NSHT, x=out.dt,colNames=T,rowNames=T,startCol=1+(7*(REP-1)),startRow=LINE+1+(floor((REP-1)/3)*7) )
 REP <- REP+1
 if ( REP>3 ) {
  REP  <- 1
  LINE <- LINE+7
 }

}

saveWorkbook(    wb=out.wb,file=paste0( OUTNAME,"-",pepDate(),".MutMatrices.xlsx" ),overwrite=T )


library(data.table)

cat("\n @ pepap-functions loaded : \"clipped_seq\"\n\n",sep="")

clipped_seq <- function(xBAM,MAX.HI=NULL,collapsed=F) {

if ( sum(names(elementMetadata(xBAM)) %in% "seq")   == 0 ) { stop("\n\nNo \"seq\" slot found in input xBAM\n\n") }
if ( sum(names(elementMetadata(xBAM)) %in% "qname") == 0 ) { stop("\n\nNo \"qname\" slot found in input xBAM\n\n") }
if ( sum(names(elementMetadata(xBAM)) %in% "HI")    == 0 ) { stop("\n\nNo \"HI\" slot found in input xBAM\n\n") }
if ( sum(names(elementMetadata(xBAM)) %in% "NH")    == 0 ) { stop("\n\nNo \"NH\" slot found in input xBAM\n\n") }

cat( "\n ++ GenomicAlignment => data.table\n",sep="" )
if ( collapsed ) {

xDT <-
 data.table(qname      = as.character(xBAM@elementMetadata[,"qname"]),
            rCount     = as.numeric( gsub( pattern = "^[0-9]*[-]", replacement = "", x = mcols(xBAM)[ , "qname" ] ) ),
            cigar      = as.character(cigar(xBAM)),
            red_cigar  = as.character(gsub(pattern = "[0-9]", replacement = "", x = cigar(xBAM))),
            aln.width  = width(xBAM),
            chromosome = as.character(seqnames(xBAM)),
            strand     = as.character(  strand(xBAM)),
            rStart     = start(xBAM),
            rEnd       =   end(xBAM),
            clipped_5p = as.character(""),
            clipped_3p = as.character(""),
            seq        = as.character(xBAM@elementMetadata[,"seq"]),
            HI         = as.integer(  xBAM@elementMetadata[,"HI"]),
            NH         = as.integer(  xBAM@elementMetadata[,"NH"]))
setkey(x = xDT, "qname")

if ( is.null(MAX.HI) ) {
xDT <- xDT[ , { minHI = min(HI)
               list(rCount     = rCount,
                    cigar      = cigar,
                    red_cigar  = red_cigar,
                    aln.width  = aln.width,
                    chromosome = chromosome,
                    strand     = strand,
                    rStart     = rStart,
                    rEnd       = rEnd,
                    clipped_5p = clipped_5p,
                    clipped_3p = clipped_3p,
                    seq        = seq,
                    HI         = HI,
                    NH         = NH,
                    minHI      = minHI)
              } , by = qname ]
#! no HI filter !#xDT <- xDT[ HI == minHI, ]
} else                 {
xDT <- xDT[ HI <= MAX.HI , ]
}

} else           {

xDT <-
 data.table(qname      = as.character(xBAM@elementMetadata[,"qname"]),
            cigar      = as.character(cigar(xBAM)),
            red_cigar  = as.character(gsub(pattern = "[0-9]", replacement = "", x = cigar(xBAM))),
            aln.width  = width(xBAM),
            chromosome = as.character(seqnames(xBAM)),
            strand     = as.character(  strand(xBAM)),
            rStart     = start(xBAM),
            rEnd       =   end(xBAM),
            clipped_5p = as.character(""),
            clipped_3p = as.character(""),
            seq        = as.character(xBAM@elementMetadata[,"seq"]),
            HI         = as.integer(  xBAM@elementMetadata[,"HI"]),
            NH         = as.integer(  xBAM@elementMetadata[,"NH"]))
setkey(x = xDT, "qname")

if ( is.null(MAX.HI) ) {
xDT <- xDT[, { minHI = min(HI)
               list(cigar      = cigar,
                    red_cigar  = red_cigar,
                    aln.width  = aln.width,
                    chromosome = chromosome,
                    strand     = strand,
                    rStart     = rStart,
                    rEnd       = rEnd,
                    clipped_5p = clipped_5p,
                    clipped_3p = clipped_3p,
                    seq        = seq,
                    HI         = HI,
                    NH         = NH,
                    minHI      = minHI)}, by = qname]
#! no HI filter !#xDT <- xDT[ HI == minHI, ]
} else                 {
xDT <- xDT[ HI <= MAX.HI, ]
}

}

cat( " ++ Get 5'/3' coordinates\n",sep="" )
xDT[["coord5p"]] <- as.integer("")
xDT[["coord3p"]] <- as.integer("")
cat( "   > Forward strand\n",sep="" )
if ( nrow(xDT[ strand == "+" ]) != 0 ) {
cat( "   >> 5' coordinates\n",sep="" )
xDT[ strand == "+" ][["coord5p"]] <- xDT[ strand == "+" ][["rStart"]]
cat( "   >> 3' coordinates\n",sep="" )
xDT[ strand == "+" ][["coord3p"]] <- xDT[ strand == "+" ][["rEnd"]]
}
cat( "   > Reverse strand\n",sep="" )
if ( nrow(xDT[ strand == "-" ]) != 0 ) {
cat( "   >> 5' coordinates\n",sep="" )
xDT[ strand == "-" ][["coord5p"]] <- xDT[ strand == "-" ][["rEnd"]]
cat( "   >> 3' coordinates\n",sep="" )
xDT[ strand == "-" ][["coord3p"]] <- xDT[ strand == "-" ][["rStart"]]
}

cat( " ++ Get clipped sequences\n",sep="" )
cat( "   > Forward strand\n",sep="" )
#= forward Sx
cat( "   >> 5' clipped only\n",sep="" )
tmpDT<-xDT[  grepl(pattern="^S",x=red_cigar) & !grepl(pattern="S$",x=red_cigar) & strand=="+",]
if ( nrow(tmpDT) != 0 ) {
       xDT[  grepl(pattern="^S",x=red_cigar) & !grepl(pattern="S$",x=red_cigar) & strand=="+",][["clipped_5p"]] <-
 as.character(                   unlist(extractAt(x=DNAStringSet(tmpDT[,seq]),at=cigarRangesAlongQuerySpace(cigar=tmpDT[,cigar],ops="S"))))
}
#= forward  xS
cat( "   >> 3' clipped only\n",sep="" )
tmpDT<-xDT[ !grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="+",]
if ( nrow(tmpDT) != 0 ) {
       xDT[ !grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="+",][["clipped_3p"]] <-
 as.character(                   unlist(extractAt(x=DNAStringSet(tmpDT[,seq]),at=cigarRangesAlongQuerySpace(cigar=tmpDT[,cigar],ops="S"))))
}
#= forward SxS
cat( "   >> both ends clipped\n",sep="" )
tmpDT<-xDT[  grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="+",]
if ( nrow(tmpDT) != 0 ) {
       xDT[  grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="+",][["clipped_5p"]] <-
 as.character(                   unlist(extractAt(x=DNAStringSet(tmpDT[,seq]),at=cigarRangesAlongQuerySpace(cigar=tmpDT[,cigar],ops="S"))))[seq(from=1,to=(2*nrow(tmpDT)),by=2)]
}
tmpDT<-xDT[  grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="+",]
if ( nrow(tmpDT) != 0 ) {
       xDT[  grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="+",][["clipped_3p"]] <-
 as.character(                   unlist(extractAt(x=DNAStringSet(tmpDT[,seq]),at=cigarRangesAlongQuerySpace(cigar=tmpDT[,cigar],ops="S"))))[seq(from=2,to=(2*nrow(tmpDT)),by=2)]
}

cat( "   > Reverse strand\n",sep="" )
#= reverse Sx
cat( "   >> 3' clipped only\n",sep="" )
tmpDT<-xDT[  grepl(pattern="^S",x=red_cigar) & !grepl(pattern="S$",x=red_cigar) & strand=="-",]
if ( nrow(tmpDT) != 0 ) {
       xDT[  grepl(pattern="^S",x=red_cigar) & !grepl(pattern="S$",x=red_cigar) & strand=="-",][["clipped_3p"]] <-
 as.character(reverseComplement(unlist(extractAt(x=DNAStringSet(tmpDT[,seq]),at=cigarRangesAlongQuerySpace(cigar=tmpDT[,cigar],ops="S")))))
}
#= reverse  xS
cat( "   >> 5' clipped only\n",sep="" )
tmpDT<-xDT[ !grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="-",]
if ( nrow(tmpDT) != 0 ) {
       xDT[ !grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="-",][["clipped_5p"]] <-
 as.character(reverseComplement(unlist(extractAt(x=DNAStringSet(tmpDT[,seq]),at=cigarRangesAlongQuerySpace(cigar=tmpDT[,cigar],ops="S")))))
}
#= reverse SxS
cat( "   >> both ends clipped\n",sep="" )
tmpDT<-xDT[  grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="-",]
if ( nrow(tmpDT) != 0 ) {
       xDT[  grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="-",][["clipped_3p"]] <-
 as.character(reverseComplement(unlist(extractAt(x=DNAStringSet(tmpDT[,seq]),at=cigarRangesAlongQuerySpace(cigar=tmpDT[,cigar],ops="S")))))[seq(from=1,to=(2*nrow(tmpDT)),by=2)]
}
tmpDT<-xDT[  grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="-",]
if ( nrow(tmpDT) != 0 ) {
       xDT[  grepl(pattern="^S",x=red_cigar) &  grepl(pattern="S$",x=red_cigar) & strand=="-",][["clipped_5p"]] <-
 as.character(reverseComplement(unlist(extractAt(x=DNAStringSet(tmpDT[,seq]),at=cigarRangesAlongQuerySpace(cigar=tmpDT[,cigar],ops="S")))))[seq(from=2,to=(2*nrow(tmpDT)),by=2)]
}

#Sranges <- cigarRangesAlongQuerySpace(cigar = cigar(xBAM), ops = "S")
#Sseqs   <- extractAt(x = xBAM@elementMetadata[,"seq"], at = Sranges)

cat( " ++ Format output table\n",sep="" )
if ( collapsed ) {
 xDT <- xDT[ , c("qname","rCount","cigar","aln.width","chromosome","strand","coord5p","clipped_5p","coord3p","clipped_3p","seq","NH","HI") , with = F ]
# xDT <- xDT[ , c("qname","rCount","chromosome","coord5p","clipped_5p","coord3p","clipped_3p","NH") , with = F ]
} else           {
 xDT[["rCount"]] <- as.integer("1")
 xDT <- xDT[ , c("qname","rCount","cigar","aln.width","chromosome","strand","coord5p","clipped_5p","coord3p","clipped_3p","seq","NH","HI") , with = F ]
# xDT <- xDT[ , c("qname,"chromosome"","coord5p","clipped_5p","coord3p","clipped_3p","NH") , with = F ]
}

cat( " ++ All done\n\n",sep="" )
return(xDT)

}


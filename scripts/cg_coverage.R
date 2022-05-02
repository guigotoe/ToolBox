library(tidyverse)
library(GenomicAlignments)
library(GenomicFeatures)
library(genomation)
library(Gviz)
library(Biostrings)
library(rtracklayer)
library(ggbio)
library(grid)
library(plotly)
library(gridExtra)

path <- '/Users/guillermotorres/Dropbox/Presentations/Courses/2022_tracingDiseases'
system(paste0('cat ',path,'/data/NC_0031* > ',path,'/data/Y_pestis.gff3'))

## NC_003131.gff3 = pCD1 = 70305
## NC_003132.gff3 = pPCP1 = 9612 
## NC_003134.gff3 = pMT1 = 96210
## NC_003143.gff3 = Y_pestis_genome =  4653728 
ypall <- makeTxDbFromGFF(paste0(path,'/data/Y_pestis.gff3'))
ypseq <- readDNAStringSet(paste0(path,'/data/Y_pestis_genome.fna'))
names(ypseq)<- lapply(names(ypseq),function(x) x%>%str_split(' ')%>%unlist()%>%.[[1]])%>%unlist()
isCircular(ypseq) <- c(TRUE,TRUE,TRUE,TRUE)
genome(ypseq) <- c('Y.pestis','Y.pestis','Y.pestis','Y.pestis')
tx <- transcripts(ypall)
isCircular(tx) <- c(TRUE,TRUE,TRUE,TRUE)
seqlengths(tx) <- c(70305,9612,96210,4653728)
genome(tx) <- c('Y.pestis','Y.pestis','Y.pestis','Y.pestis')
elementMetadata(tx)$id <- elementMetadata(tx)[,2]
elementMetadata(tx)$group <- elementMetadata(tx)[,2]
atrack <- AnnotationTrack(tx, name="CDSs",chromosome=options(ucscChromosomeNames=FALSE))
strack <- SequenceTrack(ypseq, genome='Y.pestis', name="Y.pestis")
bams <- list.files(paste0(path,'/data/'),pattern = ".bam$",full.names = T)
altrackx <- AlignmentsTrack(bams[1],isPaired = FALSE,name='Riga (PostBD)')
altracky <- AlignmentsTrack(bams[2],isPaired = FALSE,name='Laishevo (BD)')
gtrack <- GenomeAxisTrack(add53 = TRUE,add35 = TRUE)
dtrackx <- DataTrack(range = bams[1],genome='Y.pestis',type = "l", name = "Riga (PostBD)", window = -1,chromosome=options(ucscChromosomeNames=FALSE))
dtracky <- DataTrack(range = bams[2],genome='Y.pestis',type = "l", name = "Laishevo (BD)", window = -1,chromosome=options(ucscChromosomeNames=FALSE))

plotTracks(list(gtrack,atrack,dtrackx,dtracky), from=1,to=96210,chromosome="NC_003134.1",featureAnnotation='id',fontcolor.feature = 1,cex.feature = 0.6,shape='box',min.height = 8)

plotTracks(list(gtrack,atrack,altrackx,altracky,strack),from=56000,to=61000,chromosome="NC_003134.1",featureAnnotation='id',fontcolor.feature = 1,cex.feature = 0.6,shape='box',min.height = 8)

plotTracks(list(gtrack,atrack,dtrackx,dtracky), from=1,to=70305,chromosome="NC_003131.1",featureAnnotation='id',fontcolor.feature = 1,cex.feature = 0.6,shape='box',min.height = 8)

plotTracks(list(gtrack,atrack,altrackx,altracky,strack),from=9000,to=11000,chromosome="NC_003131.1",featureAnnotation='id',fontcolor.feature = 1,cex.feature = 0.6,shape='box',min.height = 8)





plotTracks(list(gtrack,atrack,altrackx,strack),from=1,to=9610,chromosome="NC_003132.1",featureAnnotation='id',fontcolor.feature = 1,cex.feature = 0.6,shape='box',min.height = 8)

plotTracks(list(gtrack,atrack,altracky,strack),from=1,to=9610,chromosome="NC_003132.1",featureAnnotation='id',fontcolor.feature = 1,cex.feature = 0.6,shape='box',min.height = 8)

plotTracks(list(gtrack,atrack,altrackx,altracky,strack),from=1,to=9610,chromosome="NC_003132.1",featureAnnotation='id',fontcolor.feature = 1,cex.feature = 0.6,shape='box',min.height = 8)

#plotTracks(list(gtrack,atrack,altrack),from=1,to=seqlengths(tx)[4],chromosome="NC_003143.1",featureAnnotation='id',fontcolor.feature = 1,cex.feature = 0.6,shape='box')
plotTracks(list(gtrack,atrack,dtrackx,dtracky), from=1,to=seqlengths(tx)['NC_003131.1'],chromosome="NC_003131.1",shape='box',min.height = 8)


plotTracks(list(gtrack,atrack,altracky,altrackx,strack),from=5000,to=10000,chromosome="NC_003131.1",featureAnnotation='id',fontcolor.feature = 1,cex.feature = 0.6,min.height = 8)


plotTracks(list(gtrack,atrack,altracky,altrackx,strack),from=20000,to=30000,chromosome="NC_003143.1",featureAnnotation='id',fontcolor.feature = 1,cex.feature = 0.6,min.height = 8)

fragmentseq <- getSeq(ypseq,GRanges(seqnames="NC_003143.1",ranges=IRanges(start=12000,end=17400),strand='+'))



seqlengths(tx)['NC_003143.1']

z <- GRanges('NC_003132.1',IRanges(1,seqlengths(tx)['NC_003143.1']))
xnum <- as.numeric(coverage(readGAlignments(bams[1]))$NC_003143.1[ranges(z)])#Uncompress the coverage
ynum <- as.numeric(coverage(readGAlignments(bams[2]))$NC_003143.1[ranges(z)])#Uncompress the coverage
plot(xnum[20000:100000],type='l',col='blue',lwd=2)
lines(ynum[20000:100000], col="red", lwd=2)
pcp <- tibble(seqname='NC_003143.1',start=1:length(xnum),x=xnum,y=ynum)%>%mutate(logdiff=(((x+1)/(y+1))-1))
dd <- pcp%>%filter(logdiff>5 | logdiff<(-5))
GRanges(seqnames=Rle(dd$seqname,1),ranges=IRanges(dd$start-1,end=dd$start))


xnum <- as.numeric(xcov$NC_003134.1[ranges(z)])#Uncompress the coverage
plot(xnum,type='l',col='blue',lwd=2)
xnum <- as.numeric(xcov$NC_003131.1[ranges(z)])#Uncompress the coverage
plot(xnum,type='l',col='blue',lwd=2)
xnum <- as.numeric(xcov$NC_003143.1[ranges(z)])#Uncompress the coverage
plot(xnum,type='l',col='blue',lwd=2)


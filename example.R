library(ComplexHeatmap)
library(rtracklayer)
library(BSgenome.Crobusta.HT.KY)
library(CrobustaTFs)
library(tfenrichr)
library(peakToGene)
library(dirfns)

# motifs are expressed as position weight matrices
# 
CrobustaMotifs
Matrix(CrobustaMotifs)

seqLogo::seqLogo(Matrix(CrobustaMotifs[[1]]))

# Redundancy
# CrobustaTFs::CrobustaMotifs includes only the most specific motifs for each gene (motifs mapped to the fewest genes)
# But I don't have a good way of filtering multiple motifs for the same TF
tags(CrobustaMotifs)

# extract fields from tags(CrobustaMotifs)
meta <- tfTags(CrobustaMotifs, c('KHID', 'KYID', 'Family_Name', 'Source'))
unique(meta[,"KYID"])

write.table(meta, 'CrobustaMotifs.tsv', sep='\t', quote=F)

# write PWMS in MEME format
writeMeme(CrobustaMotifs, 'CrobustaMotifs.meme')

selex <- CrobustaMotifs[meta[,"Source"]=="ANISEED"]

selex <- meta[,"Source"]=="ANISEED"
homer <- meta[,"Source"]=="HOMER" & !meta[,"KYID"]%in%meta[selex,"KYID"]
#there are some redundancies because HOMER includes human & mouse TFs
names(CrobustaMotifs[homer])

cisbp <- meta[,"Source"]=="CisBP" & !meta[,"KYID"]%in%meta[selex|homer,"KYID"]
sum(cisbp)

unique(meta[cisbp,"KYID"])

motifsub <- CrobustaMotifs[selex|homer|cisbp]

#split motifs by KYID
motifsub <- split(motifsub, tfTags(motifsub,"KYID"))

#select only longest motif for each gene
motifsub <- do.call(PWMatrixList, sapply(motifsub, function(x) x[[which.max(tfWidth(x))]]))

# read genomic features
htky <- getFeatures('HT.Gene.gff3')

#read peakome
peaks <- import('accessomeKY.bed')
overlaps <- lapply(htky, getOverlaps, peaks)

# match score is a function of background NT frequency
peakbg <- letterFrequency(Views(Crobusta, peaks), c("A","C","G","T"))
peakbg <- apply(peakbg,2,sum)/sum(peakbg)

#getMatches doesn't directly calculate p-values for each match
#To save time, it just finds a score cutoff for a motif given an alpha value 
library(TFMPvalue)
ptm <- proc.time()
TFMpv2sc(Matrix(CrobustaMotifs[[1]]), 1e-05, peakbg)
proc.time()-ptm

# I wanted to demonstrate the time difference of calculating a p-value for a 10nt motif versus a 22nt motif, but it just runs forever.
# ptm <- proc.time()
# TFMpv2sc(Matrix(CrobustaMotifs[[2]]), 1e-05, peakbg)
# proc.time()-ptm
# 
matches <- getMatches(
	peaks, Crobusta, 'position', motifsub, p.cutoff=1e-05
)

# count number of matches in each peak
countMatches(matches, peaks)
export(unlist(matches), 'matches.bed')

comb <- read.delim('combMotif.tsv')
mespComb <- read.delim('mespComb.tsv')
handrComb <- read.delim('handrComb.tsv')
timeComb <- read.delim('timeComb.tsv')
tissueComb <- read.delim('mesenchymeComb.tsv')

mat <- cbind(mespComb[,1:4], handrComb[,1:3], timeComb[,1:3], tissueComb[,1,drop=F])
row.names(mat) <- make.unique(name(motifsub)[row.names(comb)])
fam <- tfTags(motifsub,c('KYID','Family_Name'))[row.names(comb),]

dir.pdf('DAmotifs',height=36,width=8)
Heatmap(mat, cluster_columns=F,split=fam[,'Family_Name'],row_title_rot=0)
dev.off()

ct <- read.delim('motifct.tsv')
mesp <- read.delim('mesp.tsv')
handr <- read.delim('handr.tsv')
lacz <- read.delim('time.tsv')
gfp <- read.delim('mesenchyme.tsv')

mat <- cbind(mesp[,1:4], handr[,1:3], lacz[,1:3], gfp[,1,drop=F])
sel <- mesp$padj<.05&!is.na(mesp$padj)
mat <- mat[sel,]
fam <- tfTags(motifsub,c('KYID','Family_Name'))[ct$KYID[sel],]

dir.pdf('DAmotifSites',height=20,width=8)
Heatmap(
	mat, cluster_columns=F, show_row_names = F,
	split=fam[,'Family_Name'],row_title_rot=0, gap=unit(0.01,'npc'),
	height=unit(0.8,'npc'), width=unit(0.6, 'npc')
)
dev.off()

## This script is for generation of mock .NA file for the following data training
setwd("/Users/syniu/ncbi/public/sra/")
library(Rsamtools)
library(reshape2)
library(plyr)
# load bam file
#bamfile <- "Ecoli_alignment.BAM"
#bf <- BamFile(bamfile)
# sortBam(file = bamfile, "sorted.bam")
# indexBam(file = "sorted.bam.bam")
bamfile <- "sorted.bam.bam"
bf <- BamFile(bamfile)
seqlengths(bf)
param <- ScanBamParam(which=GRanges("NC_000913.3",
                                    IRanges(start=1, end=seqlengths(bf))))

p_param <- PileupParam(min_mapq = 15, min_base_quality = 10)
pileup(bf, scanBamParam=param, pileupParam = p_param)
res <- pileup(bf, scanBamParam=param, pileupParam=p_param)
plot(res$nucleotide, res$count ,pch=".", log="y", ylab="count (log scale)")

df <- dcast(res, seqnames+pos+strand ~ nucleotide, value.var="count", fun.aggregate=sum)
df_pos <- df[which(df$strand == "+"), ]
pos_hit <- rowSums(df_pos[,4:7])
df_pos <- cbind(df_pos, pos_hit)
df_pos <- df_pos[,c("pos","pos_hit")]
new_df <- data.frame(c(1:seqlengths(bf)))
colnames(new_df) <- "pos"
zz <- join(new_df, df_pos, type = "left", by = "pos")
zz[is.na(zz)] <- 0
df_pos <- zz

df_neg <- df[which(df$strand == "-"), ]
neg_hit <- rowSums(df_neg[,4:7])
df_neg <- cbind(df_neg, neg_hit)
df_neg <- df_neg[,c("pos","neg_hit")]
new_df <- data.frame(c(1:seqlengths(bf)))
colnames(new_df) <- "pos"
zz <- join(new_df, df_neg, type = "left", by = "pos")
zz[is.na(zz)] <- 0
df_neg <- zz
df_na <- join(df_pos, df_neg, type = "left", by = "pos")
df_na <- df_na[,2:3]
write.table(df_na, file="Ecoli.NA", col.names = FALSE, row.names = FALSE, sep="\t")

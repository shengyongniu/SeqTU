setwd("/Users/syniu/ncbi/public/sra/")
library(Rsubread)
library(Rsamtools)
fastq <- "SRR400619.fastq"
#buildindex(basename = "Ecoli", reference = "GCF_000005845.2_ASM584v2_genomic.fna")
## unique is FALSE by default, nBestLocation is 1 here
align(index="Ecoli", readfile1 = "SRR400619.fastq", input_format="FASTQ", output_file="Ecoli_alignment_nbL8.BAM", nBestLocations = 8)
bam.files <- "Ecoli_alignment_nbL8.BAM"
asSam(bam.files, destination=sub("\\.bam", "", bam.files))
sam.file <- "Ecoli_alignment.sam"
# props <- propmapped(files=bam.files)
# props

# Extract quality scores
qs <- qualityScores(filename="SRR400619.fastq",nreads=100)
# Check dimension of qs
dim(qs)
# Check first few elements of qs with head
head(qs)
boxplot(qs)

## FeatureCounts
## non-strand speciifc (duplicate column one nad col two // and strand specific 
fc <- featureCounts(bam.files, annot.ext = "mocked.GTF", isGTFAnnotationFile=TRUE,
                    GTF.featureType  = "CDS", GTF.attrType="position_ID", countMultiMappingReads = TRUE, strandSpecific = 0)
# See what slots are stored in fc
names(fc)
## Take a look at the featurecounts stats
fc$stat
## Take a look at the dimensions to see the number of genes
dim(fc$counts)
## Take a look at the first 6 lines
head(fc$counts)
head(fc$annotation)

save(fc, "fc.RData")

## Gerneration of mocked GTF file for data traning
genomelength=4641628
baseGTFforward<-data.frame(NC="NC123",source="Refseq",feature="CDS",start=c(1:genomelength), end=c(1:genomelength), annotation1=".", strand="+", annotation2=".",note="positionID")
head(baseGTFforward)
baseGTFreverse<-data.frame(NC="NC123",source="Refseq",feature="CDS",start=c(1:genomelength), end=c(1:genomelength), annotation1=".", strand="-", annotation2=".",note="positionID")
head(baseGTFreverse)
baseGTF <- rbind(baseGTFforward, baseGTFreverse)
write.table(baseGTF, "mocked.GTF", row.names = F, col.names = F, sep="\t", quote = F)


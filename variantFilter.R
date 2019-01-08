### exam and resize the referenc files 
library("GenomicRanges")
library("VariantAnnotation")

#define a function to screen variants in CCLE files based on ActionSeq gene list, 
# which contains official symbols from 148 genes, but only 141 have full coverage of exons. 

AcGenFilter <- function(file) {
    dataF1 <- read.csv(file)
    actGen <- read.csv("~/Documents/BioinfoProjects/CCLE_data/ActionSeq_genes.csv")
    filtered <- subset(dataF1, dataF1$Hugo.Symbol %in% actGen$GeneName)
    return(filtered)
}

actGen <- read.csv("OVCAR8_OVARY Mutations.csv")
dim(actGen)


ActionSeqGenes <- read.csv("ActionSeq_genes.csv")
head(ActionSeqGenes)

SKOV3 <- AcGenFilter("SKOV3_OVARY Mutations.csv")
SKOV3[, 1:12]
SKOV3[, 1:9]





OVCAR8 <- AcGenFilter("OVCAR8_OVARY Mutations.csv")
dim(OVCAR8)

OVCAR5 <- AcGenFilter("OVCAR5_OVARY Mutations.csv")
MDAMB468 <- AcGenFilter("MDAMB468_BREAST Mutations.csv")
MDAMB468[, 1:9]
MDAMB231 <- AcGenFilter("MDAMB231_BREAST Mutations-3.csv")
MDAMB231

HCC38 <- AcGenFilter("HCC38_BREAST Mutations.csv")
HCC38
HCC1086 <- AcGenFilter("HCC1806_BREAST Mutations.csv")
HCC1086

## test a VCF file for finding both symbol and start postion 
vignette('VariantAnnotation')
HCC38vcf <- readVcf("HCC38_S7.hg19.mutect2.vcf.gz", "hg19")
HCC38vcf
header(HCC38vcf)
head(info(HCC38vcf))
geno(HCC38vcf)
head(rowRanges(HCC38vcf), 30)

HCC38
## make the above into a GRanges
gr.HCC38 <- makeGRangesFromDataFrame(HCC38,
                                     keep.extra.columns=FALSE,
                                     ignore.strand=FALSE,
                                     seqinfo=NULL,
                                     seqnames.field=c("Chr"),
                                     start.field="Start.Position",
                                     end.field=c("End.Position"))
#seqnames are not matched up with 
newStyle <- mapSeqlevels(seqlevels(gr.HCC38), "UCSC" )
gr.HCC38.ucsc <- renameSeqlevels(gr.HCC38, newStyle)


sum(countOverlaps(gr.HCC38.ucsc, rowRanges(HCC38vcf)))/length(gr.HCC38.ucsc)

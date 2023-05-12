library(Biostrings)
library(hgu95av2probe)

x2 <- sample(hgu95av2probe$sequence, 100, replace=TRUE)
dna <- DNAStringSet(x2)


matches <- vmatchPattern("C", dna)
df <- DataFrame(
  Seqnames = rep(seq_along(dna), lengths(matches)),
  Position = unlist(start(matches))
)

df$Start <- df$Position - 6L
df$End <- df$Position + 6L
df$Pattern <- DNAStringSet( substr(dna[df$Seqnames], df$Start, df$End) )

df[1:10,]

df <- subset(df, width(Pattern) == 13L)

library(GenomicRanges)
GRanges(df)


si <- Seqinfo(as.character(seq_along(dna)), width(dna))
gr <- GRanges(ranges = vmatchPattern("C", dna)[[1]], seqinfo=si[1:7])
window.size <- 11L
windows <- subset(trim(resize(gr, window.size, fix="center")),
                  width == window.size)
gr$pattern <- dna[windows]



query <- DNAString("GAG")
query2 <- DNAString("GGG")
max.mismatch <- 1
myseq <- DNAString("ACCATTGATTAT")
myseq2 <- DNAString("ACCATTGATTAA")
fwd <- vmatchPattern(query, dna, max.mismatch=max.mismatch)
fwd

### the following line does not work
fwd <- as(fwd, "GRanges")

library(BSgenome)
library(BSgenome.Celegans.UCSC.ce2)

matches <- vmatchPattern(pattern = "TCAG", subject = BSgenome.Celegans.UCSC.ce2, exclude = c("M", "_"))



library(ggplot2)
library(Rsamtools)
library(GenomicAlignments)

bam_dir = "/media/disco/work_folder/Downloads/LACES/interface_vcf/bam/"
bam_dir = "/media/disco/work_folder/Downloads/lace_run_test/picard20/"
vcf_dir = "/media/disco/work_folder/Downloads/vcf_20_orgs/vcf_organoids"
file.info(list.files(bam_dir, pattern = ".bam",  full.names = TRUE))["size"] / 1024 /1024
file.info(list.files(vcf_dir, pattern = ".vcf", full.names = TRUE))["size"] / 1024 /1024

bam_filename = list.files(bam_dir, pattern = ".bam$",  full.names = TRUE)[3]
bam_filename
vcf_filename


#load BAM file with Rsamtools
utils:::format.object_size(file.info(bam_filename)$size, "auto")

bam_file <- BamFile(bam_filename)
bam_file

seqinfo(bam_file)


df <- as.data.frame(seqinfo(bam_file))
head(df)
rownames(df)

sequences_on_interest = c("1",   "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9", "10",
                          "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                          "21", "22", "X", "Y", "MT")

df_selection <- df[sequences_on_interest, , drop = FALSE]
df_selection




ggplot(df_selection, aes(rownames(df_selection), seqlengths)) +
  # Y values are given
  geom_bar(stat='identity') +
  
  # Commas at thousand places
  scale_y_continuous(labels = function(x) format(x, big.mark = ",",scientific = FALSE)) +
  
  # Rotate axes for better readability
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  
  # Don't sort labels on X axis alphabetically
  scale_x_discrete(limits=rownames(df_selection)) +
  
  # Plot and axis titles are useful
  ggtitle('Sequence lenghts per chromosome') +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Chromosome') + 
  ylab('Sequence lenght')

quickBamFlagSummary(bam_file)


aln <- scanBam(BamFile(bam_filename, yieldSize=50))
length(aln)
class(aln)
aln <- aln[[1]]
names(aln)
aln


yieldSize(bam_file) <- 1
open(bam_file)
scanBam(bam_file)[[1]]$seq
scanBam(bam_file)[[1]]$seq
close(bam_file)
yieldSize(bam_file) <- NA



gr <- GRanges(seqnames = "1", ranges = IRanges(start = c(0, 10000), end = c(15000, 20000000)))
params <- ScanBamParam(which = gr, what = scanBamWhat())

aln <- scanBam(bam_file, param = params)
names(aln)

head(aln[[1]]$pos)
aln

countBam(bam_file, param = params)




fls <- list.files(system.file("extdata", package="GenomicAlignments"), recursive=TRUE, pattern="*bam$", full=TRUE)
fls


features <- GRanges(seqnames = c(rep("chr2L", 4), rep("chr2R", 5), rep("chr3L", 2)), ranges = IRanges(c(1000, 3000, 4000, 7000, 2000, 3000, 3600, 4000, 7500, 5000, 5400), width=c(rep(500, 3), 600, 900, 500, 300, 900, 300, 500, 500)), "-", group_id=c(rep("A", 4), rep("B", 5), rep("C", 2)))

features                       

olap <- summarizeOverlaps(features, fls)
olap@assays@data$counts

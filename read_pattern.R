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

GenomicRanges::GRanges(df)

library(GenomicRanges)

si <- Seqinfo(as.character(seq_along(dna)), width(dna))
gr <- GRanges(ranges = vmatchPattern("C", dna)[[1]], seqinfo=si[1:7])
window.size <- 11L
windows <- subset(trim(resize(gr, window.size, fix="center")),
                  width == window.size)
gr$pattern <- dna[windows]



query <- DNAString("AGGAGGT")
query2 <- DNAString("AGGAGAT")
max.mismatch <- 1
myseq <- DNAString("ACCATTGATTAT")
myseq2 <- DNAString("ACCATTGATTAA")
fwd <- vmatchPattern(query, dna, max.mismatch=max.mismatch)
fwd <- as(fwd, "GRanges")

library(BSgenome)
library(BSgenome.Celegans.UCSC.ce2)

matches <- vmatchPattern(pattern = "TCAG", subject = BSgenome.Celegans.UCSC.ce2, exclude = c("M", "_"))

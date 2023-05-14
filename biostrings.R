## ---------------------------------------------------------------------
## A. FROM "query" TO "reference" SPACE
## ---------------------------------------------------------------------

library(GenomicAlignments)

## Load read sequences from a BAM file (they will be returned in a
## GAlignments object):
bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools")
param <- ScanBamParam(what="seq")
(gal <- readGAlignments(bamfile, param=param))
(qseq <- mcols(gal)$seq)  # the read sequences (aka query sequences)

## Lay the query sequences alongside the reference space. This will
## remove the substrings associated with insertions to the reference
## (I operations) and soft clipping (S operations), and will inject new
## substrings (filled with "-") where deletions from the reference (D
## operations) and skipped regions from the reference (N operations)
## occurred during the alignment process:
(qseq_on_ref <- sequenceLayer(qseq, cigar(gal)))

## A typical use case for doing the above is to compute 1 consensus
## sequence per chromosome. The code below shows how this can be done
## in 2 extra steps.

## Step 1: Compute one consensus matrix per chromosome.
(qseq_on_ref_by_chrom <- splitAsList(qseq_on_ref, seqnames(gal)))
(pos_by_chrom <- splitAsList(start(gal), seqnames(gal)))

cm_by_chrom <- lapply(names(pos_by_chrom),
                      function(seqname)
                        consensusMatrix(qseq_on_ref_by_chrom[[seqname]],
                                        as.prob=TRUE,
                                        shift=pos_by_chrom[[seqname]]-1,
                                        width=seqlengths(gal)[[seqname]]))

names(cm_by_chrom) <- names(pos_by_chrom)
cm_by_chrom$seq1[,1:10]

## 'cm_by_chrom' is a list of consensus matrices. Each matrix has 17
## rows (1 per letter in the DNA alphabet) and 1 column per chromosome
## position.

## Step 2: Compute the consensus string from each consensus matrix.
## We'll put "+" in the strings wherever there is no coverage for that
## position, and "N" where there is coverage but no consensus.
(cs_by_chrom <- lapply(cm_by_chrom,
                      function(cm) {
                        ## Because consensusString() doesn't like consensus matrices
                        ## with columns that contain only zeroes (and you will have
                        ## columns like that for chromosome positions that don't
                        ## receive any coverage), we need to "fix" 'cm' first.
                        idx <- colSums(cm) == 0
                        cm["+", idx] <- 1
                        DNAString(consensusString(cm, ambiguityMap="N"))
                      }))


## consensusString() provides some flexibility to let you extract
## the consensus in different ways. See '?consensusString' in the
## Biostrings package for the details.

## Finally, note that the read quality strings can also be used as
## input for sequenceLayer():
param <- ScanBamParam(what="qual")
gal <- readGAlignments(bamfile, param=param)
(qual <- mcols(gal)$qual)  # the read quality strings
(qual_on_ref <- sequenceLayer(qual, cigar(gal)))
## Note that since the "-" letter is a valid quality code, there is
## no way to distinguish it from the "-" letters inserted by
## sequenceLayer().

## ---------------------------------------------------------------------
## B. FROM "query" TO "query-after-soft-clipping" SPACE
## ---------------------------------------------------------------------

## Going from "query" to "query-after-soft-clipping" simply removes
## the substrings associated with soft clipping (S operations):
(qseq <- DNAStringSet(c("AAAGTTCGAAAANNN", "TTACGATTANAANNN", "GGATANGTNTTAANN")))
(cigar <- c("3H10M", "2S7M1S2H", "2D3M1I1M1I3M4S"))
(clipped_qseq <- sequenceLayer(qseq, cigar,
                              from="query", to="query-after-soft-clipping"))

sequenceLayer(clipped_qseq, cigar,
              from="query-after-soft-clipping", to="query")

sequenceLayer(clipped_qseq, cigar,
              from="query-after-soft-clipping", to="query",
              S.letter="-")

## ---------------------------------------------------------------------
## C. BRING QUERY AND REFERENCE SEQUENCES TO THE "pairwise" or
##    "pairwise-dense" SPACE
## ---------------------------------------------------------------------

## Load read sequences from a BAM file:
library(RNAseqData.HNRNPC.bam.chr14)
bamfile <- RNAseqData.HNRNPC.bam.chr14_BAMFILES[1]
param <- ScanBamParam(what="seq",
                      which=GRanges("chr14", IRanges(1, 25000000)))
gal <- readGAlignments(bamfile, param=param)
qseq <- mcols(gal)$seq  # the read sequences (aka query sequences)

## Load the corresponding reference sequences from the appropriate
## BSgenome package (the reads in RNAseqData.HNRNPC.bam.chr14 were
## aligned to hg19):
library(BSgenome.Hsapiens.UCSC.hg19)
rseq <- getSeq(Hsapiens, as(gal, "GRanges"))  # the reference sequences

## Bring 'qseq' and 'rseq' to the "pairwise" space.
## For 'qseq', this will remove the substrings associated with soft
## clipping (S operations) and inject substrings (filled with "-")
## associated with deletions from the reference (D operations) and
## skipped regions from the reference (N operations). For 'rseq', this
## will inject substrings (filled with "-") associated with insertions
## to the reference (I operations).
qseq2 <- sequenceLayer(qseq, cigar(gal),
                       from="query", to="pairwise")
rseq2 <- sequenceLayer(rseq, cigar(gal),
                       from="reference", to="pairwise")

## Sanity check: 'qseq2' and 'rseq2' should have the same shape.
stopifnot(identical(elementNROWS(qseq2), elementNROWS(rseq2)))

## A closer look at reads with insertions and deletions:
cigar_op_table <- cigarOpTable(cigar(gal))
head(cigar_op_table)

I_idx <- which(cigar_op_table[ , "I"] >= 2)  # at least 2 insertions
qseq2[I_idx]
rseq2[I_idx]

D_idx <- which(cigar_op_table[ , "D"] >= 2)  # at least 2 deletions
qseq2[D_idx]
rseq2[D_idx]

## A closer look at reads with skipped regions:
N_idx <- which(cigar_op_table[ , "N"] != 0)
qseq2[N_idx]
rseq2[N_idx]

## A variant of the "pairwise" space is the "pairwise-dense" space.
## In that space, all indels and skipped regions are removed from 'qseq'
## and 'rseq'.
qseq3 <- sequenceLayer(qseq, cigar(gal),
                       from="query", to="pairwise-dense")
rseq3 <- sequenceLayer(rseq, cigar(gal),
                       from="reference", to="pairwise-dense")

## Sanity check: 'qseq3' and 'rseq3' should have the same shape.
stopifnot(identical(elementNROWS(qseq3), elementNROWS(rseq3)))

## Insertions were removed:
qseq3[I_idx]
rseq3[I_idx]

## Deletions were removed:
qseq3[D_idx]
rseq3[D_idx]

## Skipped regions were removed:
qseq3[N_idx]
rseq3[N_idx]

## ---------------------------------------------------------------------
## D. SANITY CHECKS
## ---------------------------------------------------------------------
SPACES <- c("reference",
            "reference-N-regions-removed",
            "query",
            "query-before-hard-clipping",
            "query-after-soft-clipping",
            "pairwise",
            "pairwise-N-regions-removed",
            "pairwise-dense")

cigarWidth <- list(
  function(cigar) cigarWidthAlongReferenceSpace(cigar),
  function(cigar) cigarWidthAlongReferenceSpace(cigar,
                                                N.regions.removed=TRUE),
  function(cigar) cigarWidthAlongQuerySpace(cigar),
  function(cigar) cigarWidthAlongQuerySpace(cigar,
                                            before.hard.clipping=TRUE),
  function(cigar) cigarWidthAlongQuerySpace(cigar,
                                            after.soft.clipping=TRUE),
  function(cigar) cigarWidthAlongPairwiseSpace(cigar),
  function(cigar) cigarWidthAlongPairwiseSpace(cigar,
                                               N.regions.removed=TRUE),
  function(cigar) cigarWidthAlongPairwiseSpace(cigar, dense=TRUE)
)

cigar <- c("3H2S4M1D2M2I1M5N3M6H", "5M1I3M2D4M2S")

seq <- list(
  BStringSet(c(A="AAAA-BBC.....DDD", B="AAAAABBB--CCCC")),
  BStringSet(c(A="AAAA-BBCDDD", B="AAAAABBB--CCCC")),
  BStringSet(c(A="++AAAABBiiCDDD", B="AAAAAiBBBCCCC++")),
  BStringSet(c(A="+++++AAAABBiiCDDD++++++", B="AAAAAiBBBCCCC++")),
  BStringSet(c(A="AAAABBiiCDDD", B="AAAAAiBBBCCCC")),
  BStringSet(c(A="AAAA-BBiiC.....DDD", B="AAAAAiBBB--CCCC")),
  BStringSet(c(A="AAAA-BBiiCDDD", B="AAAAAiBBB--CCCC")),
  BStringSet(c(A="AAAABBCDDD", B="AAAAABBBCCCC"))
)

stopifnot(all(sapply(1:8,
                     function(i) identical(width(seq[[i]]), cigarWidth[[i]](cigar))
)))

sequenceLayer2 <- function(x, cigar, from, to)
  sequenceLayer(x, cigar, from=from, to=to, I.letter="i")

identical_XStringSet <- function(target, current)
{
  ok1 <- identical(class(target), class(current))
  ok2 <- identical(names(target), names(current))
  ok3 <- all(target == current)
  ok1 && ok2 && ok3
}

res <- sapply(1:8, function(i) {
  sapply(1:8, function(j) {
    target <- seq[[j]]
    current <- sequenceLayer2(seq[[i]], cigar,
                              from=SPACES[i], to=SPACES[j])
    identical_XStringSet(target, current)
  })
})
stopifnot(all(res))
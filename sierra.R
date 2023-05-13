extdata_path <- system.file("extdata",package = "Sierra")
reference.file <- paste0(extdata_path,"/Vignette_cellranger_genes_subset.gtf")
bamfile <- c(paste0(extdata_path,"/Vignette_example_TIP_sham.bam"), paste0(extdata_path,"/Vignette_example_TIP_mi.bam") )
whitelist.bc.file <- paste0(extdata_path,"/example_TIP_sham_whitelist_barcodes.tsv")
peak.merge.output.file = paste0(extdata_path, "/TIP_merged_peaks.txt")



importFrom magrittr "%>%"
importFrom foreach "%dopar%"
importFrom Matrix writeMM
import utils


peak.sites.file = peak.merge.output.file 
gtf.file = reference.file
bamfile = bamfile[1]
whitelist.file = whitelist.bc.file[1]
output.dir = count.dirs[1] 
countUMI = TRUE
ncores = 1
chr.names = NULL
filter.chr = FALSE 
gene.symbol.ref = 'gene_name'
CBtag='CB'
UMItag='UB'




lock <- tempfile()
whitelist.bc <- read.table(whitelist.file, stringsAsFactors = FALSE)
whitelist.bc <- whitelist.bc[,1]
n.bcs <- length(whitelist.bc)
message("There are ", n.bcs, " whitelist barcodes.")

n.columns <- n.bcs + 1

# read in gene reference
genes.ref <- make_reference(gtf_file = gtf.file, chr.names = chr.names, filter.chr = filter.chr, gene.symbol.ref = gene.symbol.ref)
chr.names <- as.character(unique(genes.ref$chr))
n.genes <- nrow(genes.ref)

peak.sites <- read.table(peak.sites.file, header = T, sep = "\t", quote = '',
                         stringsAsFactors = FALSE)

# Count the peaks
n.total.sites <- nrow(peak.sites)
message("There are ", n.total.sites, "  sites")
message("Doing counting for each site...")

# Set up multiple workers
system.name <- Sys.info()['sysname']
new_cl <- FALSE
if (system.name == "Windows") {
  new_cl <- TRUE
  cluster <- parallel::makePSOCKcluster(rep("localhost", ncores))
  doParallel::registerDoParallel(cluster)
} else {
  doParallel::registerDoParallel(cores=ncores)
}
#print(chr.names)
mat.to.write <- foreach::foreach(each.chr = chr.names, .combine = 'rbind', .packages=c("magrittr")) %dopar% {
  mat.per.chr <- c()
  message("Processing chr: ", each.chr)
  
  
  for(strand in c(1, -1) ) {
    message(" and strand ", strand)
    isMinusStrand <- if(strand==1) FALSE else TRUE
    
    peak.sites.chr <- dplyr::filter(peak.sites, Chr == each.chr & Strand == strand) %>%
      dplyr::select(Gene, Chr, Fit.start, Fit.end, Strand)
    
    peak.sites.chr$Fit.start <- as.integer(peak.sites.chr$Fit.start)
    peak.sites.chr$Fit.end <- as.integer(peak.sites.chr$Fit.end)
    peak.sites.chr <- dplyr::filter(peak.sites.chr, Fit.start < Fit.end)
    
    # If there are no sites in this range, then just keep going
    if(nrow(peak.sites.chr) == 0) {
      next
    }
    
    isMinusStrand <- if(strand==1) FALSE else TRUE
    which <- GenomicRanges::GRanges(seqnames = each.chr, ranges = IRanges::IRanges(1, max(peak.sites.chr$Fit.end) ))
    
    param <- Rsamtools::ScanBamParam(tag=c(CBtag, UMItag),
                                     which = which,
                                     flag=Rsamtools::scanBamFlag(isMinusStrand=isMinusStrand))
    
    aln <- GenomicAlignments::readGAlignments(bamfile, param=param)
    
    nobarcodes <- which(unlist(is.na(GenomicRanges::mcols(aln)[CBtag])))
    noUMI <- which(unlist(is.na(GenomicRanges::mcols(aln)[UMItag])))
    
    
    to.remove <- dplyr::union(nobarcodes, noUMI)
    if (length(to.remove) > 0) {
      aln <- aln[-to.remove]
    }
    
    whitelist.pos <- which(unlist(GenomicRanges::mcols(aln)[CBtag]) %in% whitelist.bc)
    aln <- aln[whitelist.pos]
    
    # Check that at least one whitelisted barcode was identified. 
    if (length(aln) == 0) {
      next
    }
    
    # For de-duplicating UMIs, let's just remove a random read
    # when there is a duplicate
    if(countUMI) {
      GenomicRanges::mcols(aln)$CB_UB <- paste0(unlist(GenomicRanges::mcols(aln)[CBtag]), 
                                                "_", unlist(GenomicRanges::mcols(aln)[UMItag]))
      uniqUMIs <- which(!duplicated(GenomicRanges::mcols(aln)$CB_UB))
      aln <- aln[uniqUMIs]
    }
    
    aln <- GenomicRanges::split(aln, unlist(GenomicRanges::mcols(aln)[CBtag]))
    
    polyA.GR <- GenomicRanges::GRanges(seqnames = peak.sites.chr$Chr,
                                       IRanges::IRanges(start = peak.sites.chr$Fit.start,
                                                        end = as.integer(peak.sites.chr$Fit.end)))
    n.polyA <- length(polyA.GR)
    barcodes.gene <- names(aln)
    res <- sapply(barcodes.gene, function(x) GenomicRanges::countOverlaps(polyA.GR, aln[[x]]))
    
    # Reorder the columns of the res matrix to match the whitelist barcodes
    res.mat <- matrix(0L, nrow = n.polyA, ncol = n.bcs)
    res.mat[,match(barcodes.gene, whitelist.bc)] <- res
    
    # Return a sparse matrix
    mat.per.strand <- Matrix::Matrix(res.mat, sparse = TRUE)
    #mat.to.write <- matrix(0L, nrow = n.polyA, ncol = n.bcs)
    #mat.to.write[,match(barcodes.gene, whitelist.bc)] <- res
    polyA.ids <- paste0(peak.sites.chr$Gene, ":", peak.sites.chr$Chr, ":", peak.sites.chr$Fit.start,
                        "-", peak.sites.chr$Fit.end, ":", peak.sites.chr$Strand )
    rownames(mat.per.strand) <- polyA.ids
    
    
    # Need to combine the two matrices from each strand
    if(is.null(mat.per.chr)) {
      mat.per.chr <- mat.per.strand
    } else {
      mat.per.chr <- rbind(mat.per.chr, mat.per.strand)
    }
  } # Loop for strand
  
  # Return sparse matrix for each chromosome for combining across all threads
  return(mat.per.chr)
} # Loop for chr

if (new_cl) { ## Shut down cluster if on Windows
  ## stop cluster
  parallel::stopCluster(cluster)
}

if (!dir.exists(output.dir)){
  dir.create(output.dir)
}
Matrix::writeMM(mat.to.write, file = paste0(output.dir, "/matrix.mtx"))
write.table(whitelist.bc, file = paste0(output.dir, "/barcodes.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(mat.to.write), file = paste0(output.dir, "/sitenames.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)

## Compress the output files
R.utils::gzip(paste0(output.dir, "/matrix.mtx"), overwrite = TRUE)
R.utils::gzip(paste0(output.dir, "/barcodes.tsv"), overwrite = TRUE)
R.utils::gzip(paste0(output.dir, "/sitenames.tsv"), overwrite = TRUE)




## ========================
## ========================
library(import)

extdata_path <- system.file("extdata",package = "Sierra")
reference.file <- paste0(extdata_path,"/Vignette_cellranger_genes_subset.gtf")
junctions.file <- paste0(extdata_path,"/Vignette_example_TIP_sham_junctions.bed")
bamfile <- c(paste0(extdata_path,"/Vignette_example_TIP_sham.bam"), paste0(extdata_path,"/Vignette_example_TIP_mi.bam") )
peak.output.file <- c("Vignette_example_TIP_sham_peaks.txt",  "Vignette_example_TIP_MI_peaks.txt")

#FindPeaks(output.file=peak.output.file[1], gtf.file = reference.file, bamfile=bamfile[1], junctions.file=junctions.file)


import::from(magrittr, "%>%")
import::from(foreach, "%dopar%")
import::here(GenomicRanges)


output.file = peak.output.file[1]
gtf.file = reference.file
bamfile = bamfile[1]
junctions.file ==junctions.file
min.jcutoff=50 
min.jcutoff.prop = 0.05
min.cov.cutoff = 500
min.cov.prop = 0.05 
min.peak.cutoff=200 
min.peak.prop = 0.05 
ncores = 1
chr.names = NULL 
filter.chr = FALSE 
gene.symbol.ref = 'gene_name'
fit.method = "NLS"
  
  if(!fit.method %in% c("NLS", "MLE") ) { stop("fit.method needs to be either NLS or MLE. ")}
  
  lock <- tempfile()
  #genes.ref <- read.table(reference.file,
  #                        header = TRUE, sep = ",", stringsAsFactors = FALSE)
  #chr.names <- as.character(unique(genes.ref$chr))
  #genes.ref <- subset(genes.ref, chr %in% chr.names)
  
  ## Read in the gtf file
  genes.ref = Sierra:::make_reference(gtf_file = gtf.file, chr.names = chr.names, filter.chr = filter.chr, gene.symbol.ref = gene.symbol.ref)
  n.genes = nrow(genes.ref)
  message(paste(n.genes, "gene entries to process"))
  
  # Initiate the output file
  write("Gene\tChr\tStrand\tMaxPosition\tFit.max.pos\tFit.start\tFit.end\tmu\tsigma\tk\texon/intron\texon.pos\tLogLik", 
        file = output.file)
  
  ## Read in junctions from either bed file or SJ file from STAR aligner
  # step1 determine file type by extension
  
  idx.star <- grep(pattern="SJ.out.tab(.gz|)", x=junctions.file)
  idx.bed <- grep(pattern=".bed(.gz|)", x=junctions.file)
  
  junctions <- read.table(junctions.file, sep = "\t",header = FALSE)
  
  if (length(idx.star) > 0)
  {
    junctions.GR <- GenomicRanges::GRanges(seqnames = junctions$V1, 
                                           IRanges::IRanges(start =junctions$V2,  end = junctions$V3), 
                                           counts = junctions$V8)
  }    
  else if(length(idx.bed) > 0)
  {
    junctions <- cbind(junctions,
                       reshape2::colsplit(junctions$V11, ",", c("blocks1","blocks2")))
    junctions$start <- junctions$V2+junctions$blocks1
    junctions$end <- junctions$V3-junctions$blocks2
    junctions.GR <- GenomicRanges::GRanges(seqnames = junctions$V1,
                                           IRanges::IRanges(start = junctions$start,
                                                            end = junctions$end), counts = junctions$V5)
  }
  else # unrecognisable format
  {
    warning(paste0("Unrecognisable junction file format.\n",
                   "File must have either bed or SJ.out.tab extension.\n",
                   "Cannot continue to findPeaks\n"))
    return(NULL)
  }
  
  # Set up multiple workers
  system.name <- Sys.info()['sysname']
  new_cl <- FALSE
  if (system.name == "Windows") {
    new_cl <- TRUE
    cluster <- parallel::makePSOCKcluster(rep("localhost", ncores))
    doParallel::registerDoParallel(cluster)
  } else {
    doParallel::registerDoParallel(cores=ncores)
  }
  
  i=1
  #foreach::foreach(i = 1:n.genes, .packages = c("GenomicRanges")) %dopar% {
    gene.name <- genes.ref[i, "Gene"]
    seq.name <- genes.ref[i,"chr"]
    gene.start <- genes.ref[i,"start"]
    gene.end <- genes.ref[i, "end"]
    strand <- genes.ref[i,"strand"]
    
    #message(i, " :", gene.name)
    isMinusStrand <- if(strand==1) FALSE else TRUE
    which <- GenomicRanges::GRanges(seqnames = seq.name, ranges = IRanges::IRanges(gene.start, gene.end))
    param <- Rsamtools::ScanBamParam(which = which,
                                     flag=Rsamtools::scanBamFlag(isMinusStrand=isMinusStrand))
    
    aln <- GenomicAlignments::readGAlignments(bamfile, param=param)
    aln_cov <- GenomicRanges::coverage(aln)[seq.name][[1]]
    
    if (length(aln_cov@values) > 1) { ## check for 0 read coverage for this gene
      
      data <- data.frame(pos = seq(gene.start, gene.end),
                         coverage = S4Vectors::runValue(aln_cov)[S4Vectors::findRun(gene.start:gene.end, aln_cov)])
      
      # Find the junction which overlaps this gene
      j.cutoff <- max(min.jcutoff,min.jcutoff.prop*max(data$coverage))
      hits <- GenomicRanges::findOverlaps(which, junctions.GR)
      this.junctions.GR <- junctions.GR[hits@to]
      this.junctions.GR <- this.junctions.GR[this.junctions.GR$counts > j.cutoff]
      #this.junctions.GR <- IRanges::subset(this.junctions.GR, counts > j.cutoff)
      n.junctions <- length(this.junctions.GR)
      data.no.juncs <- data
      
      ## Filter 
      ## This is pretty slow way to do this filtering,
      ## can definitely improve computationally!
      if(n.junctions > 0) {
        for(i in 1:n.junctions) {
          j.start <- IRanges::start(this.junctions.GR[i])
          j.end <- IRanges::end(this.junctions.GR[i])
          data.no.juncs <- data.no.juncs %>%
            dplyr::filter(pos < j.start | pos > j.end)
        }
      }
      
      ## Find peaks 
      
      totalcov <- sum(data$coverage)
      cutoff <- max(min.cov.cutoff,min.cov.prop*totalcov)
      covsum <- totalcov
      maxpeakval <- max(data$coverage)
      maxpeakcutoff <- max(min.peak.cutoff,min.peak.prop*maxpeakval )
      
      #message("Finding exonic sites...")
      n.points <- nrow(data.no.juncs)
      if(n.points > 0) {
        while(covsum > cutoff) {
          maxpeak <- which.max(data.no.juncs$coverage)
          if(data.no.juncs[maxpeak, "coverage"] < maxpeakcutoff) { break }
          start <- maxpeak - 300
          end <- maxpeak + 299
          
          if(start < 1 ) { start <- 1 }
          if(end > n.points) { end <- n.points }
          
          maxval <- max(data.no.juncs$coverage)
          #message(start, ",", end, ",", maxval)
          #message("length: ", data.no.juncs[end,"pos"] - data.no.juncs[start, "pos"])
          #message("max peak: ", data.no.juncs[maxpeak, "pos"])
          
          fit.data <- data.frame(x = seq(1,end-start+1), y = data.no.juncs[start:end,"coverage"])
          
          #nls.res <- NULL
          #tryCatch({
          #  nls.res <- nls( y ~ k*exp(-1/2*(x-mu)^2/sigma^2),
          #                  start=c(mu=300,sigma=100,k=maxval) , data = fit.data)
          #}, error = function(err) { })
          gaussian.fit <- fit_gaussian(fit.data, maxval, fit.method)
          if(!is.null(gaussian.fit)) {
            
            est_mu <- coef(gaussian.fit)["mu"]
            est_sigma <- abs(coef(gaussian.fit)["sigma"]) # this has to be positie in the Gaussian  
            est_k <- coef(gaussian.fit)["k"]
            fit.loglik <- logLik(gaussian.fit)[1] 
            fitted.peak <- maxpeak - 300 + floor(est_mu)
            from <- fitted.peak - 3*floor(est_sigma)
            to <- fitted.peak + 3*floor(est_sigma)
            
            # Handle the cases where the peak is too close to eithr the start or the end
            if(from < 1) { from = 1 }
            if(to > n.points) { to = n.points}
            
            if(fitted.peak <= 0) {
              peak.pos <- "Negative"
            } else {
              peak.pos <- data.no.juncs[fitted.peak, "pos"]
            }
            
            isGapped <- FALSE
            if(to <= 0) {
              to.pos <- "Negative"
            } else {
              to.pos <- data.no.juncs[to, "pos"]
              # Check if this is a spliced region
              pos.gaps <- diff(data.no.juncs[from:to, "pos"])
              
              if(length(which(pos.gaps > 1) > 0)) {
                isGapped <- TRUE
              }
              
            }
            
            
            if(isGapped) {
              exon.pos <- make_exons(data.no.juncs[from:to, "pos"])
              line <- paste(gene.name, seq.name, strand, data.no.juncs[maxpeak, "pos"],
                            peak.pos,
                            data.no.juncs[from, "pos"],
                            to.pos,
                            est_mu, est_sigma, est_k, 
                            "across-junctions", exon.pos, fit.loglik, sep="\t")
            } else {
              exon.pos <- "NA"
              line <- paste(gene.name, seq.name, strand, data.no.juncs[maxpeak, "pos"],
                            peak.pos,
                            data.no.juncs[from, "pos"],
                            to.pos,
                            est_mu, est_sigma, est_k, 
                            "no-junctions", exon.pos, fit.loglik, sep="\t")
            }
            
            
            #print(line)
            locked <- flock::lock(lock)
            write(line,file=output.file,append=TRUE)
            flock::unlock(locked)
          } else {
            line <- paste(gene.name, seq.name, strand, data.no.juncs[maxpeak, "pos"],
                          "NA", "NA", "NA", "NA", "NA", "NA", "no-junctions", "NA", "NA", sep="\t")
            locked <- flock::lock(lock)
            write(line,file=output.file,append=TRUE)
            flock::unlock(locked)
          }
          
          data.no.juncs[start:end, "coverage"] <- 0
          covsum <- sum(data.no.juncs$coverage)
          #print(covsum)
        }
      }
      
      ## Now let's see if there are any peaks in the introns
      ## test each intron separate
      #message("Finding intronic peaks...")
      
      reduced.junctions <- GenomicRanges::reduce(this.junctions.GR)
      n.rjunctions <- length(reduced.junctions)
      
      if(n.junctions > 0) {
        
        for(i in 1:n.rjunctions) {
          #message(i)
          j.start <- IRanges::start(reduced.junctions[i])
          j.end <- IRanges::end(reduced.junctions[i])
          intron.data <- data %>%
            dplyr::filter(pos > j.start & pos < j.end)
          
          if(nrow(intron.data) == 0) { next }
          maxpeak <- which.max(intron.data$coverage)
          maxval <- intron.data[maxpeak, "coverage"]
          
          #message(maxpeak, "    ", maxval)
          if(maxval < maxpeakcutoff) { next }
          fit.data <- data.frame(x = seq(1,nrow(intron.data)),
                                 y = intron.data[,"coverage"])
          #nls.res <- NULL
          #tryCatch({
          #  nls.res <- nls( y ~ k*exp(-1/2*(x-mu)^2/sigma^2),
          #                  start=c(mu=maxpeak,sigma=100,k=maxval) , data = fit.data)
          #}, error = function(err) { })
          gaussian.fit <- fit_gaussian(fit.data, maxval, fit.method, mu=maxpeak) 
          
          if(!is.null(gaussian.fit)) {
            # residuals <- sum(summary(nls.res)$residuals )
            #v <- summary(nls.res)$parameters[,"Estimate"]
            est_mu <- coef(gaussian.fit)["mu"]
            est_sigma <- abs(coef(gaussian.fit)["sigma"]) 
            est_k <- coef(gaussian.fit)["k"]
            fit.loglik <- logLik(gaussian.fit)[1]
            
            fitted.peak <- floor(est_mu)
            from <- fitted.peak - 3*floor(est_sigma)
            to <- fitted.peak + 3*floor(est_sigma)
            
            # Make sure position is within the intron 
            if(from < 1) { from = 1 }
            this.n.points <- nrow(intron.data)
            if(to > this.n.points) { to = this.n.points}
            
            if(fitted.peak <= 0) {
              peak.pos <- "Negative"
            } else {
              peak.pos <- intron.data[fitted.peak, "pos"]
            }
            
            if(to <= 0) {
              to.pos <- "Negative"
            } else {
              to.pos <- intron.data[to, "pos"]
            }
            
            line=paste(gene.name, seq.name, strand, intron.data[maxpeak, "pos"],
                       peak.pos,
                       intron.data[from, "pos"],
                       to.pos,
                       est_mu, est_sigma, est_k, 
                       "no-junctions", "NA", fit.loglik, sep="\t")
            
            #line=paste(gene.name, seq.name, maxpeak, v[1], v[2], v[3], "junctions", sep=",")
            #print(line)
            locked <- flock::lock(lock)
            write(line,file=output.file,append=TRUE)
            flock::unlock(locked)
          } else {
            line=paste(gene.name, seq.name, strand, intron.data[maxpeak, "pos"],
                       "NA", "NA", "NA", "NA", "NA", "NA", "no-junctions", "NA", "NA", sep="\t")
            locked <- flock::lock(lock)
            write(line,file=output.file,append=TRUE)
            flock::unlock(locked)
          }
        }
        
      }
      
    }
    
  #} # End loop for genes
  
  if (new_cl) { ## Shut down cluster if on Windows
    ## stop cluster
    parallel::stopCluster(cluster)
  }
  
  ## As a final step, read in the peak file, filter, and add Peak IDs
  peak.sites.file = output.file
  peak.sites <- read.table(peak.sites.file, header = T, sep = "\t", quote = '',
                           stringsAsFactors = FALSE)
  
  ## Filter the polyA sites
  n.total.sites <- nrow(peak.sites)
  to.filter <- which(peak.sites$Fit.max.pos == "Negative")
  to.filter <- dplyr::union(to.filter, which(peak.sites$Fit.start == "Negative"))
  to.filter <- dplyr::union(to.filter, which(peak.sites$Fit.end == "Negative"))
  to.filter <- dplyr::union(to.filter, which(is.na(peak.sites$Fit.start)))
  to.filter <- dplyr::union(to.filter, which(is.na(peak.sites$Fit.end)))
  to.filter <- dplyr::union(to.filter,  which(is.na(peak.sites$Fit.max.pos)))
  
  if (length(to.filter) > 0)
    peak.sites <- peak.sites[-to.filter,]
  
  ## Check for any examples of peaks with start before end
  sites.diffs <- as.numeric(peak.sites$Fit.end) - as.numeric(peak.sites$Fit.start)
  to.filter <- which(sites.diffs < 0)
  if (length(to.filter) > 0)
    peak.sites <- peak.sites[-to.filter,]
  
  n.filt.sites <- nrow(peak.sites)
  message("There are ", n.total.sites, " unfiltered sites and ", n.filt.sites, " filtered sites")
  
  ## Add polyA IDs to the table
  polyA.ids <- paste0(peak.sites$Gene, ":", peak.sites$Chr, ":", peak.sites$Fit.start,
                      "-", peak.sites$Fit.end, ":", peak.sites$Strand )
  peak.sites$polyA_ID = polyA.ids
  
  ## Remove any duplicates
  peak.sites %>% dplyr::distinct(polyA_ID, .keep_all = TRUE) -> peak.sites
  n.updated.sites = nrow(peak.sites)
  message("There are ", n.updated.sites, " sites following duplicate removal")
  
  ## re-write the updated table
  write.table(peak.sites, file = output.file, sep="\t", quote = FALSE, row.names = FALSE)

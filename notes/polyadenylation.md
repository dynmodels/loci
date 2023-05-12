---
title: polyadenylation
created: 2023-05-04T07:54:51.125Z
modified: 2023-05-11T10:14:49.251Z
---

# polyadenylation

The expression is proportional to the Multiplicity times the factor
$$f=\frac{c\mathrm{UTR}_{len}}{a\mathrm{UTR}_{len}+c\mathrm{UTR}_{len}}\geq0.$$
$f$ is zero mainly due to IS.




$$\mathbb{E}_n(\frac{e_{n,g,c}}{s_n}) =  \mu_g \sum_jM^{\epsilon}_{c,g,j}.$$


$u$ is the new polyadenylation position in respect to the gene intial position

If $u_{max}$ is intended as the length of the gene or the closer p(A)
$$\epsilon = \begin{cases}
u/u_{max} & u\leq u_{max}\\
0 & \text{degradation}
\end{cases}$$

If $u_{max}$ is the length of the gene till the furthest p(A)
$$\epsilon = \begin{cases}
u/u_{max} & u\leq u_{max}\\
u_{max} & \text{no degradation}\\
0 & \text{degradation}
\end{cases}$$


## parameters for the model

| k=genotype ID | g=gene | S=start     | P=PA pos    | A=APA pos   | stop-codon  | $$\epsilon_{3'}=1-\\\theta(A-P)$$ | $$\epsilon_{full}=\epsilon_{3'}\\\tfrac{A-S}{P-S}$$ |
| ------------- | ------ | ----------- | ----------- | ----------- | ----------- | --------------------------------- | --------------------------------------------------- |
| 1             | P53    | chr1:445566 | chr1:449999 | chr1:448999 | chr1:849999 | abcd                              | abcd                                                |
| 1             | AKT1   | chr1:845566 | chr1:849999 | chr1:848999 |             |                                   |                                                     |
| 2             | AKT1   | chr1:845566 | chr1:849999 | chr1:999999 |             |                                   |                                                     |





per ogni gene
	trova stops st
	trova terms standard ts
	trova terms alt ta
	t = ts U ta
	t_end_max = max(end(t))
	new_gene_end = t_end_max+200
	new_t = find_pattern(patterns, gene)
	for each 
	for each t
		if st in t then
			add there is stop



r biostrings load bam - Cerca con Google
https://www.google.com/search?q=r+biostrings+load+bam&client=firefox-b-lm&ei=3rJcZP2lCs7_7_UP96mk0AQ&ved=0ahUKEwj9jLrC9-z-AhXO_7sIHfcUCUoQ4dUDCA4&uact=5&oq=r+biostrings+load+bam&gs_lcp=Cgxnd3Mtd2l6LXNlcnAQAzIFCCEQoAEyBQghEKABOgoIABBHENYEELADOgoIABCKBRCwAxBDOgUIABCABDoTCC4QgAQQsQMQgwEQxwEQ0QMQCjoLCC4QgAQQxwEQrwE6BwguEIAEEAo6BQguEIAEOg0IABCABBCxAxCDARAKOhkILhCABBDHARCvARCXBRDcBBDeBBDgBBgBOgcIABCKBRBDOgYIABAWEB46BwgAEBMQgAQ6BwguEBMQgAQ6CAgAEBYQHhATSgQIQRgASgUIQBIBMVChD1ieQ2C5RWgFcAF4AIABhQGIAYIOkgEEMTIuN5gBAKABAcgBCcABAdoBBggBEAEYFA&sclient=gws-wiz-serp

Introduction to Bioconductor for Sequence Data
https://bioc.ism.ac.jp/packages/3.7/workflows/vignettes/sequencing/inst/doc/sequencing.html

Bioconductor - GenomicAlignments
https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html

r granges endpoint extended - Cerca con Google
https://www.google.com/search?q=r+granges+endpoint+extended&client=firefox-b-lm&ei=O7ZcZLbRFrqB9u8Pvbue0A0&oq=r+granges+endpoint+extend&gs_lcp=Cgxnd3Mtd2l6LXNlcnAQAxgAMgUIIRCgATIFCCEQoAEyBAghEBU6CggAEEcQ1gQQsAM6BQgAEKIEOgsIABDxBBCJBRCiBEoECEEYAFDpCljVF2CyJmgCcAF4AIABpwGIAdwFkgEDNC4zmAEAoAEByAECwAEB&sclient=gws-wiz-serp

granges r - Cerca con Google
https://www.google.com/search?q=granges+r&client=firefox-b-lm&ei=MrdcZLXtA9Dc7_UPxfiD-A0&ved=0ahUKEwi1zt7S--z-AhVQ7rsIHUX8AN8Q4dUDCA4&uact=5&oq=granges+r&gs_lcp=Cgxnd3Mtd2l6LXNlcnAQAzIFCAAQgAQyBggAEBYQHjIGCAAQFhAeMgYIABAWEB4yCAgAEBYQHhAPMggIABAWEB4QCjIGCAAQFhAeMgYIABAWEB4yBggAEBYQHjIGCAAQFhAeOgoIABBHENYEELADOgoIABCKBRCwAxBDOg8ILhCKBRDIAxCwAxBDGAE6EgguEIoFENQCEMgDELADEEMYAToVCC4QigUQxwEQrwEQyAMQsAMQQxgBOhUILhCKBRDHARDRAxDIAxCwAxBDGAE6BwgAEIoFEEM6BQguEIAEOgcILhCKBRBDOgcIABATEIAESgQIQRgAUJYCWK0IYKoJaAJwAXgAgAFoiAGPApIBAzIuMZgBAKABAcgBFMABAdoBBggBEAEYCA&sclient=gws-wiz-serp

An Introduction to the GenomicRanges Package
https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html

R: Intra range transformations of a Ranges, Views, RangesList,...
https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/IRanges/html/intra-range-methods.html

R: Intra range transformations of a GRanges or GRangesList...
https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/intra-range-methods.html

GenomicRanges HOWTOs - GenomicRangesHOWTOs.pdf
https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf

GenomicRanges HOWTOs - GenomicRangesHOWTOs.pdf
https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf#subsection.2.8

Bioconductor - GenomicFeatures
https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html

Making and Utilizing TxDb Objects
https://bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.html

GenomicFeatures: Conveniently import and query gene models - GenomicFeatures.pdf
https://bioconductor.org/packages/release/bioc/manuals/GenomicFeatures/man/GenomicFeatures.pdf

GenomicRanges HOWTOs - GenomicRangesHOWTOs.pdf
https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf#subsection.2.7

GenomicRanges HOWTOs - GenomicRangesHOWTOs.pdf
https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf#subsection.2.1

GenomicRanges HOWTOs - GenomicRangesHOWTOs.pdf
https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf#subsection.2.12

GenomicRanges HOWTOs - GenomicRangesHOWTOs.pdf
https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf#subsection.2.13

6.1 Operations on genomic intervals with GenomicRanges package | Computational Genomics with R
https://compgenomr.github.io/book/operations-on-genomic-intervals-with-genomicranges-package.html



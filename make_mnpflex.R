#!/bin/env Rscript

options(conflicts.policy=list(warn=FALSE))

library(optparse)
option_list <- list(
    make_option(c("-i", "--input"), help = "BED file to convert"),
    make_option(c("-r", "--ref"), help = "Illumina reference BED"),
    make_option(c("-o", "--out"), help = "Output file to create"),
    make_option(c("-m", "--methylation"), help = "Methylation basecodes to include [default 'm']", default="m")
)
opt_parser <- OptionParser(option_list = option_list)

opts <- parse_args(opt_parser)
print(str(opts))
if (is.null(opts$input)) {
    print_help(opt_parser)
    q(status = 1)
}
if (is.null(opts$ref)) {
    print_help(opt_parser)
    q(status = 1)
}
if (is.null(opts$out)) {
    print_help(opt_parser)
    q(status = 1)
}

BASE_KEEP = unlist(strsplit(opts$methylation, ""))
cat(sprintf("Keeping methylation basecodes: [%s]\n", paste(BASE_KEEP, collapse=",")))


library(readr)
library(GenomicRanges, warn.conflicts=F)
library(plyranges, warn.conflicts=FALSE)

ilmn = read_tsv(opts$ref, col_names=F, progress=T, lazy=T, show_col_types=F)
colnames(ilmn) = c("chr", "start", "end", "ref_basecode", "coverage", "strand", "IlmnID")
ilmn = ilmn[, c("chr", "start", "end", "IlmnID")]


bed = read_tsv(opts$input, col_names = F, progress=T, lazy=T, show_col_types=F)
colnames(bed) = c("chr", "start", "end", "basecode", "score", "strand",
  "compat_start", "compat_end", "color", "coverage", "methylation_percentage",
  "Nmod", "Ncanonical", "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall")
# Filter out basecode=h
bed = bed[bed$basecode %in% BASE_KEEP, ]
bed = bed[, c("chr", "start", "end", "coverage", "methylation_percentage")]

bedR = GRanges(
    seqnames = bed$chr,
    ranges   = IRanges(bed$start, bed$end),
    strand   = "*",
    coverage = bed$coverage,
    methylation_percentage = bed$methylation_percentage
)

ilmnR = GRanges(
    seqnames = ilmn$chr,
    ranges   = IRanges(ilmn$start, ilmn$end),
    strand   = "*",
    IlmnID   = ilmn$IlmnID
)

intersectR <- as.data.frame(find_overlaps(bedR, ilmnR))
intersectR <- intersectR[, c(
    "seqnames", "start", "end",
    "coverage", "methylation_percentage", "IlmnID"
)]
colnames(intersectR)[1] <- "chr"

write.table(intersectR, opts$out, col.names=T, row.names=F, sep="\t", quote=TRUE)

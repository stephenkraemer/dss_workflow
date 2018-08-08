# snakemake script to perform DMR calling with DSS::callDMR
#
# Args (provided by snakemake):
#     input:
#         RDS of a DMLtest result for a pairwise comparison
#     output:
#         BED3+X with DMR intervals and additional metadata columns
#     params:
#         dmr_pvalue
#         dmr_delta
#         dmr_minlen
#         dmr_minCG
#         dmr_distance_merge_bp
#         dmr_pct_significant
#
# params correspond to the similarly named callDMR arguments

library(bsseq)
library(GenomicRanges)
library(DSS)
library(purrr)

print(snakemake)

dml_df = readRDS(snakemake@input[[1]])

# R won't throw error for missing list elements
expected_cols = c('dmr_pvalue', 'dmr_delta', 'dmr_minlen', 'dmr_minCG', 'dmr_distance_merge_bp', 'dmr_pct_significant')
stopifnot(all(sapply(expected_cols, function(x) x %in% names(snakemake@params))))

dmrs = callDMR(dml_df,
               p.threshold = snakemake@params[['dmr_pvalue']],
               delta = snakemake@params[['dmr_delta']],
               minlen = snakemake@params[['dmr_minlen']],
               minCG = snakemake@params[['dmr_minCG']],
               dis.merge = snakemake@params[['dmr_distance_merge_bp']],
               pct.sig = snakemake@params[['dmr_pct_significant']]
               )

dmrs = dmrs[order(dmrs[[1]], dmrs[[2]], dmrs[[3]]), ]

# assert that sorting order is as defined for workflow
existing_chromosomes_in_original_order = keep(snakemake@params[['chromosomes']], function(x) x %in% unique(dmrs[[1]]))
stopifnot(all(unique(dmrs[[1]]) == existing_chromosomes_in_original_order))

# convert to BED
# --------------
# for merged CpG input, the output intervals are [1-based start of first CpG, 1-based start of last CpG]
# calculation for the end position:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's say we have a CpG with 1-based coordinates [301, 302]
# zero based coordinates would be [300, 302]
# correct bed coordinate for the dmr interval is then: 302
# the end position given by DSS would be the start of the last contained CpG, ie 301
# so we need to do DSS_end + 1
# length
# ~~~~~~
# DSS takes the distance between the start of the first and last CpG in the DMR.
# We want to take the distance between the start of the first CpG and the *end* of the last CpG.

dmrs['start'] = dmrs['start'] - 1
dmrs['end'] = dmrs['end'] + 1
dmrs['length'] = dmrs['end'] - dmrs['start']

write.table(dmrs, file = snakemake@output[[1]], sep = "\t",
            row.names = FALSE, col.names = TRUE, quote=FALSE)

## setClass("SnakemakeInput", representation = list(input = "list", output = "list", params = "list"))
## snakemake <- new("SnakemakeInput",
##                  input = list("/home/kraemers/temp/dmr_calling_test/dmls/dmls_hsc_vs_mpp2.rds"),
##                   output = list("/home/kraemers/temp/test_dmrs.rds"),
##                   params = list(
##                       dmr_pvalue = 0.05,
##                       dmr_delta = 0.1,
##                       dmr_minlen=50,
##                       dmr_minCG=3,
##                       dmr_distance_merge_bp=50,
##                       dmr_pct_significant=0.5,
##                       chromosomes = list('1', '10', '2', '3')
##                   ))

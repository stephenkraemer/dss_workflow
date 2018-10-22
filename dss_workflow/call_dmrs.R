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
print('Testing that all DSS parameters are there')
expected_parameters = c('pvalue', 'delta', 'minlen', 'minCG', 'merge_dist', 'pct_sign')
dmr_test_parameters = snakemake@params[['dmr_test_parameters']]
stopifnot(all(sapply(expected_parameters, function(x) x %in% names(dmr_test_parameters))))

print('Calling DMRs')
dmrs = callDMR(dml_df,
               p.threshold = dmr_test_parameters[['pvalue']],
               delta = dmr_test_parameters[['delta']],
               minlen = dmr_test_parameters[['minlen']],
               minCG = dmr_test_parameters[['minCG']],
               dis.merge = dmr_test_parameters[['merge_dist']],
               pct.sig = dmr_test_parameters[['pct_sign']]
               )

if (is.null(dmrs)) {
    print('No DMRs found, creating empty dataframe')
    # get header - retrieve programatically in case the DSS columns change
    # create dataframe instead of just writing out header to BED
    # in case I need an rds/feather/... later
    pseudo_dmrs = callDMR(dml_df[1:10000, ], p.threshold=0.5)
    stopifnot(nrow(pseudo_dmrs) > 0)
    dmrs = data.frame(matrix(nrow=0, ncol=ncol(pseudo_dmrs)))
    colnames(dmrs) = colnames(pseudo_dmrs)
} else {
    print('Sorting DMRs')
    dmrs = dmrs[order(dmrs[[1]], dmrs[[2]], dmrs[[3]]), ]

    print('Assert chromosome order')
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

    print('Convert to BED')
    dmrs['start'] = dmrs['start'] - 1
    dmrs['end'] = dmrs['end'] + 1
    dmrs['length'] = dmrs['end'] - dmrs['start']
}

colnames(dmrs)[1] = '#chrom'

write.table(dmrs, file = snakemake@output[[1]], sep = "\t",
            row.names = FALSE, col.names = TRUE, quote=FALSE)

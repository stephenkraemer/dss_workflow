print(snakemake)
chrom_dfs = lapply(snakemake@input, readRDS)
full_df = do.call(rbind, chrom_dfs)
saveRDS(full_df, snakemake@output[[1]])

## For testing
## ===========
## setClass("SnakemakeInput", representation = list(input = "list", output = "list", params = "list"))
## snakemake <- new("SnakemakeInput",
##                   input = list(
##                       "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/pairwise_dmr_based_segmentation/dmls/cmop_vs_cdp_chr17_dmls.rds",
##                       "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/cohort_results/pairwise_dmr_based_segmentation/dmls/cmop_vs_cdp_chr18_dmls.rds"),
##                   output = list('/home/kraemers/temp/test.rds'),
##                   params = list())


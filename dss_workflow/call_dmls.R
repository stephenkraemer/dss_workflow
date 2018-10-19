# Snakemake script to call DMLs between two groups using DSS
#
# Arguments provided by snakemake:
#     input:
#         group1_calls, group2_calls: list of absolute bed file paths
#         group1_sample_ids, group2_sample_ids: list of sample IDs
#     output:
#         dml_df_rds: absolute file path.
#
# Output:
# RDS of DMLtest result. Contains columns chr and pos. For merged CpGs, pos is the 1-based start coordinate of the CpG. There is no end column.
#
# Notes:
# - expected input file format
#   - must have header column
#   - must have columns named chrom (or similar, see below), start, end, n_meth and n_total
#   - may have additional columns, which will be ignored
#   - start, end, n_meth and n_total will be parsed as integer columns,
#     chrom will be parsed as character
# - The chromosome column is allowed to have any name, including names
#   starting with '#'. The name is determined from the first line of
#   each bed file

library(DSS)
library(GenomicRanges)
library(bsseq)
library(rtracklayer)
library(stringr)
library(readr)

main <- function() {

    print(snakemake)

    group1_calls = snakemake@input[['group1_calls']]
    group1_sample_ids = snakemake@params[['group1_sample_ids']]
    group2_calls = snakemake@input[['group2_calls']]
    group2_sample_ids = snakemake@params[['group2_sample_ids']]
    output_fp <- snakemake@output[[1]]

    print('Creating individual bsseq objects')
    all_bsseq_objs = list()
    for (i in seq(1, length(group1_calls))) {
        all_bsseq_objs[[i]] = bed_to_gr(group1_calls[[i]], group1_sample_ids[[i]])
    }
    for (i in seq_along(group2_calls)) {
        j = i + length(group1_calls)
        all_bsseq_objs[[j]] = bed_to_gr(group2_calls[[i]], group2_sample_ids[[i]])
    }


    print('Merging bsseq objects')
    full_bs = combineList(all_bsseq_objs)
    rm('all_bsseq_objs')

    print("Starting dmlTest")
    dml_df = DMLtest(full_bs,
                     group1 = group1_sample_ids,
                     group2 = group2_sample_ids,
                     equal.disp = FALSE,
                     smoothing = snakemake@params[['dml_test_parameters']][['smooth']],
                     smoothing.span = snakemake@params[['dml_test_parameters']][['smooth_span']])

    print("Saving RDS")
    saveRDS(dml_df, output_fp)

    warnings()
    print("Done")

}

bed_to_gr <- function(fp, pid) {

    cat("Reading bed: ", fp)
    con  <- file(fp)
    first_line  <- readLines(con, n=1)
    close(con)
    chrom_name  <- str_split(first_line, "\t")[[1]][1]

    required_cols  <- list(start = 'i', end = 'i', n_meth = 'i', n_total = 'i')
    required_cols[[chrom_name]] <- "c"

    bed_df <- read_tsv(fp, col_types=do.call(cols_only, required_cols), progress=FALSE)
    gr <- GRanges(seqnames=bed_df[[chrom_name]], ranges=IRanges(bed_df[['start']] + 1, bed_df[['end']]))
    gr

    meth_mat <- matrix(bed_df$n_meth, dim(bed_df)[1], 1)
    cov_mat <- matrix(bed_df$n_total, dim(bed_df)[1], 1)
    bs <- BSseq(gr = gr, M = meth_mat, Cov = cov_mat, sampleNames = pid)

    return(bs)
}

main()

## For testing
## ===========
## setClass("SnakemakeInput", representation = list(
##   input = "list", output = "list", params = "list", wildcards = "list"))
## snakemake <- new("SnakemakeInput",
##                  input = list(
##                    group1_calls = list('/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/results_per_pid_july15/hsc_1/meth/meth_calls/mcalls_hsc_1_CG_16.bed.gz',
##                                      '/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/results_per_pid_july15/hsc_2/meth/meth_calls/mcalls_hsc_2_CG_16.bed.gz',
##                                      '/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/results_per_pid_july15/hsc_3/meth/meth_calls/mcalls_hsc_3_CG_16.bed.gz',
##                                      '/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/results_per_pid_july15/hsc_4/meth/meth_calls/mcalls_hsc_4_CG_16.bed.gz',
##                                      '/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/results_per_pid_july15/hsc_5/meth/meth_calls/mcalls_hsc_5_CG_16.bed.gz'),
##                    group2_calls = list('/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/results_per_pid_july15/mpp2_1/meth/meth_calls/mcalls_mpp2_1_CG_16.bed.gz',
##                                      '/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/results_per_pid_july15/mpp2_2/meth/meth_calls/mcalls_mpp2_2_CG_16.bed.gz',
##                                      '/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox/results_per_pid_july15/mpp2_3/meth/meth_calls/mcalls_mpp2_3_CG_16.bed.gz')
##                  ),
##                  output = list("/home/kraemers/temp/dmls.rds"),
##                  params = list(group2_sample_ids=c("mpp2_1", "mpp2_2", "mpp2_3"),
##                                group1_sample_ids=c("hsc_1", "hsc_2", "hsc_3", "hsc_4", "hsc_5")),
##                  wildcards = list()
## )


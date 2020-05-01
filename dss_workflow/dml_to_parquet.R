library(arrow)

print(snakemake)

dml_df = readRDS(snakemake@input[[0]])
write_parquet(dml_df, snakemake@output[[0]])


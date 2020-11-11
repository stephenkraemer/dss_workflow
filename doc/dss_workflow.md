# # Implementation notes

# DML parquet output: currently not possible from within the workflow
# (see comments on rule). Files are still present in metadata table.
# Use utils.dml_test_results_to_parquet to create parquet files,
# from env with r-arrow (not the workflow env, which is anyway automatically
# created and handled by the workflow)

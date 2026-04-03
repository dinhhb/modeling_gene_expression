DATASETS USED:
------------------------
1. gene_expression_matrix.csv
2. gene_expression.csv
3. clinical.csv

SUMMARY OF THE PIPELINES:
------------------------
1. Create sample-source metadata for batch correction (gene_expression_matrix, gene_expression_matrix -> sample_source)
2. Perform batch-correction (gene_expression_matrix, sample_source -> gene_expression_matrix_corrected)
3. Merge clinical+GEX dataset (clinical, gene_expression_matrix_corrected -> clinical_gex)
4. Feature selection and run model (clinical_gex)

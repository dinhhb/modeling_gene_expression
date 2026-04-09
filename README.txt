DATASETS USED:
------------------------
1. gene_expression_matrix.csv
2. gene_expression.csv
3. clinical.csv

SUMMARY OF THE PIPELINES:
------------------------
1_gex_mat.ipynb : Preprocess GEX data, create gex_mat (gene_expression -> gex, gex_mat)
2_clinical.ipynb : Preprocess Clinical data (clinical -> clinical)
3_sample_source.ipynb : Prepare metadata for batch correction (gex-mat, gex, clinical -> sample_source)
r-code : Perform batch correction (gex_mat, sample_source -> gex_mat_corrected)
4_clinical_gex.ipynb : Merge clinical and GEX dataset (clinical, gex_mat_corrected -> clinical_gex)
5_full_pipeline.ipynb : Feature selection and run model (clinical_gex)

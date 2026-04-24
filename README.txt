DATASETS USED:
------------------------
1. gene_expression_matrix.csv
2. gene_expression.csv
3. clinical.csv

SUMMARY OF THE PIPELINES:
------------------------
1_clinical.ipynb                : Preprocess Clinical data          (clinical -> clinical)
2_clinical_gex.ipynb            : Merge clinical and GEX dataset    (clinical, gex -> clinical_gex)
3_gex_mat.ipynb                 : Prepare data for batch correction (clinical_gex
                                                                    -> gex_mat, sample_source)
batch-correction/main.r         : Perform batch correction          (gex_mat, sample_source 
                                                                    -> gex_mat_corrected)
4_leakage_pipe.ipynb       : Merge dataset, feature selection and run model   
                                                                    (gex_mat_corrected, clinical_gex 
                                                                    -> clinical_gex_corrected)

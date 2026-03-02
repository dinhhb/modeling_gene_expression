import pandas as pd
import numpy as np
import json

TRAIN_SOURCES = ['doi:10.1200/PO.16.00054', 'doi:10.1158/2159-8290.CD-13-0617'] # 'Federica Catalanotti, David B. Solit', 'Eliezer M. Van Allen, Dirk Schadendorf' (45+66=111 patients) - full sequencing
VAL_SOURCES = ['doi:10.3390/cancers12082224', 'doi:10.3390/cancers11081203']    # 'Pauline Blateau, Jerome Solassol', 'Baptiste Louveau, Samia Mourah' (24+53=77 patients) - partial sequencing

def load_clinical_raw(split=None):
    if split == 'validation':
        df_patients = pd.read_csv("dataset/clinical.csv")
        df_patients = df_patients.loc[df_patients['source'].isin(VAL_SOURCES)]
    elif split == 'train':
        df_patients = pd.read_csv("dataset/clinical.csv")
        df_patients = df_patients.loc[df_patients['source'].isin(TRAIN_SOURCES)]
    else:
        raise ValueError("Please specify data split: 'train' or 'validation'")
    # drop irrelevant columns
    df_patients.drop(columns=['id', 'creation_datetime'], inplace=True)
    df_patients.replace(to_replace=['<NA>', 'nan', 'N.E.', None], value=np.nan, inplace=True)
    df_patients.reset_index(drop=True, inplace=True)
    return df_patients

# Split drug categories
def split_drug_categories(split=None):
    """ Split drug categories to mono and combined drug therapy categories. """
    
    drug_pats = load_clinical_raw(split=split)[['patientID', 'drug']].copy()
    unique_drugs = ['vemurafenib', 'cobimetinib', 'dabrafenib', 'trametinib']
    drug_cats = dict()
    for _, (id, pat_drugs) in drug_pats.iterrows():
        if pat_drugs is np.nan:
            # pass
            drug_cats[id] = pd.Series([0 for _ in unique_drugs], index=unique_drugs)
        elif ' + ' in pat_drugs:
            drugs = pat_drugs.split(' + ')
            drug_cats[id] = pd.Series([1 for _ in drugs], index=drugs)
        else:
            drug_cats[id] = pd.Series([1], index=[pat_drugs])
    
    return pd.DataFrame(drug_cats).T.fillna(0)

def load_clinical_clean(split=None):
    df_patients = load_clinical_raw(split=split)
    # Drop brain_metastasis, immunotherapy_treatment, M_stage (many missing values)
    df_patients.drop(columns=['M_stage', 'brain_metastasis', 'immunotherapy_treatment', 'OS_status', 'OS_month'], inplace=True)
    # Remove missing values in PFS_status
    df_patients.dropna(subset=['PFS_status'], inplace=True) # (2 rows)
    # change str dtype events to numeric 
    df_patients['PFS_status'] = pd.to_numeric(df_patients.PFS_status)
    # remove patients with only post treatment snp data (6 rows)
    # df_patients = df_patients[df_patients.patientID.isin(['Pat_02', 'Pat_21', 'Pat_27', 'Pat_28', 'Pat_36', 'Pat_37', 'Pat_49']) == False]
    # drop feature AJCC_stage - high class imbalance
    df_patients.drop(columns=['AJCC_stage'], inplace=True)
    # drop few other class
    df_patients.drop(columns=['original_patientID', 'CNA_data', 'SNV_data', 'GEX_data', 'source'], inplace=True)
    # split drug categories
    drug_table = split_drug_categories(split=split)
    df_patients = df_patients.join(drug_table, on='patientID')
    df_patients.drop(columns=['drug'], inplace=True)
    # clean dataset
    df_patients.reset_index(drop=True, inplace=True)
    return df_patients

def load_snps(split = None):
    if split == 'validation':
        df_snps = pd.read_csv("dataset/snvs.csv")
        df_snps = df_snps.loc[df_snps['source'].isin(VAL_SOURCES)]
    elif split == 'train':
        df_snps = pd.read_csv("dataset/snvs.csv")
        df_snps = df_snps.loc[df_snps['source'].isin(TRAIN_SOURCES)]
    else:
        raise ValueError("Please specify data split: 'train' or 'validation'")
    
    # get list of Catalanotti's sequenced genes
    bait = pd.ExcelFile("datafiles/catalanotti_supplement2.xlsx").parse(0)
    list_genes = list(bait['Gene Symbol'])
    df_snps = df_snps[df_snps.HGNC.isin(list_genes)]
    
    # drop irrelevant columns
    # df_snps.drop(columns=['id', 'creation_datetime', 'other_prelevements'], inplace=True)
    df_snps.drop(columns=['id', 'creation_datetime', 'HGVSp', 'other_prelevements', 'uniprot_id'], inplace=True)
    
    # patients with pre-treatment snp data
    df_snps = df_snps[df_snps.temporality == 'pre treatment']
    df_snps.drop(columns=['temporality'], inplace=True)
    
    # fill missing values
    snp_fills = {'consequence' : 'missing', 'variant_classification' : 'missing'}
    df_snps.fillna(value=snp_fills, inplace=True)
    df_snps.drop_duplicates(inplace=True)
    df_snps.reset_index(drop=True, inplace=True)
    return df_snps

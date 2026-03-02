### Module to create and load different datasets
## additional features are added to clinical dataset



# import required packages
import pandas as pd
import numpy as np
from load_data import load_clinical_clean
import integrate_pathway_data as ipd
import os

TRAIN_SOURCES = ['doi:10.1200/PO.16.00054', 'doi:10.1158/2159-8290.CD-13-0617'] # 'Federica Catalanotti, David B. Solit', 'Eliezer M. Van Allen, Dirk Schadendorf' (45+66=111 patients) - full sequencing
VAL_SOURCES = ['doi:10.3390/cancers12082224', 'doi:10.3390/cancers11081203']    # 'Pauline Blateau, Jerome Solassol', 'Baptiste Louveau, Samia Mourah' (24+53=77 patients) - partial sequencing

clinical_data = load_clinical_clean(split='train')

full_genes="CAT"

if full_genes == 'ALL':
    # Melanoma network (KEGG)
    binary_kegg = ipd.map_melanoma_kegg(split='train')
    # Melanoma network (Virtual cell)
    binary_vcells = ipd.map_melanoma_vcells(split='train')
    # Pathway pvalues
    pathway_pvalues = ipd.get_pathway_pvalues(split='train')
elif full_genes == 'CAT':
    #### For Catalanotti Analysis (CAT)
    # Melanoma network (KEGG)
    binary_kegg = ipd.map_melanoma_kegg('CAT_Binary_kegg.csv', split='train')
    # Melanoma network (Virtual cell)
    binary_vcells = ipd.map_melanoma_vcells('CAT_Binary_vcells.csv', split='train')
    # Pathway pvalues
    pathway_pvalues = ipd.get_pathway_pvalues('CAT_Pathway_pvalues.csv', split='train')
elif full_genes == 'CA':
    # #### For Cscape Analysis (CA)
    # Melanoma network (KEGG)
    binary_kegg = ipd.map_melanoma_kegg('CA_Binary_kegg.csv', split='train')
    # Melanoma network (Virtual cell)
    binary_vcells = ipd.map_melanoma_vcells('CA_Binary_vcells.csv', split='train')
    # Pathway pvalues
    pathway_pvalues = ipd.get_pathway_pvalues('CA_Pathway_pvalues.csv', split='train')
    
    
########## MELANOMA NETWORK (KEGG) ##########

# join clinical and network data on patient id
clinical_kegg = clinical_data.join(binary_kegg, on='patientID')
# # Drop empty columns
# #clinical_kegg.dropna(how='all', axis=1, inplace=True)
# drop_kegg = threshold_level(clinical_kegg, binary_kegg.columns)
# clinical_kegg.drop(columns=drop_kegg, inplace=True)
# # fill na with 0
# cols_added = list(set(clinical_kegg) - set(clinical_data))
clinical_kegg[binary_kegg.columns] = clinical_kegg[binary_kegg.columns].fillna(value=0)

# join clinical and network data on patient id
clinical_vcells = clinical_data.join(binary_vcells, on='patientID')
# # Drop empty columns
# #clinical_vcells.dropna(how='all', axis=1, inplace=True)
# drop_vcells = threshold_level(clinical_vcells, binary_vcells.columns)
# clinical_vcells.drop(columns=drop_vcells, inplace=True)
# # fill na with 0
# cols_added = list(set(clinical_vcells) - set(clinical_data))
clinical_vcells[binary_vcells.columns] = clinical_vcells[binary_vcells.columns].fillna(value=0)

########## PATHWAYS P-VALUE ##########

# join clinical and pathways pvalues on patient id
clinical_pp = clinical_data.join(pathway_pvalues, on='patientID')
# # fill na with 0
# cols_added = list(set(clinical_pp) - set(clinical_data))
clinical_pp[pathway_pvalues.columns] = clinical_pp[pathway_pvalues.columns].fillna(value=0)

########## MELANOMA PATHWAY P-VALUE (KEGG) ##########

# add pvalues of melanoma pathway from kegg
clinical_pp_kegg = clinical_data.join(pathway_pvalues['Melanoma'], on='patientID')
clinical_pp_kegg['Melanoma'] = clinical_pp_kegg['Melanoma'].fillna(value=0)

########## MELANOMA PATHWAY P-VALUE (VCELLS) ##########

# add pvalues of melanoma pathway from vcells
clinical_pp_vcells = clinical_data.join(pathway_pvalues['Melanoma Map (Curated)'], on='patientID')
clinical_pp_vcells['Melanoma Map (Curated)'] = clinical_pp_vcells['Melanoma Map (Curated)'].fillna(value=0)



########## COMBINED DATASET ##########

## Combine all additional features
# common genes between kegg and vcells network
common_genes = np.intersect1d(binary_kegg.columns, binary_vcells.columns)
# add kegg network data
combined_c_k = clinical_data.join(binary_kegg, on='patientID')
# add vcells network data ( dropped common genes )
combined_c_k_v = combined_c_k.join(binary_vcells.drop(columns=common_genes), on='patientID')
# add pathway pvalues
combined_c_k_v_p = combined_c_k_v.join(pathway_pvalues, on='patientID')
# # Drop empty columns
# #combined_c_k_v_p.dropna(how='all', axis=1, inplace=True)
# drop_cols = threshold_level(combined_c_k_v_p, combined_c_k_v_p.drop(columns=clinical_data.columns).columns)
# combined_c_k_v_p.drop(columns=drop_cols, inplace=True)
# fill na with 0
cols_added = list(set(combined_c_k_v_p) - set(clinical_data))
combined_c_k_v_p[cols_added] = combined_c_k_v_p[cols_added].fillna(value=0)

#############################################################
#################### VALIDATION DATASETS ####################
#############################################################


# ########## CLINICAL DATASET ##########

# cleaned clinical data
clinical_data_validation = load_clinical_clean(split='validation')
#clinical_data = load_clinical_clean_validation(binary_BRAF=False, censor_val=6.0, drop_dcr=True)

# spns data
binary_kegg_validation = ipd.map_melanoma_kegg(split='validation')
binary_vcells_validation = ipd.map_melanoma_vcells(split='validation')
pathway_pvalues_validation = ipd.compute_pathway_pvalues(split='validation')

########## MELANOMA NETWORK (KEGG) ##########
# join clinical and network data on patient id
clinical_kegg_validation = clinical_data_validation.join(binary_kegg_validation, on='patientID')
clinical_kegg_validation[binary_kegg_validation.columns] = clinical_kegg_validation[binary_kegg_validation.columns].fillna(value=0)


# ########## MELANOMA NETWORK (VCELLS) ##########
# # join clinical and network data on patient id
clinical_vcells_validation = clinical_data_validation.join(binary_vcells_validation, on='patientID')
clinical_vcells_validation[binary_vcells_validation.columns] = clinical_vcells_validation[binary_vcells_validation.columns].fillna(value=0)


# ########## PATHWAYS P-VALUE ##########
# # join clinical and pathways pvalues on patient id
clinical_pp_validation = clinical_data_validation.join(pathway_pvalues_validation, on='patientID')
clinical_pp_validation[pathway_pvalues_validation.columns] = clinical_pp_validation[pathway_pvalues_validation.columns].fillna(value=0)

########## MELANOMA PATHWAY P-VALUE (KEGG) ##########
# add pvalues of melanoma pathway from kegg
clinical_pp_kegg_validation = clinical_data_validation.join(pathway_pvalues_validation['Melanoma'], on='patientID')
clinical_pp_kegg_validation['Melanoma'] = clinical_pp_kegg_validation['Melanoma'].fillna(value=0)


########## MELANOMA PATHWAY P-VALUE (VCELLS) ##########

# add pvalues of melanoma pathway from vcells
clinical_pp_vcells_validation = clinical_data_validation.join(pathway_pvalues_validation['Melanoma Map (Curated)'], on='patientID')
clinical_pp_vcells_validation['Melanoma Map (Curated)'] = clinical_pp_vcells_validation['Melanoma Map (Curated)'].fillna(value=0)


########## LOAD ALL DATASETS ##########

def load_datasets(split=None):
    """ Load all datasets. """
    if split=='validation':
        datasets = {
        'clinical_val' : clinical_data_validation,
        'clinical_kegg_val' : clinical_kegg_validation,
        'clinical_vcells_val' : clinical_vcells_validation,
        'clinical_pp_val' : clinical_pp_validation,
        'clinical_pp_kegg_val' : clinical_pp_kegg_validation,
        'clinical_pp_vcells_val' : clinical_pp_vcells_validation,
    }
    elif split=='train':
        datasets = {
            'clinical' : clinical_data,
            'clinical_kegg' : clinical_kegg,
            'clinical_vcells' : clinical_vcells,
            'clinical_pp' : clinical_pp,
            'clinical_pp_kegg' : clinical_pp_kegg,
            'clinical_pp_vcells' : clinical_pp_vcells,
            'combined_data' : combined_c_k_v_p,
            # 'clinical_fathmm' : clinical_fathmm,
            # 'clinical_cscape' : clinical_cscape,
            # 'clinical_cscape_high' : clinical_cscape_high
        }
    else:
        raise ValueError("Please specify data split: 'train' or 'validation'")
    return datasets

os.makedirs('preprocessed_datasets', exist_ok=True)
    
# Save train datasets
train_datasets = load_datasets(split='train')
for name, df in train_datasets.items():
    path = f'preprocessed_datasets/{name}.csv'
    df.to_csv(path, index=True)
    print(f"Saved: {path} — shape: {df.shape}")

# Save validation datasets
val_datasets = load_datasets(split='validation')
for name, df in val_datasets.items():
    path = f'preprocessed_datasets/{name}.csv'
    df.to_csv(path, index=True)
    print(f"Saved: {path} — shape: {df.shape}")

print("\nAll datasets saved to 'preprocessed_datasets/' folder.")
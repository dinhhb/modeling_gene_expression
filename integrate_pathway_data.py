from scipy.stats import hypergeom
import os
import numpy as np
import pandas as pd
from load_data import load_snps

############ CHECK REQUIRED FILES ############

# datafiles directory
DATAFILES = 'datafiles/'
# files with extracted pathways information
REQUIRED_FILES = [
    'Melanoma_vcells.csv',
    'Melanoma_kegg.csv',
    'Pathways_kegg.csv'
]

TRAIN_SOURCES = ['doi:10.1200/PO.16.00054', 'doi:10.1158/2159-8290.CD-13-0617'] # 'Federica Catalanotti, David B. Solit', 'Eliezer M. Van Allen, Dirk Schadendorf' (45+66=111 patients) - full sequencing
VAL_SOURCES = ['doi:10.3390/cancers12082224', 'doi:10.3390/cancers11081203']    # 'Pauline Blateau, Jerome Solassol', 'Baptiste Louveau, Samia Mourah' (24+53=77 patients) - partial sequencing

# files present in the directory
PRESENT_FILES = os.listdir(DATAFILES)
# check for required files in the directory
if len(np.intersect1d(REQUIRED_FILES, PRESENT_FILES)) < len(REQUIRED_FILES):
    print('Not enough required files')
    # print('Extracting data from resources.')
    # from extract_pathway_data import extract_data
    # extract_data()
    # print(f'Extracted data stored in "{DATAFILES}" directory')

############ MUTATION DATA ############

# function to get list of all mutated genes in patient from clinical and snps data
# used to integrate network/pathways data with patients data

def patients_mut_genes(split = None):
    """ Returns mapping of all mutated genes to respective patients. """
    snps_data = load_snps(split=split).copy()
    
    # if (split=='validation'): # Load mut genes for Blateau & Louveau patients
    #     snps_data = snps_data.loc[snps_data['source'].isin(VAL_SOURCES)]
    # elif split == 'train':
    #     snps_data = snps_data[snps_data.source != 'doi:10.3390/cancers12082224'] # Pauline Blateau, Jerome Solassol
    # else:
    #     raise ValueError("Please specify data split: 'train' or 'validation'")
    
    snps_data = snps_data[['patientID', 'HGNC']]

    return snps_data # cscape_genes_df

############ MELANOMA NETWORK ############

# create binary table of mutated genes respective to patients
def mut_genes_table(mutation_data, network_genes):
    """
    Function to create binary table of mutated genes in patient.

    arguments:
    df_snps : dataframe with patient id, HGNC symbols and mutation status
    network_genes : genes common between the patients and preferred pathway

    returns:
    pd.DataFrame : binary table of mutated genes

    """
    # get list of mutated genes for each patient
    mutated_df = mutation_data[mutation_data['HGNC'].isin(network_genes)].copy(deep=True)
    mutated_df.drop_duplicates(inplace=True)
    mutated_df = mutated_df.groupby(['patientID'], as_index=False).agg({'HGNC': lambda x: x.tolist()})
    # create binary table of mutated genes for patients
    mutation_dict = dict()
    for _, (id, genes) in mutated_df.iterrows():
        if id not in mutation_dict:
            mutation_dict[id] = pd.Series([1 for i in genes], index=genes)
    binary_table = pd.DataFrame(mutation_dict).T#.fillna(0)
    binary_table = pd.DataFrame(binary_table, columns=network_genes)
    binary_table = binary_table[sorted(binary_table.columns)]
    return binary_table


# map patients snp data to melanoma network from kegg
def map_melanoma_kegg(filename = None, split = None):
    mutation_data = patients_mut_genes(split=split).copy()

    if split=='validation':
        filename = 'Binary_kegg_validation.csv'
    elif split == 'train':
        filename = filename if filename is not None else 'Binary_kegg.csv'
    else:
        raise ValueError("Please specify data split: 'train' or 'validation'")

    if filename in PRESENT_FILES:
        return pd.read_csv(DATAFILES+filename, index_col=0)
    else:
        melanoma_kegg = pd.read_csv(DATAFILES+'Melanoma_kegg.csv')
        nodes_kegg = pd.concat([melanoma_kegg['src'], melanoma_kegg['dest']]).unique()  # unique values from src+dest (70 KEGG proteins) = ['AKT3', 'BRAF', etc]
        binary_kegg = mut_genes_table(mutation_data, nodes_kegg)    # like one hot encoding: patient - 70 proteins
        binary_kegg.to_csv(DATAFILES+filename, index=True)
        return binary_kegg
    
# map patients snp data to curated melanoma network from virtual cell
def map_melanoma_vcells(filename = None, split=None):
    mutation_data = patients_mut_genes(split=split).copy()

    if split=='validation':
        filename = 'Binary_vcells_validation.csv'
    elif split == 'train':
        filename = filename if not None else 'Binary_vcells.csv'
    else:
        raise ValueError("Please specify data split: 'train' or 'validation'")

    if filename in PRESENT_FILES:
        return pd.read_csv(DATAFILES+filename, index_col=0)
    else:
        melanoma_vcells = pd.read_csv(DATAFILES+'Melanoma_vcells.csv')
        nodes_vcells = pd.concat([melanoma_vcells['node1'], melanoma_vcells['node2']]).unique()
        binary_vcells = mut_genes_table(mutation_data, nodes_vcells)
        binary_vcells.to_csv(DATAFILES+filename, index=True)
        return binary_vcells
    
    
############ HYPERGEOMETRIC TEST ############

# Perform hypergeometric test of mutated genes in all pathways
def perform_hypergeometric_test(patient_mutations, pathway_mapping, genes_population):
    # population of protein coding genes
    M = len(genes_population)

    patients_pvals = {}  # store info
    for i, (patient, pat_genes) in patient_mutations.iterrows():
        # mutated genes in patient
        #pat_genes = df_snps_comp[df_snps_comp.patientID == patient].HGNC.unique()
        n = np.intersect1d(pat_genes, genes_population)
        
        pathway_pvals = {} # store pvals
        for pathway in pathway_mapping.pathway.unique():
            # total no. of genes in the pathway
            path_nodes = pathway_mapping[pathway_mapping.pathway == pathway].nodes.unique()
            N = np.intersect1d(path_nodes, genes_population)
            # no. of genes mutated in pathway
            #path_mut_genes = np.intersect1d(pat_genes, path_nodes)
            k = len(np.intersect1d(n, N))
            # perform hypergeometric test
            pval = hypergeom.pmf(k=k, M=M, n=len(n), N=len(N))
            value = -np.log10(pval)
            #lens = (len(pat_genes), len(n), len(path_nodes), len(N), len(path_mut_genes), len(k))
            #print(lens, pval, value)
            pathway_pvals[pathway] = value
        
        patients_pvals[patient] = pathway_pvals

    hg_results = pd.DataFrame(patients_pvals).T
    #print(hg_results.shape)
    return hg_results


# Calculate the pvalues for all pathways
def compute_pathway_pvalues(split=None):
    """ Perform hypergeometric test on mutated gene for all pathways and return pvalues """

    ## All protein coding genes fron HGNC database
    genes_population = pd.read_csv(DATAFILES+'hgnc_complete_set_2022-09-01.txt', sep='\t')
    genes_population = genes_population[genes_population.locus_group == 'protein-coding gene'].symbol
    #genes_population = genes_population[['symbol', 'alias_symbol', 'prev_symbol',]]

    ## Genes mapped to pathways
    # Load file with genes mapped to all pathways in kegg
    pathways_kegg = pd.read_csv(DATAFILES+'Pathways_kegg.csv')
    # Load file with genes mapped to Vcells melanoma pathway
    edges_vcells = pd.read_csv(DATAFILES+'Melanoma_vcells.csv')
    nodes_vcells = pd.concat([edges_vcells['node1'],edges_vcells['node2']]).unique()
    pathway_vcells = pd.DataFrame({'pathway': 'Melanoma Map (Curated)', 'nodes': nodes_vcells})
    # dataframe with genes mapped to all pathways (kegg + vcells)
    pathway_mapping = pd.concat([pathways_kegg, pathway_vcells], ignore_index=True)
    
    snps = load_snps(split=split)

    # # Patients mutation data
    # if split=='validation':
    #     # only patients with complete sequencing information
    #     snps = snps.loc[snps['source'].isin(VAL_SOURCES)]
    # elif split=='train':
    #     # only patients with complete sequencing information
    #     snps = snps[snps.source != 'doi:10.3390/cancers12082224'] # Pauline Blateau, Jerome Solassol
    # else:
    #     raise ValueError("Please specify data split: 'train' or 'validation'")
    
    patient_mutations = snps[['patientID', 'HGNC']].drop_duplicates()
    patient_mutations = patient_mutations.groupby(['patientID'], as_index=False).agg({'HGNC': lambda x: x.tolist()})

    ## Perform hypergeometric test for all pathways
    test_results = perform_hypergeometric_test(patient_mutations, pathway_mapping, genes_population)
    # test_results = perform_hypergeometric_test(cscape_genes_df, pathway_mapping, genes_population) # Cscape analysis

    ## Scale the values
    from sklearn.preprocessing import MinMaxScaler
    scaled_results = MinMaxScaler().fit_transform(test_results)
    scaled_results = pd.DataFrame(scaled_results, columns=test_results.columns, index=test_results.index)

    return scaled_results


# Load computed pathway pvalues
def get_pathway_pvalues(filename = None, split = None):
    
    if split=='validation':
        filename='Pathway_pvalues_validation.csv'
    elif split=='train':
        filename = filename if not None else 'Pathway_pvalues.csv'
    else:
        raise ValueError("Please specify data split: 'train' or 'validation'")
    
    if filename in PRESENT_FILES:
        return pd.read_csv(DATAFILES+filename, index_col=0)
    else:
        print('computing pvalues ....')
        pval_df = compute_pathway_pvalues(split=split)
        pval_df.to_csv(DATAFILES+filename, index=True)
        return pval_df

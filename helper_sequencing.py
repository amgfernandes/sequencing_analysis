__authors__ = 'fernandes Dec 2020'
#TODO finish helper class. Add metadata possibility

import pandas as pd
import scanpy as sc
def load_samples(data_location=None,load_method=None,samplelist=None,cache=False):
    '''
    """[summary]
    loads list of samples
    Returns:
    Anndata [list]: [adds batch number automatically based on sample order]
    """
    '''
    samples=[]
    n=0
    for sample in samplelist:
        n+=1
        if load_method==sc.read_h5ad:
            data=load_method(data_location + sample +"/adata.h5ad")
            print (n)
            print (data)
            data.obs['batch']=n
            samples.append(data)
            print ('loading with h5ad: check that gene names are recognized')
        if load_method==sc.read_10x_mtx:
            data=load_method(data_location + sample, # the directory with the `.mtx` file
            var_names='gene_ids',
            cache=cache)
            print (n)
            data.obs['batch']=n
            samples.append(data)
        if load_method==sc.read_mtx:
            data=load_method(data_location + sample + '/cells_x_genes.mtx') # the directory with the `.mtx` file
            data.obs.index = pd.read_csv(data_location + sample +  '/cells_x_genes.barcodes.txt', header=None)[0].values
            data.var.index = pd.read_csv(data_location + sample + '/cells_x_genes.genes.txt', header=None)[0].values
            data.var_names_make_unique() #Makes the index unique by appending a number string to each duplicate index element
            print ('batch number:',n)
            print (data)
            data.obs['batch']=n
            samples.append(data)
            print ('loading with read_mtx: check that gene names are recognized') 
    print ('loading data finished')
    return samples

class Sequencing():
    """[summary]
       Class with functions to help sequencing analysis
        """
    def __init__(self, seqdata,metadata=None):
        """ Loads the experiment from sequencing data
        :param seqdata: sequencing data
        :param metadata: metadata for the sequencing
        """
        self.seqdata = seqdata
        self.metadata = metadata

    def return_indices_not_in_list(self,a, b):
        b_set = set(b)
        return [i for i, v in enumerate(a) if v not in b_set]

    def return_indices_inside_list(self,a, b):
        b_set = set(b)
        return [i for i, v in enumerate(a) if v in b_set]

    def remove_gene_list(self,dgc_mat, gene_list):
        '''Function receives a list of gene matrices and a list of genes and returns
        a gene matrix without the list gene rows.'''
        #The loaded 10x data is a dgc matrix.
        dgc_corr_res=[]
        for dgc in dgc_mat:
            idx=self.return_indices_not_in_list(dgc.var_names,gene_list.gene.values)
            dgc_mat_corr= dgc[:,idx]
            dgc_corr_res.append(dgc_mat_corr)
        return (dgc_corr_res)

    def keep_gene_list(self,dgc_mat, gene_list):
        '''Function receives a list of gene matrices and a list of genes and
        returns a gene matrix containing the list of gene rows.'''
        #The loaded 10x data is a dgc matrix.
        dgc_corr_res=[]
        for dgc in dgc_mat:
            idx=self.return_indices_inside_list(dgc.var_names,gene_list.gene.values)
            dgc_mat_corr= dgc[:,idx]
            dgc_corr_res.append(dgc_mat_corr)
        return (dgc_corr_res)

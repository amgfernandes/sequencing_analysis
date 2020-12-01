#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scanpy.external as sce
from matplotlib import rcParams


# In[2]:


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


# In[3]:


results_file = '/home/fernandes/sample_data/scanpy_test.h5ad'  # the file that will store the analysis results


# #### Read in the count matrix into an `AnnData <https://anndata.readthedocs.io/en/latest/anndata.AnnData.html>`__ object, which holds many slots for annotations and different representations of the data. It also comes with its own HDF5 file format: .h5ad.

# In[4]:


IEG_list=pd.read_csv('/home/fernandes/sample_data/IEG_list.csv', header=None)
IEG_list.columns=['gene']
IEG_list.gene.values


# In[5]:


def return_indices_of_a(a, b):
    b_set = set(b)
    return [i for i, v in enumerate(a) if v not in b_set]

def remove_genes(dgc_mat, gene_list):
#Function receives a  gene matrix and a list of genes and returns a  gene matrix without the gene rows.
#The loaded 10x data is a dgc matrix.
  
    idx=return_indices_of_a(dgc_mat.var_names,gene_list.gene.values)
    dgc_mat_corr= dgc_mat[:,idx] 
  
    return (dgc_mat_corr)


# In[262]:


adata1 = sc.read_10x_mtx(
    '/home/fernandes/sample_data/sample_1/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata1.obs['batch']='1'
adata1=remove_genes(adata1, IEG_list)


# In[263]:


adata2 = sc.read_10x_mtx(
    '/home/fernandes/sample_data/sample_2/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata2.obs['batch']='2'
adata2=remove_genes(adata2, IEG_list)


# In[264]:


adata3 = sc.read_10x_mtx(
    '/home/fernandes/sample_data/sample_3/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata3.obs['batch']='3'
adata3=remove_genes(adata3, IEG_list)


# In[265]:


adata4 = sc.read_10x_mtx(
    '/home/fernandes/sample_data/sample_4/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata4.obs['batch']='4'
adata4=remove_genes(adata4, IEG_list)


# ### Sample 5 was low quality. Not used

# In[266]:


adata6 = sc.read_10x_mtx(
    '/home/fernandes/sample_data/sample_6/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata6.obs['batch']='6'
adata6=remove_genes(adata6, IEG_list)


# In[267]:


adata7 = sc.read_10x_mtx(
    '/home/fernandes/sample_data/sample_7/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata7.obs['batch']='7'
adata7=remove_genes(adata7, IEG_list)


# In[268]:


adata8 = sc.read_10x_mtx(
    '/home/fernandes/sample_data/sample_8/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata8.obs['batch']='8'
adata8=remove_genes(adata8, IEG_list)


# In[269]:


adata9 = sc.read_10x_mtx(
    '/home/fernandes/sample_data/sample_9/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata9.obs['batch']='9'
adata9=remove_genes(adata9, IEG_list)
adata9


# In[270]:


adata4.obs['batch']


# In[271]:


adata =adata1.concatenate(adata2,adata3, adata4, adata6, adata7, adata8, adata9)


# In[272]:


'''remove not used objects'''
del (adata1,adata2,adata3, adata4, adata6, adata7, adata8, adata9)


# # Preprocessing
# Show those genes that yield the highest fraction of counts in each single cells, across all cells.

# In[273]:


sc.pl.highest_expr_genes(adata, n_top=20, )


# #### Basic filtering

# In[248]:


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


# In[249]:


adata.var['mt'] = adata.var_names.str.startswith('mt')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[250]:


sc.pl.violin(adata, ['n_genes_by_counts'],
             jitter=0.3, multi_panel=False)


# In[251]:


sc.pl.violin(adata, ['total_counts'],
             jitter=0.3, multi_panel=False)


# In[252]:


sc.pl.violin(adata, ['pct_counts_mt'],
             jitter=0.3, multi_panel=False)


# In[253]:


sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


# In[254]:


adata = adata[adata.obs.n_genes_by_counts < 3000, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]


# ### Total-count normalize (library-size correct) the data matrix 𝐗 to 10,000 reads per cell, so that counts become comparable among cells.

# In[255]:


sc.pp.normalize_total(adata, target_sum=10000)


# ### Logarithmize the data.

# In[256]:


sc.pp.log1p(adata)


# In[257]:


sc.pp.scale(adata, max_value=10)


# ### Identify highly-variable genes.

# In[261]:


sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


# In[259]:


sc.pl.highly_variable_genes(adata)


# ### Set the .raw attribute of AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.

# In[30]:


adata.raw = adata


# In[31]:


# filter
adata = adata[:, adata.var.highly_variable]


# In[32]:


#Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])


# In[33]:


#Scale each gene to unit variance. Clip values exceeding standard deviation 10
sc.pp.scale(adata, max_value=10)


# In[34]:


'''#Harmony works by adjusting the principal components, this function should be run after performing PCA but before computing the neighbor graph,'''
sc.tl.pca(adata,svd_solver='arpack') 


# In[35]:


sc.pl.pca_variance_ratio(adata, log=True)


# In[36]:


adata_corr=adata.copy()


# In[37]:


sce.pp.harmony_integrate(adata_corr, 'batch')#correct batch effect


# In[38]:


'X_pca_harmony' in adata_corr.obsm


# In[39]:


adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
adata_corr.var_names_make_unique()


# In[40]:


sc.pp.neighbors(adata_corr, n_neighbors=10, n_pcs=40)


# In[41]:


sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)


# In[42]:


sc.tl.umap(adata_corr)


# In[43]:


sc.tl.umap(adata)


# In[44]:


#rcParams['figure.figsize'] = 5, 5
sc.pl.umap(adata, color=['batch'])


# In[45]:


#rcParams['figure.figsize'] = 5, 5
sc.pl.umap(adata_corr, color=['batch'])


# In[46]:


adata_corr.write(results_file)


# In[47]:


adata_corr


# In[48]:


adata_corr.var_names


# In[55]:


sc.pl.umap(adata_corr, color=['rplp1', 'stmn1b', 'rps20'])


# In[50]:


sc.tl.leiden(adata_corr)


# In[58]:


rcParams['figure.figsize'] = 10,10
sc.pl.umap(adata_corr, color=['leiden','rplp1', 'stmn1b', 'rps20'])


# In[59]:


rcParams['figure.figsize'] = 5,5
sc.pl.umap(adata_corr, color='leiden')


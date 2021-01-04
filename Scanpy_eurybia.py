#!/usr/bin/env python
# coding: utf-8

# %%:


from anndata._core.aligned_mapping import I
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scanpy.external as sce
from matplotlib import rcParams
from helper_sequencing import Sequencing
from helper_sequencing import load_samples
import itertools as it

# %%:

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


# %%:


#results_file = '/home/fernandes/sample_data/scanpy_test.h5ad'  # the file that will store the analysis results
results_file = '/home/fernandes/tests/scanpy_tests.h5ad'  # the file that will store the analysis results

# #### Read in the count matrix into an `AnnData <https://anndata.readthedocs.io/en/latest/anndata.AnnData.html>`__ object, which holds many slots for annotations and different representations of the data. It also comes with its own HDF5 file format: .h5ad.

# In[4]:


IEG_list=pd.read_csv('/home/fernandes/sample_data/IEG_list.csv', header=None)
IEG_list.columns=['gene']
IEG_list.gene.values

# In[6]:
#data_list=["sample_1","sample_2","sample_3", "sample_4","sample_5","sample_6","sample_7","sample_8","sample_9"]
data_list=["test1"]

#folder_with_data='/home/fernandes/sample_data/'
folder_with_data='/home/fernandes/tests/'

#seq_data=load_samples(data_location=folder_with_data,load_method=sc.read_10x_mtx,samplelist=data_list)
seq_data=load_samples(data_location=folder_with_data,load_method=sc.read_h5ad,samplelist=data_list)
seq_helper=Sequencing(seqdata=seq_data)

seq_data_filtered=seq_helper.remove_gene_list(seq_data,IEG_list)

#%%
# Save figure (set to True to save)
folder_with_data_plot='/home/fernandes/tests/plots'
sc.settings.autosave = True #save figures True/False
sc.settings.figdir = folder_with_data_plot
sc.settings.set_figure_params(dpi_save=320, format="png")

outputDirectory=folder_with_data_plot
if not os.path.isdir(outputDirectory):
      os.mkdir( outputDirectory )

#%%

'''iterate over all elements on list of samples and concatenate them...only keep last result'''
adata=list(it.accumulate(seq_data_filtered, sc.AnnData.concatenate))[-1]

# list(it.accumulate([adata1, adata2, adata2], sc.AnnData.concatenate))[-1]

#%%


'''remove not used objects'''
del (seq_data_filtered,seq_data)


# # Preprocessing
# Show those genes that yield the highest fraction of counts in each single cells, across all cells.
#%%
sc.pl.highest_expr_genes(adata, n_top=20)
# #### Basic filtering

#%%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#%%
adata.var['mt'] = adata.var_names.str.startswith('mt')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

#%%
sc.pl.violin(adata, ['n_genes_by_counts','total_counts','pct_counts_mt'],
             jitter=0.3, multi_panel=True)


#%%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='pct_counts_mt.png')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='n_genes_by_counts.png')


#%%
adata = adata[adata.obs.n_genes_by_counts < 3000, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
# Total-count normalize (library-size correct) the data matrix 𝐗 to 10,000 reads per cell, so that counts become comparable among cells.
#%%
sc.pp.normalize_total(adata, target_sum=10000)



#%%
'''Logarithmize the data.'''
sc.pp.log1p(adata)


#%%

#Scale each gene to unit variance. Clip values exceeding standard deviation 10
sc.pp.scale(adata, max_value=10)
#%%
'''Identify highly-variable genes.'''
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


#%%:
sc.pl.highly_variable_genes(adata)


#%%
'''Set the .raw attribute of AnnData object to the normalized and logarithmized raw gene
 expression for later use in differential testing and visualizations of gene expression.
  This simply freezes the state of the AnnData object.'''

adata.raw = adata
#%%
'''filter only highly variable genes'''
adata = adata[:, adata.var.highly_variable]

#%%

#Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])


#%%

'''#Harmony works by adjusting the principal components, this function should be run after performing PCA but before computing the neighbor graph,'''
sc.tl.pca(adata,svd_solver='arpack') 


#%%

sc.pl.pca_variance_ratio(adata, log=True)


#%%
'''make a copy to use for batch correction'''
adata_corr=adata.copy()


#%%
sce.pp.harmony_integrate(adata_corr, 'batch')#correct batch effect


#%%
'X_pca_harmony' in adata_corr.obsm


#%%
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
adata_corr.var_names_make_unique()

#%%
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


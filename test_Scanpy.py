#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scanpy.external as sce
from matplotlib import rcParams


# In[3]:


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


# In[4]:


results_file = 'D:/scRNA-seq_demo_data/test.h5ad'  # the file that will store the analysis results


# #### Read in2 the count matrix into an `AnnData <https://anndata.readthedocs.io/en/latest/anndata.AnnData.html>`__ object, which holds many slots for annotations and different representations of the data. It also comes with its own HDF5 file format: .h5ad.

# In[5]:


adata1 = sc.read_10x_mtx(
    'D:/scRNA-seq_demo_data/sample_1/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata1.obs['batch']='1'


# In[6]:


adata2 = sc.read_10x_mtx(
    'D:/scRNA-seq_demo_data/sample_2/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata2.obs['batch']='2'


# In[7]:


adata3 = sc.read_10x_mtx(
    'D:/scRNA-seq_demo_data/sample_3/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata3.obs['batch']='3'


# In[8]:


adata4 = sc.read_10x_mtx(
    'D:/scRNA-seq_demo_data/sample_4/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata4.obs['batch']='4'


# In[9]:


adata5 = sc.read_10x_mtx(
    'D:/scRNA-seq_demo_data/sample_5/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata5.obs['batch']='5'


# In[10]:


adata6 = sc.read_10x_mtx(
    'D:/scRNA-seq_demo_data/sample_6/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata6.obs['batch']='6'


# In[11]:


adata7 = sc.read_10x_mtx(
    'D:/scRNA-seq_demo_data/sample_7/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata7.obs['batch']='7'


# In[12]:


adata8 = sc.read_10x_mtx(
    'D:/scRNA-seq_demo_data/sample_8/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading
adata8.obs['batch']='8'


# In[13]:


adata1.obs['batch']


# In[14]:


adata2.obs['batch']


# In[15]:


adata =adata1.concatenate(adata2, adata3, adata4, adata5, adata6, adata7, adata8)


# In[16]:


'''#Harmony works by adjusting the principal components, this function should be run after performing PCA but before computing the neighbor graph,'''
sc.tl.pca(adata) 


# In[17]:


sce.pp.harmony_integrate(adata, 'batch') #correct batch effect


# In[18]:


'X_pca_harmony' in adata.obsm


# In[20]:


adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`


# In[59]:


rcParams['figure.figsize'] = 10, 10
sc.pl.umap( adata, color=['batch'])


# In[21]:


adata


# In[22]:


sc.pl.highest_expr_genes(adata, n_top=20, )


# In[23]:


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


# In[24]:


adata.var['mt'] = adata.var_names.str.startswith('mt')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[25]:


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


# In[26]:


sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


# In[27]:


adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]


# In[28]:


sc.pp.normalize_total(adata, target_sum=1e4)


# In[29]:


sc.pp.log1p(adata)


# In[30]:


sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


# In[31]:


sc.pl.highly_variable_genes(adata)


# In[32]:


adata.raw = adata


# In[33]:


adata = adata[:, adata.var.highly_variable]


# In[34]:


adata


# In[35]:


sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])


# In[36]:


sc.pp.scale(adata, max_value=10)


# In[37]:


sc.tl.pca(adata, svd_solver='arpack')


# In[38]:


sc.pl.pca(adata, color='rpl39')


# In[39]:


sc.pl.pca_variance_ratio(adata, log=True)


# In[40]:


adata.write(results_file)


# In[41]:


adata


# In[42]:


sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)


# In[43]:


sc.tl.umap(adata)


# In[44]:


sc.pl.umap(adata, color=['rpl39', 'stmn1b', 'hsp70l'])


# In[45]:


sc.pl.umap(adata, color=['rpl39', 'stmn1b', 'hsp70l'], use_raw=True)


# In[46]:


sc.tl.leiden(adata)


# In[61]:


rcParams['figure.figsize'] = 10,10
sc.pl.umap(adata, color=['leiden', 'stmn1b', 'hsp70l'])


# In[63]:


rcParams['figure.figsize'] = 10,10
sc.pl.umap(adata, color='leiden')


# In[48]:


adata.write(results_file)


# In[49]:


sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[50]:


sc.settings.verbosity = 2  # reduce the verbosity


# In[51]:


sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[52]:


adata.write(results_file)


# In[53]:


sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





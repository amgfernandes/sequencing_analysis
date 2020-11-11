#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc


# In[2]:


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


# In[24]:


results_file = 'D:/scRNA-seq_demo_data/test.h5ad'  # the file that will store the analysis results


# In[61]:


adata = sc.read_10x_mtx(
    'D:/scRNA-seq_demo_data/sample_1/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading


# In[62]:


adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`


# In[63]:


adata


# In[65]:


sc.pl.highest_expr_genes(adata, n_top=30, )


# In[67]:


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


# In[68]:


adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[69]:


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


# In[70]:


sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')


# In[71]:


adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]


# In[72]:


sc.pp.normalize_total(adata, target_sum=1e4)


# In[73]:


sc.pp.log1p(adata)


# In[74]:


sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


# In[75]:


sc.pl.highly_variable_genes(adata)


# In[76]:


adata.raw = adata


# In[77]:


adata = adata[:, adata.var.highly_variable]


# In[78]:


sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])


# In[79]:


sc.pp.scale(adata, max_value=10)


# In[80]:


sc.tl.pca(adata, svd_solver='arpack')


# In[81]:


sc.pl.pca(adata, color='rpl39')


# In[82]:


sc.pl.pca_variance_ratio(adata, log=True)


# In[83]:


adata.write(results_file)


# In[84]:


adata


# In[85]:


sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)


# In[86]:


sc.tl.umap(adata)


# In[88]:


sc.pl.umap(adata, color=['rpl39', 'stmn1b', 'hsp70l'])


# In[89]:


sc.pl.umap(adata, color=['rpl39', 'stmn1b', 'hsp70l'], use_raw=True)


# In[90]:


sc.tl.leiden(adata)


# In[92]:


sc.pl.umap(adata, color=['leiden', 'stmn1b', 'hsp70l'])


# In[93]:


adata.write(results_file)


# In[94]:


sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[95]:


sc.settings.verbosity = 2  # reduce the verbosity


# In[96]:


sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[97]:


adata.write(results_file)


# In[98]:


sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os              
os.environ['PYTHONHASHSEED'] = '0'
import desc          
import pandas as pd                                                    
import numpy as np                                                     
import scanpy as sc                                                                              
from time import time                                                       
import sys
import matplotlib
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
sc.settings.set_figure_params(dpi=300)
from matplotlib import rcParams


# In[2]:


print(sys.version)


# In[3]:


sc.__version__


# In[4]:


desc.__version__


# In[5]:


import tensorflow as tf
tf.__version__


# In[6]:


adata = sc.read_10x_mtx(
    '/Users/fernandes/Documents/scRNA-seq_demo_data/sample_6/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                          # write a cache file for faster subsequent reading


# In[7]:


sc.pp.log1p(adata)
#sc.pp.filter_genes_dispersion(adata,n_top_genes=1000) #older scanpy
sc.pp.highly_variable_genes(adata,n_top_genes=1000,subset=True,inplace=True)


# In[8]:


adata.var['mt'] = adata.var_names.str.startswith('mt')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# In[9]:


sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


# In[10]:


desc.normalize_per_cell(adata, counts_per_cell_after=1e4)


# In[11]:


adata.raw=adata


# In[12]:


desc.scale(adata, zero_center=True, max_value=3)


# In[13]:


sc.pp.scale(adata,max_value=6)# if the the dataset has two or more batches you can use `adata=desc.scale(adata,groupby="BatchID")`
save_dir="test_DESC"
adata=desc.train(adata,
        dims=[adata.shape[1],64,32],
        tol=0.005,
        n_neighbors=10,
        batch_size=256,
        louvain_resolution=[1.0],# not necessarily a list, you can only set one value, like, louvain_resolution=1.0
        save_dir=str(save_dir),
        do_tsne=True,
        learning_rate=200, # the parameter of tsne
        use_GPU=False,
        num_Cores=1, #for reproducible, only use 1 cpu
        num_Cores_tsne=4,
        save_encoder_weights=False,
        save_encoder_step=3,# save_encoder_weights is False, this parameter is not used
        use_ae_weights=False,
        do_umap=True) #if do_umap is False, it will don't compute umap coordiate


# In[14]:


adata


# In[15]:


adata.obs['max.prob']=adata.uns["prob_matrix1.0"].max(1)


# In[16]:


sc.pl.scatter(adata,basis="tsne1.0",color=['desc_1.0','max.prob'])


# In[17]:


rcParams['figure.figsize'] = 10,10
sc.pl.umap(adata,color=['desc_1.0','max.prob'])


# In[18]:


sc.pp.neighbors(adata, n_pcs = 30, n_neighbors = 20)

sc.tl.leiden(adata)


# In[19]:


rcParams['figure.figsize'] = 5,5
sc.pl.umap(adata,color=['desc_1.0'])


# In[20]:


rcParams['figure.figsize'] = 5,5
sc.pl.umap(adata,color=['leiden'])


# In[21]:


rcParams['figure.figsize'] = 10,10
sc.pl.umap(adata, color='leiden')


# In[22]:


adata.write('desc_result.h5ad')


# ### One thing virtually all clustering or community detection methods have in common is some flavor of a resolution parameter. This parameter controls how fine- or coarse-grained the inferred clusters are. This parameter can have major effects on your results!
# 
# '''see https://chanzuckerberg.github.io/scRNA-python-workshop/analysis/04-clustering.html for details'''

# In[24]:


sc.tl.leiden(adata, key_added = "leiden_1.0") # default resolution in 1.0
sc.tl.leiden(adata, resolution = 0.6, key_added = "leiden_0.6")
sc.tl.leiden(adata, resolution = 0.4, key_added = "leiden_0.4")
sc.tl.leiden(adata, resolution = 1.4, key_added = "leiden_1.4")


# In[25]:


sc.pl.umap(adata, color=['leiden_0.4', 'leiden_0.6', 'leiden_1.0','leiden_1.4'])


# In[26]:


sc.tl.louvain(adata, key_added = "louvain_1.0") # default resolution in 1.0
sc.tl.louvain(adata, resolution = 0.6, key_added = "louvain_0.6")
sc.tl.louvain(adata, resolution = 0.4, key_added = "louvain_0.4")
sc.tl.louvain(adata, resolution = 1.4, key_added = "louvain_1.4")

sc.pl.umap(adata, color=['louvain_0.4', 'louvain_0.6', 'louvain_1.0','louvain_1.4'])


# In[27]:


adata


# In[28]:


from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

# extract pca coordinates
X_pca = adata.obsm['X_pca'] 

# kmeans with k=5
kmeans = KMeans(n_clusters=5, random_state=0).fit(X_pca) 
adata.obs['kmeans5'] = kmeans.labels_.astype(str)

# kmeans with k=10
kmeans = KMeans(n_clusters=10, random_state=0).fit(X_pca) 
adata.obs['kmeans10'] = kmeans.labels_.astype(str)

# kmeans with k=15
kmeans = KMeans(n_clusters=15, random_state=0).fit(X_pca) 
adata.obs['kmeans15'] = kmeans.labels_.astype(str)

sc.pl.umap(adata, color=['kmeans5', 'kmeans10', 'kmeans15'])


# In[29]:


from sklearn.cluster import AgglomerativeClustering

cluster = AgglomerativeClustering(n_clusters=5, affinity='euclidean', linkage='ward')
adata.obs['hclust_5'] = cluster.fit_predict(X_pca).astype(str)

cluster = AgglomerativeClustering(n_clusters=10, affinity='euclidean', linkage='ward')
adata.obs['hclust_10'] = cluster.fit_predict(X_pca).astype(str)

cluster = AgglomerativeClustering(n_clusters=15, affinity='euclidean', linkage='ward')
adata.obs['hclust_15'] = cluster.fit_predict(X_pca).astype(str)


sc.pl.umap(adata, color=['hclust_5', 'hclust_10', 'hclust_15'])


# In[ ]:





# In[ ]:





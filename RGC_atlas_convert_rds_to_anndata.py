#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# import os
# os.environ["R_HOME"] = 'C:/Program Files/R/R-4.0.3'
# os.environ["PATH"] = 'C:/Program Files/R/R-4.0.3/bin/x64'+ ";" + os.environ["PATH"]


# In[24]:


import scanpy as sc
import numpy as np

# import anndata2ri
# from rpy2.robjects import r
# from rpy2.robjects.conversion import localconverter


# In[47]:


adata1 = sc.read(
    'D:/RGC_atlas/larva_zFish_FINAL.rds',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)   


# In[43]:


import anndata2ri
anndata2ri.activate()
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[45]:


get_ipython().run_line_magic('R', '-i adata_paul')


# In[33]:


import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
readRDS = robjects.r['readRDS']


# In[ ]:


ls


# In[38]:


df = readRDS('D:/RGC_atlas/larva_zFish_FINAL.rds')
df = pandas2ri.ri2py(df)


# In[ ]:





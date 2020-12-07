#!/usr/bin/env python
# coding: utf-8

# In[3]:


import scanpy as sc
import numpy as np

import anndata2ri


# In[4]:


anndata2ri.activate()
get_ipython().magic(u'load_ext rpy2.ipython')


# # Run R libraries to load datasets

# In[16]:


get_ipython().run_cell_magic(u'R', u'', u'suppressPackageStartupMessages(library(Seurat))\nsuppressPackageStartupMessages(library(SingleCellExperiment))')


# In[18]:


get_ipython().run_cell_magic(u'R', u'', u'# load the 10x data:\nlarvadata<- readRDS("~/desktop/ZebrafishRGC_data/larva_zFish_FINAL.rds")')


# In[19]:


get_ipython().run_cell_magic(u'R', u'', u'adultdata<- readRDS("~/desktop/ZebrafishRGC_data/adult_zFish_FINAL.rds")')


# In[7]:


get_ipython().run_cell_magic(u'R', u'-o larva_converted', u'#convert the Seurat object to a SingleCellExperiment object\nlarva_converted <- as.SingleCellExperiment(larvadata)\n\nlarva_converted')


# In[20]:


mito_genes = larva_converted.var_names.str.startswith('mt-')
larva_converted.obs['percent_mito'] = np.sum(larva_converted[:, mito_genes].X, axis=1).A1 / np.sum(larva_converted.X, axis=1).A1

sc.pl.violin(larva_converted, ['nFeature_RNA', 'nCount_RNA', 'percent_mito'], jitter=0.4, multi_panel=True)


# In[23]:


get_ipython().run_cell_magic(u'R', u'', u'\nSaveH5Seurat(adultdataadultdata, filename = "adultdata.h5Seurat")\nConvert("adultdata.h5Seurat", dest = "h5ad")')


# In[8]:


get_ipython().run_cell_magic(u'R', u'-o adult_converted', u'#convert the Seurat object to a SingleCellExperiment object\nadult_converted <- as.SingleCellExperiment(adultdata)\n\nadult_converted')


# In[15]:


larva_converted.write('larva_converted.hdf5')


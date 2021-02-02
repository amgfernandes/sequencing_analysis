#!/usr/bin/env python
# coding: utf-8

# %%:
from anndata._core.aligned_mapping import I
from anndata._core.views import as_view
import numpy as np
from numpy.lib.npyio import save
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scanpy.external as sce
from matplotlib import rcParams
from scrublet.helper_functions import make_genes_unique
from helper_sequencing import Sequencing
from helper_sequencing import load_samples
import itertools as it
import anndata
import seaborn as sns
import scrublet as scr
import os
# %%:
# verbosity: errors (0),warnings (1),info(2),hints(3)
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor='white')
# %%:
'''The file that will store the analysis results'''
# results_file = '/home/fernandes/sample_data/scanpy_test.h5ad'

results_file = '/home/fernandes/RGC_scRNAseq_analysis/\
adult/D_rerio.GRCz11.102.h5ad'

#results_file = '/home/fernandes/RGC_scRNAseq_analysis/\
#larva/D_rerio.GRCz11.102.h5ad'

# Read in the count matrix into an `AnnData
# <https://anndata.readthedocs.io/en/latest/anndata.AnnData.html>`__ object,
#  which holds many slots for annotations and different representations of
# the data. It also comes with its own HDF5 file format: .h5ad.

# In[4]:
# List to remove immediate early genes (IEGs)'''
'''IEGs expression can be elevated in certain cells as a response to the
dissection and dissociation of the cells.This can hamper the analysis of the
cell types and therefore we will exclude IEGs from the analysis (see this
paper for additional information: https://pubmed.ncbi.nlm.nih.gov/29024657/).
We have curated a list of zebrafish known IEGs, and we’ll use this function
to exclude them:'''

IEG_list = pd.read_csv(
    '/home/fernandes/RGC_scRNAseq_analysis/IEG_list.csv', header=None)
IEG_list.columns = ['gene']
IEG_list.gene.values

# %% Load samples:
'''Load samples'''
# '''give a list of sample folders'''

'''for adult data'''
data_list = ["RGC10_S6_L001", "RGC11_S1_L001",
             "RGC11_S1_L002", "RGC12_S2_L001",
             "RGC12_S2_L002", "RGC13_S3_L001",
             "RGC13_S3_L002", "RGC14_S4_L001",
             "RGC14_S4_L002", "RGC15_S5_L001",
             "RGC15_S5_L002", "RGC16_S6_L001",
             "RGC16_S6_L002", "RGC17_S1_L008",
             "RGC18_S2_L008", "RGC19_S3_L008",
             "RGC20_S4_L008", "RGC5_S1_L001",
             "RGC5_S1_L002", "RGC6_S2_L001",
             "RGC6_S2_L002", "RGC7_S3_L001",
             "RGC7_S3_L002", "RGC8_S4_L001",
             "RGC8_S4_L002", "RGC9_S5_L001"]

#for larval data
# data_list = ["ZebraFishRGC1_S1_L001", "ZebraFishRGC2_S1_L005",
#             "ZebraFishRGC3_S2_L001","ZebraFishRGC4_S2_L005"]

''' "ZebraFishRGC4_S2_L005" maybe wetting failure'''

folder_with_data='/home/fernandes/RGC_scRNAseq_analysis/adult/' \
 'D_rerio.GRCz11.102/'

# folder_with_data = '/home/fernandes/RGC_scRNAseq_analysis/larva/' \
#     'D_rerio.GRCz11.102/'

# seq_data=load_samples(data_location=folder_with_data,
# load_method=sc.read_10x_mtx,samplelist=data_list)
seq_data = load_samples(data_location=folder_with_data,
                        load_method=sc.read_mtx, samplelist=data_list)
seq_helper = Sequencing(seqdata=seq_data)

'''remove a set of genes True or False'''
remove=False
if remove:
    print('Removing set of genes given as a list')
    seq_data_filtered = seq_helper.remove_gene_list(seq_data, IEG_list)
else:
    seq_data_filtered=seq_data
    print('No Gene list given to remove')



# %%
# Save figure (set to True to save)
folder_with_data_plot = '/home/fernandes/RGC_scRNAseq_analysis/adult' \
                        '/D_rerio.GRCz11.102/plots'

# folder_with_data_plot = '/home/fernandes/RGC_scRNAseq_analysis/larva' \
#                         '/D_rerio.GRCz11.102/plots'

sc.settings.autosave = True  # save figures True/False
sc.settings.figdir = folder_with_data_plot
sc.settings.set_figure_params(dpi_save=300, format="png")

outputDirectory = folder_with_data_plot
if not os.path.isdir(outputDirectory):
    os.mkdir(outputDirectory)

# %%

'''iterate over all elements on list of samples and concatenate them...only
keep last result '''
# Not working as it should for now TODO: Fix
# adata=list(it.accumulate(seq_data_filtered, sc.AnnData.concatenate))[-1]

'''concatenate batches'''
# for larval data
# adata = seq_data_filtered[0].concatenate(
#     seq_data_filtered[1], seq_data_filtered[2], seq_data_filtered[3])

# for adult data
adata = seq_data_filtered[0].concatenate(
    seq_data_filtered[1], seq_data_filtered[2], seq_data_filtered[3],
    seq_data_filtered[4], seq_data_filtered[5], seq_data_filtered[6],
    seq_data_filtered[7], seq_data_filtered[8], seq_data_filtered[9],
    seq_data_filtered[10], seq_data_filtered[11], seq_data_filtered[12],
    seq_data_filtered[13], seq_data_filtered[14], seq_data_filtered[15],
    seq_data_filtered[16], seq_data_filtered[17], seq_data_filtered[18],
    seq_data_filtered[19], seq_data_filtered[20], seq_data_filtered[21],
    seq_data_filtered[22], seq_data_filtered[23], seq_data_filtered[24],
    seq_data_filtered[25])


'''create a list of genes'''
list_of_genes = list(adata.var_names)
# %%
'''remove objects not used further'''

#del (seq_data_filtered, seq_data)

# %%
'''Test for library saturation'''
# Create a plot showing genes detected as a function of UMI counts.
# fig, ax = plt.subplots(figsize=(10, 7))

# x = np.asarray(adata.X.sum(axis=1))[:, 0]
# y = np.asarray(np.sum(adata.X > 0, axis=1))[:, 0]

# ax.scatter(x, y, color="green", alpha=0.25)
# ax.set_xlabel("UMI Counts")
# ax.set_ylabel("Genes Detected")
# ax.set_xscale('log')
# ax.set_yscale('log', nonposy='clip')

# ax.set_xlim((0.5, 4500))
# ax.set_ylim((0.5, 2000))
# plt.savefig(folder_with_data_plot+'/library saturation.png')

# %%
''''The vast majority of “cells” have only a few UMI detected.
Those are empty droplets. 10x claims to have cell capture rate of up to 65%,
but in practice, depending on how many cells are in fact loaded, the rate can
be much lower. A commonly used method to estimate the number of empty droplets
is barcode ranking knee and inflection points, as those are often assumed to
represent transition between two components of a distribution.
The “knee plot” is a standard single-cell RNA-seq quality control that is also
used to determine a threshold for considering cells valid for analysis in an
experiment.
More sophisticated method exist (e.g. see emptyDrops in DropletUtils) '''

knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
fig, ax = plt.subplots(figsize=(10, 7))

# for adult: 200000
# for larva : 20000
expected_num_cells=100000

ax.loglog(knee, range(len(knee)), label="kallisto", linewidth=5, color="k")
ax.axvline(x=knee[expected_num_cells], linewidth=3, color="g")
ax.axhline(y=expected_num_cells, linewidth=3, color="g")

ax.set_xlabel("UMI Counts")
ax.set_ylabel("Set of Barcodes")

plt.grid(True, which="both")
ax.legend()
plt.savefig(folder_with_data_plot+'/knee_plot.png')
print (knee[expected_num_cells])



# %%

'''Basic filtering'''

'''# removing genes that are expressed in
# fewer than 25 cells and removing cells that have fewer than
# 450 features (Harvard defaults)'''

print (adata)
#Based on kneeplot. Remove empty droplets
sc.pp.filter_cells(adata, min_counts=knee[expected_num_cells])
#Minimum number of genes expressed required for a cell to pass filtering.
sc.pp.filter_cells(adata, min_genes=450)
#Minimum number of cells expressed required for a gene to pass filtering.
sc.pp.filter_genes(adata, min_cells=25)

# %%
'''Statistics of loaded data'''
print (adata.var.describe())
print (adata.obs.describe())
n_cells=adata.var['n_cells'].count()
n_genes=adata.obs['n_genes'].count()

d={'n_cells':n_cells,'n_genes':n_genes}
df_stats=pd.DataFrame(d, index=[0])
df_stats.to_csv(folder_with_data_plot+'/stats.csv')

# %%
# Preprocessing
# Show those genes that yield the highest fraction of counts in each single
# cells, across all cells.
sc.pl.highest_expr_genes(adata, n_top=20)

# %%
# annotate the group of mitochondrial genes as 'mt'
adata.var['mt'] = adata.var_names.str.startswith('mt')
sc.pp.calculate_qc_metrics(
    adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# %%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.3, multi_panel=True, save='_counts.png')


# %%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',
              save='pct_counts_mt.png')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
              save='n_genes_by_counts.png')
# %%
'''Predict doublets'''
"""
scrub = scr.Scrublet(adata.X)
adata.obs['doublet_scores'], adata.obs['predicted_doublets']\
     = scrub.scrub_doublets()
scrub.plot_histogram()

sum(adata.obs['predicted_doublets'])

# add in column with singlet/doublet instead of True/False
adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)

sc.pl.violin(adata, 'n_genes_by_counts', jitter=0.4,
             groupby='doublet_info', rotation=45, save='_doublets.png') """

# %%
plt.figure()
sns.set_style("whitegrid")
sns.distplot(adata.obs['total_counts'], kde=False)
sns.despine()
plt.savefig(folder_with_data_plot+'/total_counts_distplot.png')


# %%
plt.figure()
sns.set_style("whitegrid")
sns.distplot(adata.var['n_cells'], kde=False)
sns.despine()
plt.savefig(folder_with_data_plot+'/ncells_distplot.png')



# %%

'''filtering based on QC'''
'''unusual number of genes (genes_by_counts), percent.mt and UMI (nCount_RNA).
What is unusual? No ground truth role here. It can vary with the tissue
examined, the efficiency the 10x or the sequencing machine. We can examine
the violin plot and get a feeling of the outliers in our data. Additionally
for mitochondrial genes, the avg. in neurons is about 12.5%, so I maybe the
threshold should be higher to capture most of the population variance. The most
important part is to check how the downstream analysis looks like, and return
and refine the threshold set here. For example, if in the downstream analysis
mitochondrial genes are highly present in the variable genes, the bar was set
too high. If some cells have marker genes of two different cluster, maybe there
are a lot duplets in the data and the threshold for gene number and UMI should
be lower.'''



#adata = adata[adata.obs.n_genes < 2500, :]
#adata = adata[adata.obs.n_genes > 450, :]

# remove cells with high mitochondrial content
adata = adata[adata.obs.pct_counts_mt < 20, :]
adata
# %%

sc.pp.normalize_total(adata, target_sum=1e4)

# %%
'''Logarithmize the data.'''
sc.pp.log1p(adata)


# %%

# Scale each gene to unit variance. Clip values exceeding standard deviation 10
sc.pp.scale(adata, max_value=10)
# %%
'''Identify highly-variable genes.'''
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
print("Highly variable genes: %d"%sum(adata.var.highly_variable))

# %%:
sc.pl.highly_variable_genes(adata)


# %%
'''Set the .raw attribute of AnnData object to the normalized and
logarithmized raw gene expression for later use in differential testing
and visualizations of gene expression.This simply freezes the state of the
AnnData object.'''

adata.raw = adata
# %%
'''filter only highly variable genes'''
adata = adata[:, adata.var.highly_variable]
adata.write(results_file)

# %%

# Regress out effects of total counts per cell and the percentage of
# mitochondrial genes expressed. Scale the data to unit variance.

#sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# %%
'''make a copy to use for batch correction'''
adata_corr = adata.copy()
# %%

'''#Harmony works by adjusting the principal components, this function should
 be run after performing PCA but before computing the neighbor graph,'''
sc.tl.pca(adata_corr, svd_solver='arpack')


# %%

sc.pl.pca_variance_ratio(adata_corr, log=True)

# %%
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
# %%
sce.pp.harmony_integrate(adata_corr, 'batch')  # correct batch effect
'X_pca_harmony' in adata_corr.obsm

# %%

# %%
# adata.var_names_make_unique()
# this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
# adata_corr.var_names_make_unique()

# %%
sc.pp.neighbors(adata_corr, n_neighbors=10, n_pcs=40)


# %%]:
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# %%

sc.tl.umap(adata)
sc.pl.umap(adata, color=['batch'],save='_no_batch_correction.png', s=10)
# %%:

sc.tl.umap(adata_corr)
sc.pl.umap(adata_corr, color=['batch'], save='_batch_corrected.png', s=10)

# %%
sc.tl.leiden(adata_corr, key_added='clusters', resolution=1.3)
sc.tl.umap(adata_corr)

# %%:
'''Embedding the neighborhood graph
We suggest embedding the graph in two dimensions using UMAP
(McInnes et al., 2018), see below. It is potentially more faithful to the
global connectivity of the manifold than tSNE, i.e., it better preserves
trajectories. In some ocassions, you might still observe disconnected
clusters and similar connectivity violations. They can usually be
remedied by running:'''

""" sc.tl.paga(adata_corr)
sc.pl.paga(adata_corr, plot=False)
# remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata_corr, init_pos='paga')
sc.tl.umap(adata_corr) """

# %%
'''if wanting tsne'''
# sc.tl.tsne(adata_corr)
# sc.pl.tsne(adata_corr,palette='Dark2',color=['clusters'])
# %%
#sc.tl.leiden(adata_corr, key_added='clusters', resolution=1.5)
# %%
rcParams['figure.figsize'] = 10, 10
sc.pl.umap(adata_corr, color=['clusters'], save='Clusters.png', s=15,
            palette='Dark2')
# %%
rcParams['figure.figsize'] = 10, 10
sc.pl.umap(adata_corr, color='clusters', add_outline=True, legend_loc='on data',
           legend_fontsize=12, legend_fontoutline=2, frameon=False,
           title='clustering of cells', palette='Dark2',save='Clusters2.png')

# %%:

sc.pl.violin(adata_corr, ['n_genes_by_counts'],groupby='clusters',
             jitter=0.3, multi_panel=True, save='_n_genes_by_counts.png')


sc.pl.violin(adata_corr, ['n_genes'],groupby='clusters',
             jitter=0.3, multi_panel=True, save='_n_genes.png')          
# %%

adata_corr.write(results_file)


# %%:

sc.pl.umap(adata_corr, color=['robo2', 'isl2b', 'rbpms2b'])

# %%
rcParams['figure.figsize'] = 20,20
sc.pl.umap(adata_corr, color=['isl2b', 'rbpms2b'], save='RGCs_markers.png',s=20)

# %%
rcParams['figure.figsize'] = 20,20
sc.pl.umap(adata_corr, color=['clusters'], legend_loc='on data', save='clusters_to_select_RGCs.png',s=20)
# %%
                                       
sc.pl.dotplot(adata_corr, var_names=['robo2', 'isl2b', 'rbpms2b'], groupby='clusters', figsize=(10, 10), save='clusters_to_select_RGCs.png')

sc.pl.violin(adata_corr, keys=['robo2', 'isl2b', 'rbpms2b'], groupby='clusters', figsize=(10, 10), save='clusters_to_select_RGCs.png')

# %%
'''Remove clusters that are not neurons. Based on low expression of 'isl2b', 'rbpms2b'''
rgcs= adata_corr[~adata_corr.obs['clusters'].isin(['5','7','10','13','18','23','24']),:]
sc.pl.umap(rgcs, color=['clusters'], legend_loc='on data', save='RGCs_clusters_only.png',s=20)
# %%
'''We will need to analyze the variable features for the new dataset and rerun the clustering analysis.'''

'''Identify highly-variable genes.'''
sc.pp.highly_variable_genes(rgcs, n_top_genes=2000)
print("Highly variable genes: %d"%sum(rgcs.var.highly_variable))


sc.tl.pca(rgcs, svd_solver='arpack')
sc.pl.pca_variance_ratio(rgcs, log=True)

sce.pp.harmony_integrate(rgcs, 'batch')  # correct batch effect
'X_pca_harmony' in rgcs.obsm

# %%
sc.pp.neighbors(rgcs, n_neighbors=10, n_pcs=40)
sc.tl.umap(rgcs)
sc.pl.umap(rgcs,color=['batch'], s=20, save='rgcs_batch.png')

# %%
sc.tl.leiden(rgcs, key_added='clusters', resolution=1.5)
sc.pl.umap(rgcs, color='clusters', s=20, save='RGCs_clusters_after_reclustering.png')

# %%
'''explore data'''
sc.tl.dendrogram(rgcs, 'clusters')
sc.tl.rank_genes_groups(
    rgcs,
    groupby='clusters',
    n_genes=rgcs.shape[1],
    method='wilcoxon')

sc.pl.rank_genes_groups_dotplot(rgcs, n_genes=2, save='rank_genes_rgcs.png')
# %%
sc.pl.umap(rgcs, color=['mafaa', 'eomesa', 'tbr1b'],
           save='selected_genes.png')
sc.pl.umap(rgcs, color=['clusters'])
# %%

sc.pl.dotplot(rgcs, var_names=['eomesa', 'tbr1b', 'mafaa', 'neurod1', 'epha7',
                                'id2b', 'tbx20', 'onecut1'],
              groupby='clusters', save='selected_genes.png', figsize=(8, 10))
# %%
sc.pl.heatmap(rgcs, var_names=['eomesa', 'tbr1b', 'mafaa', 'neurod1', 'epha7',
                                'id2b', 'tbx20'],
              groupby='clusters', save='selected_genes.png', figsize=(8, 10))


# %%
sc.pl.umap(
        rgcs,
        color=['eomesa', 'tbr1b', 'mafaa', 'neurod1', 'epha7', 'id2b', 'tbx20', 'onecut1'],
        color_map=None,
        save='_selected_genes.png', s=20)
# %%
cluster_to_check = '26'
clust_look = rgcs[rgcs.obs['clusters'].values.isin([cluster_to_check])]
clust_look.obs['cluster_to_check'] = cluster_to_check

# %%
sc.tl.rank_genes_groups(
    clust_look,
    groupby='cluster_to_check',
    n_genes=clust_look.shape[1],
    method='wilcoxon')
sc.pl.rank_genes_groups_dotplot(clust_look, n_genes=20)
# %%
sc.tl.leiden(clust_look,key_added='clusters',resolution=0.5)
sc.pl.umap(clust_look, color='clusters')
sc.tl.rank_genes_groups(clust_look, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(clust_look, n_genes=1)
# %%
sc.pl.dotplot(clust_look, var_names=['mafaa', 'epha7',
                                     'id2b'],
              groupby='leiden', save='selected_genes_subcluster.png',
              figsize=(8, 8))
# %%
marker_genes_dict = {'Prey': ['epha7', 'id2b', 'mafaa'],
                    'Phototaxis': ['eomesa', 'tbx20'],
                    'Other markers': ['tbr1b', 'onecut1', 'shisa9b']}
sc.pl.dotplot(rgcs, marker_genes_dict, groupby='clusters', dendrogram=True,
              figsize=(10, 10), save='_marker_genes.png')

# %%

sc.pl.umap(rgcs, color=['clusters'], s=5)
# %%

# %%

# %%

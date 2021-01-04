__authors__ = 'fernandes Dec 2020'
#TODO finish helper class. Add metadata possibility

def load_samples(load_method=None,samplelist=None,cache=True):
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
        data=load_method(
        '/home/fernandes/sample_data/' + sample, # the directory with the `.mtx` file
        var_names='gene_symbols',
        cache=cache)
        print (n)
        data.obs['batch']=n
        samples.append(data)
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

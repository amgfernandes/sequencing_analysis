__authors__ = 'fernandes Dec 2020'
#TODO finish helper class


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
    #Function receives a  gene matrix and a list of genes and returns a  gene matrix without the gene rows.
    #The loaded 10x data is a dgc matrix.
        idx=self.return_indices_not_in_list(dgc_mat.var_names,gene_list.gene.values)
        dgc_mat_corr= dgc_mat[:,idx]

        return (dgc_mat_corr)
    
    def keep_gene_list(self,dgc_mat, gene_list):
        #Function receives a  gene matrix and a list of genes and returns a  gene matrix without the gene rows.
    #The loaded 10x data is a dgc matrix.
        idx=self.return_indices_inside_list(dgc_mat.var_names,gene_list.gene.values)
        dgc_mat_corr= dgc_mat[:,idx]

        return (dgc_mat_corr)

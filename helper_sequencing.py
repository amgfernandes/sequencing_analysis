#TODO finish helper class
class Sequencing():
        """ Class with functions to help sequencing analysis
        """
    def __init__(self, seqdata, stimuli, metadata):
        """ Loads the experiment from sequencing data
        :param seqdata: sequencing data
        :param metadata: metadata dictionary, for the sequencing data
        """
        self.seqdata = seqdata
        self.metadata = metadata

    def return_indices_of_a(a, b):
        b_set = set(b)
        return [i for i, v in enumerate(a) if v not in b_set]

    def remove_genes(dgc_mat, gene_list):
    #Function receives a  gene matrix and a list of genes and returns a  gene matrix without the gene rows.
    #The loaded 10x data is a dgc matrix.
        idx=return_indices_of_a(dgc_mat.var_names,gene_list.gene.values)
        dgc_mat_corr= dgc_mat[:,idx]

        return (dgc_mat_corr)

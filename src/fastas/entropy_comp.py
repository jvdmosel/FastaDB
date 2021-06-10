from fastas.fasta_entropy import FastaEntropy

class EntropyComp:
    """Class which compares to FastaEntropy objects to another
    
    """

    def __init__(self, filename_1, filename_2):
        """Constructor of EntropyComp class
        
        Arguments:
            filename_1 {string} -- name of first file to compare
            filename_2 {string} -- name of second file to compare

        Attributes:
            a {list} -- entropy of first file
            b {list} -- entropy of second file
            name_a {string} -- id of first file
            name_b {string} -- id of second file
        """

        # create new FastaEntropy objects
        h1 = FastaEntropy(filename_1)
        h2 = FastaEntropy(filename_2)
        # calculate each entropy
        self.a = h1.get_entropy()
        self.b = h2.get_entropy()
        self.name_a = filename_1.split('_')[0]
        self.name_b = filename_2.split('_')[0]
    
    def sort_by_max(self, h):
        """Gives back indices of h when sorted by entropy in decreasing order
        
        Arguments:
            h {list} -- entropy of a fasta file
        
        Returns:
            indices [list] -- list of indices sorted by entropy in decreasing order
        """

        return sorted(range(len(h)), key=lambda k: h[k], reverse=True)
            
    def __repr__(self):
        """string representation of EntropyComp object
        
        Returns:
            repr [string] -- as table
        """

        a_max = self.sort_by_max(self.a)
        # make it pretty
        table = "{0:5}   {1:9}   {2:9}"
        comp = '\n'+table.format("", self.name_a, self.name_b)+'\n'
        table = "{0:5} | {1:9} | {2:9}"
        comp += table.format("INDEX", "H1", "H2")+'\n'
        # the five indices with the highest entropy, +1 to make the output 1-based
        for i in range(0,5) :
            comp += table.format((a_max[i]+1), '{:.4f}'.format(self.a[a_max[i]]), '{:.4f}'.format(self.b[a_max[i]]))+'\n'
        return comp
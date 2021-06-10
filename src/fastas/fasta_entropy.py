from .fasta_parser import FastaParser
from pathlib import Path
import math

AMINOACIDS = ['A', 'R', 'N', 'D', 
              'C', 'Q', 'E', 'G', 
              'H', 'I', 'L', 'K', 
              'M', 'F', 'P', 'S', 
              'T', 'W', 'Y', 'V',
              'X', '-']

class FastaEntropy:
    """Class which represents the Shannon entropy of an aligned fasta file
    
    """

    def __init__(self, filename):
        """Constructor of FastaEntropy class
        
        Arguments:
            filename {string} -- name of the fasta file to calculate the entropy of

        Attributes:
            freqmatrix {list<dict>} -- frequency matrix
            probmatrix {list<dict>} -- probability matrix
            entropy {list}          -- entropy of each column of the file
            numofseq {int}          -- number of sequences in the file
            seqlength {int}         -- length of sequences in the file
            filepath {Path}         -- path to file

        """

        self.freqmatrix = []
        self.probmatrix = []
        self.entropy = []
        self.numofseq = 0
        self.seqlength = 0
        self.filepath = Path('out') / filename

    def init_matrices(self, n):
        """initializes the matrices, entropy and sequence length
        
        Arguments:
            n {int} -- length of sequences in the file
        """

        self.freqmatrix = [dict.fromkeys(AMINOACIDS, 0) for i in range(n)]
        self.probmatrix = [dict.fromkeys(AMINOACIDS, 0.0) for i in range(n)]
        self.entropy = [0.0 for i in range(n)]
        self.seqlength = n

    def get_freqMatrix(self):
        """Calculates the frequency of each aminoacid at each column in the file
 
        """

        parser = FastaParser(str(self.filepath))
        init = True
        for f in parser:
            self.numofseq += 1
            seq = f.get_sequence()
            if init:
                self.init_matrices(len(seq))
                init = False
            i = 0
            for amino in seq:
                self.freqmatrix[i][amino] += 1
                i += 1

    def get_probMatrix(self):
        """Calculates the probability of each aminoacid at each column in the file

        """

        i = 0
        for freq in self.freqmatrix:
            for amino in AMINOACIDS:
                self.probmatrix[i][amino] = freq[amino]/self.numofseq
            i += 1
            
    def normalize(self):
        """Normalizes all entropies 

        """

        norm = math.log2(self.numofseq)
        for i in range(len(self.entropy)):
            self.entropy[i] = self.entropy[i]/norm

    def calc_h(self, prob):
        """Helper function to calculate the shannon entropy 

        Returns:
            entropy {float} -- shannon entropy of a column
        """
        return -(prob*(math.log2(prob) if prob>0 else 0))

    def get_entropy(self):
        """ Calculates the entropy for each column 

        """

        # get matrices
        self.get_freqMatrix()
        self.get_probMatrix()
        i = 0
        # calculate entropy, ignore gaps
        for prob in self.probmatrix:
            for amino in AMINOACIDS:
                if amino != '-':
                    self.entropy[i] += self.calc_h(prob[amino])
            i += 1
        # normalize
        self.normalize()
        return self.entropy

    def __repr__(self):
        """String representation of FastaEntropy
        
        Returns:
            repr [string] -- representation
        """

        self.get_entropy()
        h = ""
        for i in self.entropy:
            # round to two digits and seperate by pipes
            h += (str('%.2f' % i) + ' | ')
        return h
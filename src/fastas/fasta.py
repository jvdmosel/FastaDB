class Fasta:
    """Fasta class represents a .fasta datum 

    """
    
    def __init__(self, gen, spec, subs, fact, clas):
        """Constructor of Fasta class, initializes new
        Fasta datum with description and sequence lines
        
        Arguments:
            gen {str} -- genus e.g. Canis
            spec {str} -- species e.g. lupus
            subs {str} -- sub-species e.g. familiaris might be NONE
            fact {str} -- transcription factor
            clas {str} -- class e.g. mammalia(ma)

        Attributes:
            genus {str}          -- gen
            species {str}        -- spec
            subspecies {str}     -- subs
            factor {str}         -- fact
            classification {str} -- clas
            sequences {str list} -- sequence of the fasta datum, might be multiple lines long
        """

        self.genus = gen
        self.species = spec        
        # check whether instance of Fasta object has a sub-species defined 
        self.subspecies = '' if subs == None else subs
        self.factor = fact 
        # cut the newline symbol
        clas = clas[:len(clas)-len("\n")]
        self.classification = clas[:2]
        self.sequences = []

    def add_sequence(self, seq):
        """Appends line to the sequence list
        
        Arguments:
            seq {str} -- sequence to append
        """

        self.sequences.append(seq)

    def get_sequence(self):
        """Concatenates own sequence list to string and returns it
        
        Returns:
            seq {str} -- own sequence as string
        """

        seq = ''.join(self.sequences)
        seq = seq[:len(seq)-len("\n")]
        return seq

    def get_fullspecies(self):
        """Concatenates genus,species,sub-species seperated with 
        underscores and returns it
         
        Returns:
            str -- genus_species(_subspecies)
        """

        # make sure there won't be a double underscore if there is no subspecies
        subspecies = '_'+self.subspecies if self.subspecies != '' else self.subspecies
        return self.genus+'_'+self.species+subspecies

    def get_factor(self):
        """Returns own transcription factor
        
        Returns:
            str -- transcription factor
        """

        return self.factor

    def get_class(self):
        """Returns own classification
        
        Returns:
            str -- class
        """

        return self.classification
    
    def __repr__(self):
        """Unique string representation of Fasta datum object
        
        Returns:
            str -- description line and sequence line
        """

        description = self.get_fullspecies()+'_'+self.factor+'_'+self.classification
        # make list of all attributes and sequences
        fasta_repr = [description, ]
        fasta_repr.extend(self.sequences)
        # list to string
        return ''.join(fasta_repr)
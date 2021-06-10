from .fasta import Fasta

class FastaParser:
    """Iterable parser class for .fasta files which yields
    Fasta objects
    """
    
    def __init__(self, path):
        """FastaParser constructor 
        
        Arguments:
            path {Path} -- path to file
        """

        # string representation of Path object
        self.path = str(path)

    def __iter__(self):
        """Makes FastaParser iterable, reads file specified by path attribute
        and yields Fasta objects

        Yields:
            fasta_obj {Fasta} -- object of Fasta class
        """

        fasta_obj = None
        # read file from path
        with open(self.path, 'r') as fasta_reader:
            # iterate over every line
            for line in fasta_reader:
                # declaration line is marked with '>' symbol
                if line[0] == '>':
                    # make sure yield is not called when fasta_obj is not initialized yet
                    if fasta_obj:
                        yield fasta_obj
                    # get rid of the '>' symbol
                    line = line[1:]
                    # split line on undermarks into several strings 
                    attr = line.split('_')
                    # if there are more than five strings then there is a sub-species specified
                    # make sure all arguments are given correctly on constructor call
                    has_subs = (len(attr) >= 5 and attr[2].islower())
                    subs = attr[2] if has_subs else None 
                    prot = attr[3] if has_subs else attr[2]
                    if (len(attr) == 5):
                        clas = attr[4]
                    elif (len(attr) == 6):
                        clas = attr[5]
                    else:
                        clas = attr[3]
                    # initialize new Fasta object
                    fasta_obj = Fasta(attr[0], attr[1], subs, prot, clas)
                # not marked with '>' as decription line, must be a sequence
                else:
                    fasta_obj.add_sequence(line)
            yield fasta_obj


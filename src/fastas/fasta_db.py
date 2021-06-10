import sqlite3
import os
from .fasta_parser import FastaParser
from pathlib import Path

SUPER_ID = 0
CLASS_ID = 1
FAMILY_ID = 2
SUBFAM_ID = 3
GENUS_ID = 4
SPECIES = 5
FACTOR = 6
CLASSIFICATION = 7
SEQUENCE = 8
CLASS_SEQ = 9
FAMILY_SEQ = 10
SUBFAM_SEQ = 11

class FastaDB:
    """Database class which creates a sqlite database from fasta files 
    and provides access to it 

    """

    def __init__(self, path):
        """Constructor of FastaDB class, creates new sqlite database or connects to
        existing one
        
        Arguments:
            path {Path} -- path to the location where the database is or will be stored
                           (if not created yet)

        Attributes:
            connection {connection object} -- represents the database
            cursor {cursor object}         -- cursor to call execute methods on to perform SQL commands 
            map {dict}                     -- representation of the name2ID.txt file, maps tf name to ID
        """

        self.connection = sqlite3.connect(str(path))
        self.cursor = self.connection.cursor()
        self.map = {}
        # initialize map from name2ID file
        with open('src/fastas/name2ID.txt', 'r') as map_reader:
            for line in map_reader:
                nameToID = line.split(';')
                # make sure all names have the same case
                nameToID[0] = nameToID[0].upper()
                self.map[nameToID[0]] = nameToID[1][:len(nameToID[1])-len("\n")]
    
    def build_table(self):
        """Creates the database table         
        
        Table consists of:
            ID {integer}          -- 5 values/columns
            species {text}        -- species name
            factor {text}         -- tf name
            classification {text} -- class e.g mammalia
            sequence {text}       -- unaligned sequence of species tf
            class_seq {text}      -- level 2 aligned
            family_seq {text}     -- level 3 aligned
            subfam_seq {text}     -- level 4 aligned
            PRIMARY KEY           -- full ID + species name (unique entry)
        """

        sql_name = """CREATE TABLE IF NOT EXISTS fastas(
            super_id integer,
            class_id integer,
            family_id integer,
            subfam_id integer,
            genus_id integer,
            species text,
            factor text,
            classification text,
            sequence text,
            class_seq text,
            family_seq text,
            subfam_seq text,
            PRIMARY KEY(super_id, class_id, family_id, subfam_id, genus_id, species))"""
        self.cursor.execute(sql_name)
        # fill database with values
        self.populate()
        # save database via commit
        self.connection.commit()
    
    def factorToID(self, factor):
        """gets corresponding ID to tf name
        
        Arguments:
            factor {string} -- name of transcription factor
        
        Returns:
            ids [list] -- full ID as array
        """
       
        # fasta files are not consistent in their naming of transcription factors
        # make sure case matches 
        factor = factor.upper() 
        # try to match a factor in map
        while(factor not in self.map and len(factor) != 1):
            # reduce name each time by one place and try again until match is found
            factor = factor[:len(factor)-1]
        ids = self.map[factor].split('.')
        return ids

    def insert_query(self, fasta):
        """takes fasta datum object and inserts its values into the database
        
        Arguments:
            fasta {Fasta} -- fasta datum
        """

        arguments = "(super_id, class_id, family_id, subfam_id, genus_id, species, factor, sequence, classification)"
        # only insert if primary key entry does not exist yet
        query = "INSERT OR IGNORE into fastas " + arguments + " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)"
        # use tf name to get full ID including level 5
        ids = self.factorToID(fasta.get_factor())
        # execute insert query
        self.cursor.execute(query, (ids[0], ids[1], ids[2], ids[3], ids[4], fasta.get_fullspecies(), fasta.get_factor(), fasta.get_sequence(), fasta.get_class()))

    def update_query(self, fasta, size):
        """update entries with sequences from aligned fasta files
        
        Arguments:
            fasta {Fasta} -- aligned fasta datum
            size {int}    -- length of the id, used to get corresponding column (level) 
        """

        # make sure the corresponding entry does exist, if not create a new entry
        self.insert_query(fasta)
        begin = "UPDATE fastas SET"
        mid = "WHERE super_id = ? AND class_id = ? AND family_id = ? AND subfam_id = ? AND genus_id = ?"
        end = "AND species = ?"
        ids = self.factorToID(fasta.get_factor())
        # get correct column
        if size == 2:
            query = begin + " class_seq = ? " + mid + " " + end
        elif size == 3:
            query = begin + " family_seq = ? " + mid + " " + end
        elif size == 4: 
            query = begin + " subfam_seq = ? " + mid + " " + end
        # execute update query
        self.cursor.execute(query, (fasta.get_sequence(), ids[0], ids[1], ids[2], ids[3], ids[4], fasta.get_fullspecies()))
        

    def fillTable(self, path, aligned):
        """iterates over all fasta files in given directory and inserts them
        or updates their entries one by one 
        
        Arguments:
            path {Path} -- path where fasta files are stored
            aligned {bool} -- true if aligned, false if not
        """

        # iterate over all files in path directory
        for filename in path.iterdir():
            fn = str(filename)
            # creates new FastaParser object for file
            parser = FastaParser(fn)
            fn = filename.name
            fn = fn[:len(fn)-len(".fasta")]
            ids = fn.split('.')
            size = len(ids)
            # for every fasta datum in file
            for f in parser:
                # update aligned files, insert non-aligned
                if aligned:
                    self.update_query(f, size)
                else:
                    self.insert_query(f)

    def populate(self):
        """small helper method that calls fillTable two times: non-aligned and
        aligned

        """

        self.fillTable(Path('src/fastas/files'), False)
        self.fillTable(Path('src/fastas/files_aligned'), True)

    def writeToFile(self, fpath, column):
        """creates a new fasta file and fills it with data where the cursor points to
        
        Arguments:
            fpath {Path} -- path where output file will be created
            column {int} -- column of sequence which is asked for (8-11 are possible values)
        
        Returns:
            fetch [bool] -- true if file was successfully created, false if not
        """

        # delete file if it already exists
        if fpath.is_file():
            fpath.unlink()
        filename = str(fpath)
        fetch = False
        # write to file row for row
        for row in self.cursor:
            if not fetch:
                print("Generating file " + filename[4:] + " at ./" + filename[:3])
            fetch = True
            with open(filename, 'a') as f:
                # description line
                f.write('>'+row[SPECIES]+"_"+row[FACTOR]+"_"+row[CLASSIFICATION]+'\n')
                # sequence line
                f.write(row[column]+'\n')
        return fetch

    def alignedQuery(self, ids, id_len):
        """Takes and ID and creates an aligned query from it, which can then
        be used to get the corresponding data from the database
        
        Arguments:
            ids {list}   -- ID of the node which is asked for
            id_len {int} -- length of the id, used to get the correct column
        
        Returns:
            query [list] -- aligned query 
        """

        path = ids[0]+'.'+ids[1]
        query = " IS NOT NULL AND class_id=?"
        # two digits
        if id_len == CLASS_ID:
            query = " AND class_seq" + query
            return (query, (ids[0],ids[1]), CLASS_SEQ, path)
        path += '.'+ids[2]
        query += " AND family_id=?"
        # three digits
        if id_len == FAMILY_ID:
            query = " AND family_seq" + query
            return (query, (ids[0],ids[1],ids[2]), FAMILY_SEQ, path)
        path += '.'+ids[3]
        query += " AND subfam_id=?"
        # four digits
        if id_len == SUBFAM_ID:
            # fourth digit is not a zero -> level 4 aligned
            if ids[3] != '0':
                query = " AND subfam_seq" + query
                return (query, (ids[0],ids[1],ids[2],ids[3]), SUBFAM_SEQ, path)
            else:
                # fourth digit is a zero -> level 3 aligned
                query = " AND family_seq" + query
                return (query, (ids[0],ids[1],ids[2],ids[3]), FAMILY_SEQ, path)
        path += '.'+ids[4]
        query += " AND genus_id=?"
        # five digits
        if id_len == GENUS_ID:
            # fourth digit is not a zero -> level 4 aligned
            if ids[3] != '0':
                query = " AND subfam_seq" + query
                return (query, (ids[0],ids[1],ids[2],ids[3],ids[4]), SUBFAM_SEQ, path)
            else:
                # fourth digit is a zero -> level 3 aligned
                query = " AND family_seq" + query
                return (query, (ids[0],ids[1],ids[2],ids[3],ids[4]), FAMILY_SEQ, path)

    def unalignedQuery(self, ids, id_len):
        """Takes and ID and creates a non-aligned query from it, which can then
        be used to get the corresponding data from the database
        
        Arguments:
            ids {list}   -- ID of the node which is asked for
            id_len {int} -- length of the id, used to get the correct column
        
        Returns:
            query [list] -- non-aligned query 
        """

        path = ids[0]+""
        if id_len == SUPER_ID:
            return ("", (ids[0],), path)
        path += '.'+ids[1]
        query = " AND class_id=?"
        if id_len == CLASS_ID:
            return (query, (ids[0],ids[1]), path)
        path += '.'+ids[2]
        query += " AND family_id=?"
        if id_len == FAMILY_ID:
            return (query, (ids[0],ids[1],ids[2]), path)
        path += '.'+ids[3]
        query += " AND subfam_id=?"
        if id_len == SUBFAM_ID:
            return (query, (ids[0],ids[1],ids[2],ids[3]), path)
        path += '.'+ids[4]
        query += " AND genus_id=?"
        if id_len == GENUS_ID:
            return (query, (ids[0],ids[1],ids[2],ids[3],ids[4]), path)

    def get_node(self, node, aligned):
        """retrieves data for given node and creates output for it
        
        Arguments:
            node {string} -- node which is asked for
            aligned {bool} -- true if aligned, false if non-aligned
        
        Returns:
            fetch [bool] -- true if successful, false if not
        """

        output_path = Path('./out')
        ids = node.split('.')
        query = "SELECT * FROM fastas WHERE super_id=?"
        if aligned and len(ids) > 1:
            args = self.alignedQuery(ids, len(ids)-1)
            column = args[2]
            output_path = output_path / (args[3]+"_mammalia_aligned_fasta.fasta")
        else:
            args = self.unalignedQuery(ids, len(ids)-1)
            column = SEQUENCE
            output_path = output_path / (args[2]+"_mammalia_fasta.fasta")
        self.cursor.execute(query+args[0], args[1])
        return self.writeToFile(output_path, column)            

    def get_species(self, species):
        """retrieves data for given species and creates output for it
        
        Arguments:
            species {string} -- species which is asked for
        
        Returns:
            fetch [bool] -- true if successful, false if not
        """

        output_path = Path('./out')
        query = "SELECT * FROM fastas WHERE species=?"
        self.cursor.execute(query,(species,))
        output_path = output_path / (species+"_mammalia_fasta.fasta")
        return self.writeToFile(output_path, SEQUENCE)
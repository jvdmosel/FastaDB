from fastas.fasta import Fasta
from fastas.fasta_parser import FastaParser
from fastas.fasta_db import FastaDB
from fastas.entropy_comp import EntropyComp
from pathlib import Path
import sys
import re

def help():
    """Prints available options and valid queries for help
    
    """

    return """
options:
    -help     for valid options and queries
    -exit     to exit program
    -comp     to compare two nodes to another (see examples in report)

valid queries:
    node:     ID seperated by dots
    aligned:  default non-aligned, for aligned ID -a
    species:  full biological name
"""

def success(val, query):
    """Prints if file could be created

    """

    if val:
        print("File was successfully created!")
    elif query[:5] != "-comp":
        print("Error: No data could be fetched for query: " + query)

def has_level_4(length):
    """Returns if level of node is four
    
    """

    return 1 if length == 7 else 0

def file_exists(node):
    """Returns if file for given node exists already

    """

    path = Path('./out/'+node+"_mammalia_aligned_fasta.fasta")
    return path.is_file()

def node_query(db, query):
    """Queries a node to FastaDB object
    
    Arguments:
        db {FastaDB} -- database
        query {string} -- node query
    
    Returns:
        fetch [bool] -- true if success, false if not
    """

    fetch = False
    # check if aligned should be true
    if query[-1] == 'a':
        query = query[:len(query)-len(" -a")]
        fetch = db.get_node(query, True)
    else:
        fetch = db.get_node(query, False)
    return fetch

def spec_query(db, query):
    """Queries a species to FastaDB object
    
    Arguments:
        db {FastaDB} -- database
        query {string} -- node query
    
    Returns:
        fetch [bool] -- true if species could be found, false if not
    """
    species = "_".join(query.split())
    return db.get_species(species)

def comp_query(db, query):
    """Compare query
    
    Arguments:
        db {FastaDB} -- database
        query {string} -- node query
    
    Returns:
        fetch [bool] -- true if success, false if not
    """

    fetch = False
    nodes = query.split(' ')
    a = has_level_4(len(nodes[1]))
    b = has_level_4(len(nodes[2]))
    if not (a and b):
        try:
            if not file_exists(nodes[1]):
                fetch = db.get_node(nodes[1], True)
            if not file_exists(nodes[2]):
                fetch = db.get_node(nodes[2], True)
            comp = EntropyComp(nodes[1]+"_mammalia_aligned_fasta.fasta",nodes[2]+"_mammalia_aligned_fasta.fasta")
            print(comp)
            return fetch
        except IndexError:
            print("Error: Sequence lengths differ")
            return False
        except FileNotFoundError:
            print("Error: No data could be fetched for query: " + query)
            return False
    else:
        print("Error: To compare, atleast one node must be of level 5")
        return fetch

def main():
    """Loops over user input

    """

    db = FastaDB(Path('src/db/fasta.db'))
    print("Creating database...")
    db.build_table()
    print("Successful.")
    # regex to filter out queries
    id_regx = re.compile("^((([0-9]{1,2}\.){1,4}[0-9]{1,2})|^[0-9])(\s-a){0,1}$")
    spec_regx = re.compile("^[A-Z]([a-z]*)(\s([a-z]*)){1,2}$")
    comp_regx = re.compile("^\-comp(\s(([0-9]{1,2}\.){3,4}[0-9]{1,2})){2}$")
    print("Please type in queries as specified by the report. For help type in -help, to exit the program use -exit.")
    try:
        while True:
            query = input()
            if query == "-help":
                print(help())
            elif query == "-exit":
                sys.exit(0)
            elif id_regx.match(query):
                success(node_query(db, query), query)
            elif spec_regx.match(query):
                success(spec_query(db, query), query)
            elif query[:5] == "-comp":
                if comp_regx.match(query):
                    success(comp_query(db, query), query)
                else:
                    print("Error: Node(s) are not of level 4/5")
            else:
                print("Error: " + query + " is not a valid input")
    # catch interrupts, exit
    except (KeyboardInterrupt, SystemExit, EOFError):
        print("\nExiting program")
        sys.exit(0)    

if __name__ == "__main__":
    main()
"""
I want to see what the distance between the ends of the peptide are for a bunch of pMHC complexes. I want to know if the ends should be treated as fixed.

Given: A list of PDB entries, presumably pMHC complexes. I assume the smallest chain in each is the peptide.

Print distance between first and last alpha atom of each peptide.
"""

import sys
import os
from Bio.PDB import *
from kabsch import *

if len(sys.argv) == 3:
    pdb_id = sys.argv[1]
    pdb_l = PDBList()
    parser = PDBParser()

    to_print = list()



    pdb_file_name = pdb_l.retrieve_pdb_file(pdb_id)
    structure = parser.get_structure('B', pdb_file_name)
    
    ppb = PPBuilder()
    #store a list of tuples like (sequence, polypeptide). 
    polypeptides = list()

    for polypeptide in ppb.build_peptides(structure):
        sequence = polypeptide.get_sequence()
        polypeptides.append((sequence, polypeptide))
        sys.stdout = sys.__stdout__
    #peptides are between 6 and 12 residues long. I'm being really generous here
    possible_peptides = list(filter(lambda x: len(x[0]) >= 6 and len(x[0]) <= 12, polypeptides))
    #just print out the peptide we're gonna use.
    to_print.append('Will use peptide: ' + str(possible_peptides[0][0]) + ' for PDB entry: ' + pdb_id)
    pMHC_c_alpha_atoms = possible_peptides[0][1].get_ca_list()


    #now, open up the structure of our new peptide
    new_peptide_struct = parser.get_structure('C', sys.argv[2])
    chain = new_peptide_struct.get_chains()[0]
    print('chain id: ' + str(chain.get_id()))
    new_pp = ppb.build_peptides(new_peptide_struct)[0]
    new_pp_atoms = new_pp.get_ca_list()
    #get the atom positions of the new peptide's first and last residue, run kabsch on it.
    positions = kabsch([pMHC_c_alpha_atoms[0].get_coord(), pMHC_c_alpha_atoms[-1].get_coord()], [new_pp_atoms[0].get_coord(), new_pp_atoms[-1].get_coord()])
    model = structure.get_models()[0]
    model.add(chain)
    
    
    
        
            

    for x in to_print:
        print(x)

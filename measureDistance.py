"""
I want to see what the distance between the ends of the peptide are for a bunch of pMHC complexes. I want to know if the ends should be treated as fixed.

Given: A list of PDB entries, presumably pMHC complexes. I assume the smallest chain in each is the peptide.

Print distance between first and last alpha atom of each peptide.
"""

import sys
import os
from Bio.PDB import *

if len(sys.argv) > 1:
    pdb_ids = sys.argv[1::]
    pdb_l = PDBList()
    parser = PDBParser()
    i = 1
    to_print = list()
    for pdb_id in pdb_ids:
        #I'm going to silence the output of the download
        sys.stdout = open(os.devnull, 'w')
        pdb_file_name = pdb_l.retrieve_pdb_file(pdb_id)
        structure = parser.get_structure('A' + str(i), pdb_file_name)
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
        c_alpha_atoms = possible_peptides[0][1].get_ca_list()
        to_print.append('distance: ' + str(c_alpha_atoms[-1] - c_alpha_atoms[0]))
        
        
            
        i += 1
    for x in to_print:
        print(x)

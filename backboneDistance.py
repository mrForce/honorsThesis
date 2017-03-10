"""
I want to see what the distance between the ends of the peptide are for a bunch of pMHC complexes. I want to know if the ends should be treated as fixed.

Given: A list of PDB entries, presumably pMHC complexes. I assume the smallest chain in each is the peptide.

Print distance between first and last alpha atom of each peptide.
"""

import sys
import os
from statistics import mean
from Bio.PDB import *

if len(sys.argv) > 1:
    pdb_id_file = sys.argv[1]
    pdb_ids = list()
    with open(pdb_id_file, 'r') as f:
        pdb_ids = list(filter(lambda x: len(x)  == 4, f.read().split('\n')))
    print('pdb ids')
    print(pdb_ids)
    pdb_l = PDBList()
    parser = PDBParser()
    i = 1
    bond_lengths_two = list()
    for pdb_id in pdb_ids:
        #I'm going to silence the output of the download
        bond_lengths = list() 
        print('pdb_id')
        print(pdb_id)
#        sys.stdout = open(os.devnull, 'w')
        pdb_file_name = pdb_l.retrieve_pdb_file(pdb_id)
        structure = parser.get_structure('A' + str(i), pdb_file_name)
        ppb = PPBuilder()
        #store a list of tuples like (sequence, polypeptide). 
        polypeptides = list()

        for polypeptide in ppb.build_peptides(structure):
            sequence = polypeptide.get_sequence()
            polypeptides.append((sequence, polypeptide))
       # sys.stdout = sys.__stdout__
        print('polypeptides')
        print(polypeptides)
        #peptides are between 6 and 12 residues long. I'm being really generous here
        possible_peptides = list(filter(lambda x: len(x[0]) >= 6 and len(x[0]) <= 20, polypeptides))
        print('possible peptides')
        print(possible_peptides)
        #just print out the peptide we're gonna use.
        c_alpha_atoms = possible_peptides[0][1].get_ca_list()
        for j in range(0, len(c_alpha_atoms) - 1):
            bond_lengths.append(c_alpha_atoms[j + 1] - c_alpha_atoms[j])
        bond_lengths_two.append(bond_lengths)
        print('min bond length: ' + str(min(bond_lengths)))
        print('max bond length: ' + str(max(bond_lengths)))
#        print('mean bond length: ' + str(mean(bond_lengths)))
        
            
        i += 1
print('bond lengths')
for x in bond_lengths_two:
    print(x)

        

from Bio.PDB import *

parser = PDBParser()
structure = parser.get_structure('thing', '5hgh.pdb')

for model in structure:
    for chain in model:
        print(list(chain.get_residues()))

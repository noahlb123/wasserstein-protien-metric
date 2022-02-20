from Bio.PDB.ic_rebuild import structure_rebuild_test
from Bio.PDB.internal_coords import IC_Chain
from Bio.PDB.PDBParser import PDBParser

#simple script to get dihedral angles
parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure("placeholder", './pdb-files/pdb6s7o.ent')

print(IC_Chain(structure).ordered_aa_ic_list)
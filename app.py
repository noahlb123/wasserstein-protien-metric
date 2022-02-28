pdb_path = './pdb-files/pdb6s7o.ent'

#Use Bio to get phi psi
from Bio.PDB.ic_rebuild import structure_rebuild_test
from Bio.PDB.internal_coords import IC_Chain
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser(PERMISSIVE=1)
structure = parser.get_structure("placeholder", pdb_path)

#phi psi angles
#print(IC_Chain(structure).ordered_aa_ic_list)


#Use RamachanDraw (superset of Bio) to get phi psi and plot
from RamachanDraw import phi_psi, plot

#phi psi dict
phi_psi_dict = phi_psi(pdb_path, return_ignored=False)

#plot
plot(pdb_path, out='ram-plots/plot.png')
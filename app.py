pdb_path = './pdb-files/4axp.pdb'

#Use Bio to get phi psi
from Bio.PDB.ic_rebuild import structure_rebuild_test
from Bio.PDB.internal_coords import IC_Chain
from Bio.PDB.PDBParser import PDBParser

#get structure with Bio
#parser = PDBParser(PERMISSIVE=1)
#structure = parser.get_structure("placeholder", pdb_path)

#phi psi angles
#print(IC_Chain(structure).ordered_aa_ic_list)


#Use RamachanDraw (superset of Bio) to get phi psi and plot
from RamachanDraw import phi_psi, plot

#phi psi dict
phi_psi_dict = phi_psi(pdb_path, return_ignored=False)

#plot with RamachanDraw
plot(['./pdb-files/4axp.pdb', './pdb-files/2ljl.pdb'], out='ram-plots/plot.png')
import matplotlib.pyplot as plt
plt.show()

"""#plot with Pyplot from Matplotlib
phi_psi_pairs = phi_psi_dict.values()

def extractSublist(l, sublist_index):
    return [item[sublist_index] for item in l]

phi = extractSublist(phi_psi_pairs, 0)
psi = extractSublist(phi_psi_pairs, 1)

import matplotlib.pyplot as plt
figure = plt.figure(figsize=(1,1))
ax1 = figure.add_subplot(211)
ax1.scatter(phi, psi, s=1)
ax1.axis('equal')
ax1.set_ylabel('psi angles')
ax1.set_xlabel('phi angles')
ax1.set_title('Ramachandran Plot')
plt.show()
figure.savefig('pyplot.png')"""
print('done.')
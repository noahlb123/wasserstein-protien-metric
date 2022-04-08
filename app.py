#Use Bio to get phi psi
import os
from Bio.PDB.internal_coords import IC_Chain
from Bio.PDB import PDBParser, PPBuilder
from math import pi
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from pkg_resources import resource_stream

pdb_paths = ['./pdb-files/4axp.pdb', './pdb-files/2ljl.pdb']

#get structure with Bio
#parser = PDBParser(PERMISSIVE=1)
#structure = parser.get_structure("placeholder", pdb_path)

#phi psi angles
#print(IC_Chain(structure).ordered_aa_ic_list)

#Use RamachanDraw (superset of Bio) to get phi psi and plot (disabled now bec usinging manualy)
#from RamachanDraw import plot

#from inside the ramdraw package, hopefully my pull will be accepted and this code can be removed
def plot(pdb_file, cmap='viridis', alpha=0.75, dpi=100, save=True, show=False, out='plot.png'):
    batch_mode = [True if type(pdb_file) is list else False][0]

    def get_ignored_res(file: str):
        names, x, y, ignored, output = [], [], [], [], {}
        for model in PDBParser().get_structure(id=None, file=file):
            for chain in model:
                peptides = PPBuilder().build_peptides(chain)
                for peptide in peptides:
                    for aa, angles in zip(peptide, peptide.get_phi_psi_list()):
                        residue = aa.resname + str(aa.id[1])
                        output[residue] = angles

        for key, value in output.items():
            # Only get residues with both phi and psi angles
            if value[0] and value[1]:
                x.append(value[0] * 180 / pi)
                y.append(value[1] * 180 / pi)
                names.append(key)
            else:
                ignored.append((key, value))
        
        return output, ignored, x, y, names

    size = [(8.5, 5) if batch_mode else (5.5, 5)][0]
    figure = plt.figure(figsize=size, dpi=dpi)
    ax = plt.subplot(111)
    ax.set_title("".join(["Batch" if batch_mode else pdb_file]))

    # Import 'density_estimate.dat' data file
    Z = np.fromfile(resource_stream('RamachanDraw', 'data/density_estimate.dat'))
    Z = np.reshape(Z, (100, 100))

    ax.set_aspect('equal')
    ax.set_xlabel('\u03C6')
    ax.set_ylabel('\u03C8')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], minor=False)
    ax.set_yticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], minor=False)
    plt.axhline(y=0, color='k', lw=0.5)
    plt.axvline(x=0, color='k', lw=0.5)
    plt.grid(visible=None, which='major', axis='both', color='k', alpha=0.2)

    # Normalize data
    data = np.log(np.rot90(Z))
    ax.imshow(data, cmap=plt.get_cmap(cmap), extent=[-180, 180, -180, 180], alpha=alpha)

    # Fit contour lines correctly
    data = np.rot90(np.fliplr(Z))
    ax.contour(data, colors='k', linewidths=0.5,
               levels=[10 ** i for i in range(-7, 0)],
               antialiased=True, extent=[-180, 180, -180, 180], alpha=0.65)

    def start(fp, color=None):
        assert os.path.exists(fp), \
            'Unable to fetch file: {}. PDB entry probably does not exist.'.format(pdb_file)
        phi_psi_data, ignored_res, x, y, names = get_ignored_res(file=fp)
        scatter = ax.scatter(x, y, marker='.', s=3, c="".join([color if color else 'k']), label=fp)
        return phi_psi_data, ignored_res, scatter, x, y, names

    if batch_mode:
        file_output_map = {key: None for key in pdb_file}
        for _, file in enumerate(pdb_file):
            file_output_map[file] = start(fp=file, color=list(mcolors.BASE_COLORS.keys())[_])
        ax.legend(bbox_to_anchor=(1.04, 1), loc='upper left')
    else:
        output = start(fp=pdb_file)

    if save:
        plt.savefig(out)
    if show:
        plt.show()
    
    #return params
    if batch_mode:
        return ax, file_output_map, figure
    else:
        return ax, output, figure

#plot with RamachanDraw
ax, output, fig = plot(['./pdb-files/4axp.pdb', './pdb-files/2ljl.pdb'], out='ram-plots/plot.png')
x, y, names, scatters = [], [], [], []
for key in output:
    protein = output[key]
    names += protein[-1]
    y += protein[-2]
    x += protein[-3]
    scatters.append(protein[-4])

#from https://stackoverflow.com/questions/7908636/how-to-add-hovering-annotations-in-matplotlib
annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
norm = plt.Normalize(1,4)
c = np.random.randint(1,5,size=15)

def update_annot(ind, annot, sc, names):
    cmap = plt.cm.RdYlGn
    pos = sc.get_offsets()[ind["ind"][0]]
    annot.xy = pos
    text = "{}, {}".format(" ".join(list(map(str,ind["ind"]))), 
                           " ".join([names[n] for n in ind["ind"]]))
    annot.set_text(text)
    annot.get_bbox_patch().set_facecolor('#0F0F0F')
    annot.get_bbox_patch().set_alpha(0.4)

#ax might be defined wrong
def hover(event):
    sc = scatters[0]
    vis = annot.get_visible()
    if event.inaxes == ax:
        cont, ind = sc.contains(event)
        if cont:
            update_annot(ind, annot, sc, names)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", hover)
plt.show()

import sys
import os
import warnings
import matplotlib.pyplot as plt
from pymatgen.core.ion import Ion
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, IonEntry, PourbaixDiagram
from pymatgen.analysis.pourbaix_diagram import PourbaixPlotter
from pymatgen.entries.computed_entries import ComputedEntry

plotname = 'Volume_pourbaix_FeOx_pymatgen'
API_KEY = os.getenv('MAPI_KEY')
if not API_KEY:
    sys.exit("Error: MAPI_KEY environment variable not set.")
mpr = MPRester(API_KEY)

entries = mpr.get_entries_in_chemsys(['O', 'H'])

# ion_dict_Co = mpr.get_pourbaix_entries(['Co'])
# ion_dict_Fe = mpr.get_pourbaix_entries(['Fe'])
ion_dict_Fe = mpr.get_pourbaix_entries(['Fe', 'N', 'C'])

# Generate phase diagram
pd = PhaseDiagram(entries)
stable_solids = pd.stable_entries
stable_solids_minus_h2o = [entry for entry in stable_solids if
                           entry.composition.reduced_formula not in ["H2", "O2", "H2O", "H2O2"]]

# Create Pourbaix entries for solids
pbx_solid_entries = []
kJmol = 96.485
ion_dict_solids_expt = {
    'Fe': 0,
    'Fe2O3': -743.8 / kJmol,
    'Fe3O4': -1015.7 / kJmol,
    'FeOOH': -496 / kJmol
}

for entry in stable_solids_minus_h2o:
    pbx_entry = PourbaixEntry(entry)
    pbx_entry.g0_replace(pd.get_form_energy(entry))
    pbx_solid_entries.append(pbx_entry)

for key, energy in ion_dict_solids_expt.items():
    comp = Ion.from_formula(key)
    pbx_entry_ion = PourbaixEntry(ComputedEntry(comp, energy))
    pbx_entry_ion.entry_id = key
    pbx_entry_ion.conc = 1
    pbx_solid_entries.append(pbx_entry_ion)

all_entries = ion_dict_Fe

def plot_pourbaix(entries):
    """Plot and save Pourbaix diagram."""
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    plotter = PourbaixPlotter(pourbaix)

    # Generate the Pourbaix plot and get the Axes object
    ax = plotter.get_pourbaix_plot(limits=[[-2, 16], [-2, 4]])

    # Customize the plot
    for line in ax.lines:
        line.set_linewidth(1.0)  # Adjust line thickness

    for text in ax.texts:
        text.set_fontsize(14)  # Adjust phase label font size

    # Set axis labels and tick font size
    ax.set_xlabel("pH", fontsize=14)
    ax.set_ylabel("Potential (V vs SHE)", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)

    fig = ax.figure
    fig.set_size_inches((8, 7))

    plt.savefig(plotname + '.png', dpi=100, bbox_inches='tight')
    plt.savefig(plotname + '.pdf', dpi=100, bbox_inches='tight')
    plt.show()

plot_pourbaix(all_entries)

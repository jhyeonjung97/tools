from sys import argv
import numpy as np
from ase.io import read, write
from ase.data import covalent_radii as cradii
from ase.geometry import get_distances
from ase.visualize import view
from catkit.gen.utils.connectivity import get_cutoff_neighbors


def set_formal_oxidation_state(atoms, charge_O=-2, charge_N=-3,
                               charge_H=1, charge_C=-4, bond_scale=False):
    O_indices = np.array([i for i, a in enumerate(atoms)
                          if a.symbol == 'O'])
    N_indices = np.array([i for i, a in enumerate(atoms)
                      if a.symbol == 'N'])
    C_indices = np.array([i for i, a in enumerate(atoms)
                      if a.symbol == 'C'])
    H_indices = np.array([i for i, a in enumerate(atoms)
                          if a.symbol == 'H'])
    M_indices = np.array([i for i, a in enumerate(atoms)
                          if not a.symbol in ['O', 'H', 'N', 'C']])
    print(M_indices)
    # Connectivity matrix from CatKit
    con_matrix = get_cutoff_neighbors(atoms, scale_cov_radii=1.5)
    oxi_states = np.ones([len(atoms)])
    if len(O_indices) > 0:
        oxi_states[O_indices] = charge_O
    if len(N_indices) > 0:
        oxi_states[N_indices] = charge_N
    if len(C_indices) > 0:
        oxi_states[C_indices] = charge_C
    if len(H_indices) > 0:  # First correct O charge due to H
        oxi_states[H_indices] = charge_H
        for H_i in H_indices:
            H_O_connectivity = con_matrix[H_i][O_indices]
            norm = np.sum(H_O_connectivity)
            O_indices_H = O_indices[np.where(H_O_connectivity)[0]]
            oxi_states[O_indices_H] += charge_H / norm
    for metal_i in M_indices:  # Substract O connectivity
        M_O_connectivity = con_matrix[metal_i][O_indices]
        norm = np.sum(con_matrix[O_indices][:, M_indices], axis=-1)
        print(norm)
        oxi_states[metal_i] = sum(
            M_O_connectivity * -oxi_states[O_indices] / norm)
    atoms.set_initial_charges(np.round(oxi_states, 4))
    write(argv[2], atoms)
    return atoms, oxi_states

if __name__ == "__main__":
    atoms = read(argv[1])
    atoms, ox = set_formal_oxidation_state(atoms)
    print('oxidation states: ', ox)
    # view(atoms)  # oxidation states written to initial charges
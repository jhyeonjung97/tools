import subprocess
import numpy as np
from ase.io import read
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

atoms = read('restart.json')

ldau_luj = {
    'Mn': {'L': 2, 'U': 2.75, 'J': 0.0},
    'O': {'L': -1, 'U': 0.0, 'J': 0.0},
    'H': {'L': -1, 'U': 0.0, 'J': 0.0}
    }

def get_kpoints(atoms, effective_length=25, bulk=False):
    """
    Return a tuple of k-points derived from the unit cell.
    """
    l = effective_length
    cell = atoms.get_cell()
    nkx = int(round(l/np.linalg.norm(cell[0]),0))
    nky = int(round(l/np.linalg.norm(cell[1]),0))
    if bulk == True:
        nkz = int(round(l/np.linalg.norm(cell[2]),0))
    else:
        nkz = 1
    return((nkx, nky, nkz))

kpoints = get_kpoints(atoms, effective_length=25, bulk=True)

atoms.calc = vasp_calculator.Vasp(
    inimix=0,
    amix=0.1,
    amix_mag=0.1,
    bmix=1e-05,
    bmix_mag=1e-05,
    encut=600,
    sigma=0.05,
    ediff=1e-05,
    ediffg=-0.02,
    algo='Normal',
    gga='PE',
    prec='Normal',
    ibrion=2,
    isif=3,
    ismear=0,
    ispin=2,
    ldauprint=2,
    ldautype=2,
    lmaxmix=4,
    lorbit=11,
    nelm=250,
    npar=6,
    nsw=199,
    laechg=True,
    lasph=True,
    ldau=True,
    lvtot=False,
    lreal='Auto',
    ldau_luj=ldau_luj,
    kpts=kpoints,
    gamma=True,
)

energy = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

Trajectory(f'final_with_calculator.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_with_calculator.traj final_with_calculator.json', shell=True)
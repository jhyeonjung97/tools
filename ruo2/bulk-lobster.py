import subprocess
from ase.io import read
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

atoms = read('restart.json')

nbands = 0
for atom in atoms:
    if atom.symbol == 'O':
        nbands += 2
    else:
        nbands += 6

ldau_luj = {'Ti':{'L':2, 'U':3.00, 'J':0.0},
            'V': {'L':2, 'U':3.25, 'J':0.0},
            'Cr':{'L':2, 'U':3.50, 'J':0.0},
            'Mn':{'L':2, 'U':3.75, 'J':0.0},
            'Fe':{'L':2, 'U':4.30, 'J':0.0},
            'Co':{'L':2, 'U':3.32, 'J':0.0},
            'Ni':{'L':2, 'U':6.45, 'J':0.0},
            'Cu':{'L':2, 'U':9.00, 'J':0.0},
            'Ru': {'L':2, 'U':2.5, 'J': 0.0}
            }

lmaxmix = 2
for atom in atoms:
    if atom.symbol in ldau_luj:
        lmaxmix = 4
    else:
        ldau_luj[atom.symbol] = {'L': -1, 'U': 0.0, 'J': 0.0}

atoms.calc = vasp_calculator.Vasp(
    #istart=1,
    inimix=0,
    amix=0.1,
    amix_mag=0.1,
    bmix=0.0001,
    bmix_mag=0.0001,
    gga='PE',
    lreal='Auto',
    algo='Normal',
    prec='Normal',
    #kpar=8,
    ncore=4,
    ismear=0,
    sigma=0.05,
    encut=500,
    ediff=1e-06,
    ediffg=-0.01,
    isym=0,
    isif=2,
    ispin=2,
    ibrion=-1,
    nbands=nbands,
    nelm=250,
    nsw=0,
    gamma=True,
    kpts=[5,5,5],
    ldau=True,
    lmaxmix=lmaxmix,
    ldau_luj=ldau_luj,
    ldautype=2,
    ldauprint=2,
    lorbit=11,
    lasph=True,
    lvtot=False,
    laechg=True,
    )

energy = atoms.get_potential_energy()
Trajectory(f'final_with_calculator.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_with_calculator.traj final_with_calculator.json', shell=True)

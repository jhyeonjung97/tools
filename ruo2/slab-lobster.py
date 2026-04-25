import subprocess
from ase.io import read
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

atoms = read('restart.json')

nbands = 2
for atom in atoms:
    if atom.symbol == 'O':
        nbands += 4
    elif atom.symbol in ['Sc', 'Y', 'Zr', 'Nb', 'La']:
        nbands += 10
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
            'Ru':{'L':2, 'U':2.50, 'J':0.0}
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
    bix_mag=0.0001,
    gga='PE',
    lreal='Auto',
    algo='Normal',
    prec='Normal',
    kpar=9,
    ismear=0,
    sigma=0.05,
    encut=500,
    ediff=1e-05,
    ediffg=-0.02,
    isym=0,
    isif=2,
    ispin=2,
    ibrion=-1,
    nbands=nbands,
    nsw=0,
    gamma=True,
    kpts=[4,4,1],
    ldau=True,
    lmaxmix=lmaxmix,
    ldau_luj=ldau_luj,
    ldautype=2,
    ldauprint=2,
    lorbit=11,
    laechg=True,
    lasph=True,
    lvtot=False,
    #lsol=True,
    )

energy = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

Trajectory(f'final_with_calculator.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_with_calculator.traj final_with_calculator.json', shell=True)
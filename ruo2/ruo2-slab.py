import subprocess
from ase.io import read
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

atoms = read('restart.json')

ldau_luj = {
    'Ru': {'L': 2, 'U': 0.0, 'J': 0.0},
    'Re': {'L': -1, 'U': 0.0, 'J': 0.0},
    'O': {'L': -1, 'U': 0.0, 'J': 0.0},
    'H': {'L': -1, 'U': 0.0, 'J': 0.0}
    }

atoms.calc = vasp_calculator.Vasp(
     inimix=0,
     amix=0.1,
     amix_mag=0.1,
     bmix=0.0001,
     bmix_mag=0.0001,
     encut=500,
     sigma=0.05,
     ediff=1e-05,
     ediffg=-0.02,
     algo='Fast',
     gga='PE',
     prec='Normal',
     ibrion=2,
     isif=2,
     ismear=0,
     ispin=2,
     istart=1,
     ldau=True,
     lmaxmix=4,
     ldau_luj=ldau_luj,
     ldautype=2,
     ldauprint=2,
     lorbit=11,
     nelm=250,
     npar=6,
     nsw=199,
     laechg=True,
     lasph=True,
     lvtot=False,
     lreal='Auto',
     kpts=[5, 5, 1],
     gamma=True,
     lsol=True,
    )

energy = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

Trajectory(f'final_with_calculator.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_with_calculator.traj final_with_calculator.json', shell=True)
import subprocess
from ase.io import read
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

name = 'mnc1'

ldau_luj = {'Ti': {'L':2, 'U':3.00, 'J':0.0},
            'V': {'L':2, 'U':3.25, 'J':0.0},
            'Cr': {'L':2, 'U':3.50, 'J':0.0},
            'Mn': {'L':2, 'U':3.75, 'J':0.0},
            'Fe': {'L':2, 'U':4.30, 'J':0.0},
            'Co': {'L':2, 'U':3.32, 'J':0.0},
            'Ni': {'L':2, 'U':6.45, 'J':0.0},
            'Cu': {'L':2, 'U':3.00, 'J':0.0},
            'C': {'L':-1, 'U':0.0, 'J':0.0},
            'N': {'L':-1, 'U':0.0, 'J':0.0},
            'O': {'L':-1, 'U':0.0, 'J':0.0},
            'H': {'L':-1, 'U':0.0, 'J':0.0}
            }

atoms = read('restart.json')
        
atoms.calc = vasp_calculator.Vasp(
                    encut=500,
                    gga='PE',
                    ivdw=12,
                    kpts=(5,5,1),
                    kpar=8,
                    npar=4,
                    gamma=True,
                    ismear=0,
                    sigma=0.05,
                    nelm=250,
                    algo='Normal',
                    ibrion=2,
                    isif=2,
                    ediffg=-0.02,
                    ediff=1e-5,
                    prec='Normal',
                    nsw=200,
                    lvhar=True,
                    lvtot=False,
                    ispin=2,
                    setups={'base': 'recommended', 'W': '_sv'},
                    ldau=True,
                    ldautype=2,
                    ldau_luj=ldau_luj,
                    ldauprint=2,
                    lmaxmix=4,
                    lasph=True,
                    laechg=True,
                    lreal='Auto',
                    nedos=3000,
                    lorbit=11,
                    lsol=True,
                    nupdown=1
                    )

energy = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

Trajectory(f'final_{name}.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)
import subprocess
from os import path
from ase.io import read
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

name = 'gas'

if path.exists('restart.json'):
    atoms = read('restart.json')
else:
    raise ValueError('No restart.json file found')

atoms.calc = vasp_calculator.Vasp(
                    encut=500,
                    xc='PBE',
                    gga='PE',
                    kpts=(1,1,1),
                    gamma=True,
                    ismear=0,
                    sigma=0.05,
                    # ivdw=12,
                    # lsol=True,
                    nelm=250,
                    algo='Normal',
                    ibrion=2,
                    isif=2,
                    ediffg=-0.02,
                    ediff=1e-5,
                    prec='Normal',
                    nsw=200,
                    ispin=2,
                    nupdown=2,
                    setups={'base': 'recommended'},
                    lreal='False',
                    )

energy = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

Trajectory(f'final_{name}.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_{name}.traj final_with_calculator.json', shell=True)
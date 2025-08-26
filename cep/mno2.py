import subprocess
from ase.io import read
from ase.io.trajectory import Trajectory

atoms = read('k.json')
mno2_atoms = read('mno2.json')

atoms.calc = mno2_atoms.calc

energy = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

Trajectory(f'final_with_calculator.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_with_calculator.traj final_with_calculator.json', shell=True)
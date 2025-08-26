import subprocess
from ase.io import read
from ase.io.trajectory import Trajectory
import ase.calculators.vasp as vasp_calculator

atoms = read('k.json')

atoms.calc = vasp_calculator.Vasp(
     amix=0.1, 
     amix_mag=0.1, 
     bmix=0.0001, 
     bmix_mag=0.0001, 
     encut=500, 
     sigma=0.05, 
     ediff=1e-05, 
     ediffg=-0.02, 
     algo=normal, 
     gga=PE, 
     prec=Normal, 
     ibrion=2, 
     isif=2, 
     ismear=0, 
     ispin=2, 
     istart=1, 
     kpar=10, 
     ldauprint=2, 
     ldautype=2, 
     lmaxmix=6, 
     lorbit=11, 
     nelm=250, 
     npar=6, 
     nsw=199, 
     inimix=0, 
     laechg=true, 
     lasph=true, 
     ldau=true, 
     lvtot=false, 
     lreal=Auto, 
     ldau_luj={Ti={L=2, U=3.0, J=0.0}, V={L=2, U=3.25, J=0.0}, Cr={L=2, U=3.5, J=0.0}, Mn={L=2, U=2.75, J=0.0}, Fe={L=2, U=4.3, J=0.0}, Co={L=2, U=3.32, J=0.0}, Ni={L=2, U=6.45, J=0.0}, Cu={L=2, U=3.0, J=0.0}, Mo={L=2, U=4.38, J=0.0}, W={L=2, U=6.2, J=0.0}, Ce={L=3, U=4.5, J=0.0}, O={L=-1, U=0.0, J=0.0}, C={L=-1, U=0.0, J=0.0}, Au={L=-1, U=0.0, J=0.0}, Ir={L=-1, U=0.0, J=0.0}, H={L=-1, U=0.0, J=0.0}}, xc=PBE, pp=PBE, setups={base=minimal}, 
     txt=-, 
     kpts=[3, 2, 1], 
     gamma=true, 
     reciprocal=false, 
     ignore_constraints=false    
    )

energy = atoms.get_potential_energy()
print('Calculation Complete, storing the run + calculator to traj file')
subprocess.call('sh ~/bin/verve/correct-contcar.sh', shell=True)

Trajectory(f'final_with_calculator.traj', 'w').write(atoms)
subprocess.call(f'ase convert -f final_with_calculator.traj final_with_calculator.json', shell=True)
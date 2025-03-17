
# Specify the charge states
fraction_charges = [-0.5,-1.0,-1.5,-2.0,-2.5,-3.0] ## Change this to the desired charge states

# Perform the calculations for each charge state
for fraction in fraction_charges:
   
   # Calculate the number of electrons for the charge state
   nelect = (total_electrons + fraction) 

   # Create a directory for the charge state and copy the necessary files
   folder_name = 'nelect_{}'.format(nelect)
   folder_path = os.path.join(parent_dir, folder_name)
   os.makedirs(folder_path)
   shutil.copy(restart_path, folder_path)  
   shutil.copy(wavecar_path, folder_path)
   shutil.copy(run_vasp_path, folder_path)

   # Change to the directory for the charge state and set up the calculation
   os.chdir(folder_path) 
   name = 'test_relax'
   new_atoms=read('restart.json')
   new_atoms.write(name+'_init'+'.traj')
   new_atoms.write(name+'_init'+'.cif')

   # Run the charged calculation
   calc.set(istart=1, nelect=nelect)
   new_atoms.set_calculator(calc)
   new_atoms.get_potential_energy()
   new_atoms.get_forces()

   # Write the final structure of the calculation
   traj2=Trajectory('final_with_calculator.traj',  'w')
   traj2.write(new_atoms)
   subprocess.call('ase convert -f final_with_calculator.traj  final_with_calculator.json', shell=True)
   subprocess.call('ase convert -f final_with_calculator.json restart.json', shell=True)
   subprocess.call('ase convert -f OUTCAR full_relax.json', shell=True)
   subprocess.call('/global/cfs/cdirs/m2997/bin/get_restart4', shell=True)

   # Change back to the parent directory
   os.chdir('../')

# Create a directory for the uncharged calculation and copy the OUTCAR file
folder_name = 'nelect_{}'.format(total_electrons)
folder_path = os.path.join(parent_dir, folder_name)
os.makedirs(folder_path)
shutil.copy(outcar_path, folder_path)

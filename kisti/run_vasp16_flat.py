import os
exitcode = os.system('mpiexec -np 8 numactl -p 1 /home01/x2755a09/KISTI_VASP/VASP6/vasp6_vtst_beef_sol')
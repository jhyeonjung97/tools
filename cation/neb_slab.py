"""
Slab NEB (VASP + ASE NEB). Five images: endpoints from 00/restart.json and 04/restart.json;
internal images run in 01–03. See https://ase-lib.org/ase/neb.html
"""
import os
import time
import subprocess
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.mep import NEB, NEBTools
from ase.optimize import FIRE
import ase.calculators.vasp as vasp_calculator

name = 'neb_slab'
start_time = time.time()

# --- NEB: five images (endpoints 00/04, internal 01–03) ---
N_IMAGES = 5
N_INTERNAL = N_IMAGES - 2
ENDPOINT_INITIAL = os.path.join('00', 'restart.json')
ENDPOINT_FINAL = os.path.join('04', 'restart.json')
NEB_K = 0.1
CLIMB = False
FMAX = 0.05
NEB_METHOD = 'improvedtangent'
USE_IDPP = False
NEB_PARALLEL = False
# Use minimum-image convention for interpolation (periodic slab xy)
INTERPOLATE_MIC = True

# Restart: full band from a single trajectory written in a previous run
RESTART_TRAJ = 'neb_restart.traj'


def _load_endpoints():
    if not os.path.isfile(ENDPOINT_INITIAL):
        raise SystemExit(f'Missing initial structure: {ENDPOINT_INITIAL}')
    if not os.path.isfile(ENDPOINT_FINAL):
        raise SystemExit(f'Missing final structure: {ENDPOINT_FINAL}')
    return read(ENDPOINT_INITIAL).copy(), read(ENDPOINT_FINAL).copy()


def _vasp_calc(directory):
    return vasp_calculator.Vasp(
        directory=directory,
        encut=500,
        xc='PBE',
        gga='RP',
        ivdw=12,
        kpts=(3, 3, 1),
        kpar=5,
        ncore=6,
        gamma=True,
        ismear=0,
        sigma=0.05,
        algo='fast',
        lreal='auto',
        ibrion=-1,
        isif=2,
        ispin=2,
        ediff=1e-5,
        nsw=0,
        nelm=600,
        setups='recommended',
        lorbit=11,
    )


def _write_band(path, band):
    w = Trajectory(path, 'w')
    for atoms in band:
        w.write(atoms)
    w.close()


if os.path.exists(RESTART_TRAJ):
    images = list(read(RESTART_TRAJ, index=':'))
    print(f'Restart: loaded {len(images)} images from {RESTART_TRAJ}')
else:
    initial, final = _load_endpoints()
    images = [initial]
    images += [initial.copy() for _ in range(N_INTERNAL)]
    images += [final]
    neb_setup = NEB(images, k=NEB_K, climb=False, method=NEB_METHOD)
    neb_setup.interpolate(mic=INTERPOLATE_MIC)
    if USE_IDPP:
        from ase.mep.neb import idpp_interpolate

        idpp_interpolate(images, steps=100)

neb = NEB(
    images,
    k=NEB_K,
    climb=CLIMB,
    parallel=NEB_PARALLEL,
    method=NEB_METHOD,
)

for i in range(1, len(images) - 1):
    images[i].calc = _vasp_calc(directory=f'{i:02d}')

# Prefer FIRE / MDMin / BFGS over BFGSLineSearch (default QuasiNewton) for NEB
opt = FIRE(neb, trajectory=f'{name}.traj', logfile=f'{name}.log')
opt.run(fmax=FMAX)

_write_band(RESTART_TRAJ, images)
_write_band(f'final_{name}.traj', images)

try:
    tools = NEBTools(images)
    barrier, delta_e = tools.get_barrier()
    print(f'NEB barrier (fit): {barrier:.4f} eV, ΔE: {delta_e:.4f} eV')
except Exception as e:
    print(f'NEBTools analysis skipped: {e}')

subprocess.call(
    f'ase convert -f final_{name}.traj final_with_calculator.json',
    shell=True,
)

end_time = time.time()
elapsed_time = end_time - start_time

with open('time.log', 'a') as f:
    f.write(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
    f.write(f"End Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
    f.write(f"Elapsed Time: {elapsed_time:.2f} seconds\n")
    f.write("=" * 40 + "\n")

print(f"Execution time logged in 'time.log'.")

"""
Slab NEB (VASP + ASE NEB), optimizer BFGS. Five images: endpoints from
00/final_with_calculator.json and 04/final_with_calculator.json; internal 01–03.
Avoid BFGSLineSearch / QuasiNewton default for NEB. See https://ase-lib.org/ase/neb.html
"""
import os
import time
import subprocess
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.mep import NEB, NEBTools
from ase.optimize import BFGS
import ase.calculators.vasp as vasp_calculator

name = 'neb_BFGS'
start_time = time.time()

# --- NEB: five images (endpoints 00/04, internal 01–03) ---
N_IMAGES = 5
N_INTERNAL = N_IMAGES - 2
ENDPOINT_INITIAL = os.path.join('00', 'final_with_calculator.json')
ENDPOINT_FINAL = os.path.join('04', 'final_with_calculator.json')
# ASE NEB spring constant (eV/Å); ASE examples often use ~0.1.
NEB_K = 0.1
CLIMB = True
FMAX = 0.05
NEB_METHOD = 'improvedtangent'
USE_IDPP = False
# True: launch Python with MPI — one rank per internal image, e.g.
#   mpirun -np 3 python neb_BFGS.py
# Each rank runs VASP in its own directory (01/02/03). In run_vasp.py use
# mpirun -np $VASP_MPI_CORES with VASP_MPI_CORES = (node cores) / 3 so three
# VASPs do not each grab the full allocation.
NEB_PARALLEL = False
# Use minimum-image convention for interpolation (periodic slab xy)
INTERPOLATE_MIC = True

# Restart: full band from a single trajectory written in a previous run
RESTART_TRAJ = 'neb_restart.traj'

if NEB_PARALLEL:
    from ase.parallel import world as _mpi_world
else:

    class _SingleWorld:
        rank = 0
        size = 1

    _mpi_world = _SingleWorld()


def _load_endpoints():
    if not os.path.isfile(ENDPOINT_INITIAL):
        raise SystemExit(f'Missing initial structure: {ENDPOINT_INITIAL}')
    if not os.path.isfile(ENDPOINT_FINAL):
        raise SystemExit(f'Missing final structure: {ENDPOINT_FINAL}')
    return read(ENDPOINT_INITIAL), read(ENDPOINT_FINAL)


def _vasp_calc(directory):
    return vasp_calculator.Vasp(
        directory=directory,
        encut=500,
        xc='PBE',
        gga='RP',
        ivdw=12,
        kpts=(3, 3, 1),
        kpar=5,
        ncore=4,
        gamma=True,
        ismear=0,
        sigma=0.05,
        algo='normal',
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
    if _mpi_world.rank == 0:
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

if images[0].calc is None or images[-1].calc is None:
    raise SystemExit(
        'NEB endpoints need calculators (or stored SinglePoint results). '
        'Make sure 00/ and 04/ final_with_calculator.json keep calculator '
        'results; do not strip them with Atoms.copy().'
    )

neb = NEB(
    images,
    k=NEB_K,
    climb=CLIMB,
    parallel=NEB_PARALLEL,
    method=NEB_METHOD,
)

if NEB_PARALLEL:
    if _mpi_world.size != N_INTERNAL:
        raise SystemExit(
            'NEB_PARALLEL with VASP: set MPI size equal to the number of '
            f'internal images ({N_INTERNAL}), e.g. mpirun -np {N_INTERNAL} '
            f'... python neb_BFGS.py (got {_mpi_world.size} ranks).'
        )
    for k, image in enumerate(images[1:-1]):
        if k == _mpi_world.rank:
            image.calc = _vasp_calc(directory=f'{k + 1:02d}')
else:
    for i in range(1, len(images) - 1):
        images[i].calc = _vasp_calc(directory=f'{i:02d}')

_traj = f'{name}.traj' if _mpi_world.rank == 0 else None
_log = f'{name}.log' if _mpi_world.rank == 0 else None
opt = BFGS(neb, trajectory=_traj, logfile=_log)
opt.run(fmax=FMAX)

if _mpi_world.rank == 0:
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

if _mpi_world.rank == 0:
    with open('time.log', 'a') as f:
        f.write(f"Start Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}\n")
        f.write(f"End Time: {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}\n")
        f.write(f"Elapsed Time: {elapsed_time:.2f} seconds\n")
        f.write("=" * 40 + "\n")

    print("Execution time logged in 'time.log'.")

cp ../POSCAR POSCAR.phon
cp ../OUTCAR OUTCAR.phon
cp ../INCAR INCAR
cp ../POTCAR POTCAR
cp ../KPOINTS KPOINTS
cp ../run_slurm.sh run_slurm.sh

sed -i -e 's/IBRION/#IBRION/' \
-e 's/POTIM/#POTIM/' \
-e 's/NSW/#NSW/' \
-e '/LREAL/ a LEPSILON = .TRUE.' \
INCAR
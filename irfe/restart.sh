for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/5_M_top
do
    cd $dir
    find . ! -name '*top*.vasp' -delete
    cp *.vasp POSCAR
done

for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/6_M_hol
do
    cd $dir
    find . ! -name '*hol*.vasp' -delete
    cp *.vasp POSCAR
done
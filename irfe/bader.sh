for dir in /home/hyeonjung/scratch/4_IrFe3/*_*/*_*/*_*_*/
do
    cd $dir
    python ~/bin/verve/bader.py
done
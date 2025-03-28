if [ ${here} == 'slac' ]; then
    root=/Users/jiuy97/Desktop/9_pourbaixGC
fi

for dir in $root/*_vac*/
do
    cd $dir
    python ~/bin/tools/cep/analyzeGC.py --no-dz
    python ~/bin/tools/cep/analyzeGC.py --no-dz --limited
done

for dir in $root/*_*/*s*/
do
    cd $dir
    python ~/bin/tools/cep/analyzeGC.py
    python ~/bin/tools/cep/analyzeGC.py --limited
done

for dir in $root/*_*/
do
    if [[ $dir == *"vac"* ]]; then
        continue
    fi
    cd $dir
    python ~/bin/tools/cep/analyzeGC_spin.py
done
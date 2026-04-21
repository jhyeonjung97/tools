3d_dopants=('Sc' 'Ti' 'V' 'Cr' 'Mn' 'Fe' 'Co' 'Ni' 'Cu' 'Zn' 'Ga' 'Ge')
4d_dopants=('Y' 'Zr' 'Nb' 'Mo' 'Tc' 'Ru' 'Rh' 'Pd' 'Ag' 'Cd' 'In' 'Sn')
5d_dopants=('La' 'Hf' 'Ta' 'W' 'Re' 'Os' 'Ir' 'Pt' 'Au' 'Hg' 'Tl' 'Pb')

for dir in */; do
    IFS='_' read -r num dopant <<< "$dir"
    if [ $dopant in ${3d_dopants[@]} ]; then
        sed -i -e "/basisFunctions X/c\basisFunctions ${dopant} 3d 4s" $dir/lobsterin
    elif [ $dopant in ${4d_dopants[@]} ]; then
        sed -i -e "/basisFunctions X/c\basisFunctions ${dopant} 4d 5s" $dir/lobsterin
    elif [ $dopant in ${5d_dopants[@]} ]; then
        sed -i -e "/basisFunctions X/c\basisFunctions ${dopant} 5d 6s" $dir/lobsterin
    fi
    sed -i -e "s/X/${dopant}/" $dir/lobsterin
done
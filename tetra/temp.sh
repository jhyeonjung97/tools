for dir in ./*/
do
    cd $dir
    if [[ $PWD == *'Tetrahedral'* ]]; then
        n=4; python ~/bin/tools/tetra/energy.py --save -p hexa -x "Metal (MO)" -y "Hexagonal ratio [c/a]"
    elif [[ $PWD == *'Pyramidal'* ]] || [[ $PWD == *'Tetrahedral'* ]] || [[ $PWD == *'SquarePlanar'* ]]; then
        n=4; #python ~/bin/tools/tetra/energy.py --save -p hexa -x "Metal (MO)" -y "Square prism ratio [c/a]"
    elif [[ $PWD == *'Octahedral'* ]]; then
        n=6
    fi
    python ~/bin/tools/tetra/energy.py --save -p energy -x "Metal (MO)" -y "Total energy (eV)"
    python ~/bin/tools/tetra/energy.py --save -p energy -x "Metal (MO)" -y "Total energy (eV/MO)" -n m
    python ~/bin/tools/tetra/energy.py --save -p bond  -x "Metal (MO)" -y "Bond length (A)"
    python ~/bin/tools/tetra/energy.py --save -p bond -x "Metal (MO)" -y "Bond length (A/M-O)" -n $n
    python ~/bin/tools/tetra/energy.py --save -p volume -x "Metal (MO)" -y "Volume (A^3/MO)" -n m
    python ~/bin/tools/tetra/energy.py --save -p chg -e M  -x "Metal (MO)" -y "Bader charge (e-)"
    python ~/bin/tools/tetra/energy.py --save -p mag -e M -x "Metal (MO)" -y "|Magnetization|"
    python ~/bin/tools/tetra/energy.py --save -p ICOHP -x "Metal (MO)" -y "ICOHP (eV/MO)"
    python ~/bin/tools/tetra/energy.py --save -p ICOHP -x "Metal (MO)" -y "ICOHP (eV/M-O)" -n $n
    python ~/bin/tools/tetra/energy.py --save -p ICOBI -x "Metal (MO)" -y "ICOBI (/MO)"
    python ~/bin/tools/tetra/energy.py --save -p ICOBI -x "Metal (MO)" -y "ICOBI (eV/M-O)" -n $n
    python ~/bin/tools/tetra/energy.py --save -p GP_L -e M  -x "Metal (MO)" -y "Gross population (Loewdin)"
    python ~/bin/tools/tetra/energy.py --save -p Madelung_L -x "Metal (MO)" -y "Madelugn energy (Loewdin, eV/MO)" -n m
    python ~/bin/tools/tetra/formation_energy.py
    python ~/bin/tools/tetra/cohesive_energy.py
    # python ~/bin/tools/tetra/energy.py --save -p PSCENC -x "Metal (MO)" -y "PSCENC (eV/MO)" -n m
    # python ~/bin/tools/tetra/energy.py --save -p TEWEN -x "Metal (MO)" -y "TEWEN (eV/MO)" -n m
    # python ~/bin/tools/tetra/energy.py --save -p DENC -x "Metal (MO)" -y "DENC (eV/MO)" -n m
    # python ~/bin/tools/tetra/energy.py --save -p EXHF -x "Metal (MO)" -y "EXHF (eV/MO)" -n m
    # python ~/bin/tools/tetra/energy.py --save -p XCENC -x "Metal (MO)" -y "XCENC (eV/MO)" -n m
    # python ~/bin/tools/tetra/energy.py --save -p PAW_double_counting -x "Metal (MO)" -y "PAW_double_counting (eV/MO)" -n m
    # python ~/bin/tools/tetra/energy.py --save -p EENTRO -x "Metal (MO)" -y "EENTRO (eV/MO)" -n m
    # python ~/bin/tools/tetra/energy.py --save -p EBANDS -x "Metal (MO)" -y "EBANDS (eV/MO)" -n m
    # python ~/bin/tools/tetra/energy.py --save -p EATOM -x "Metal (MO)" -y "EATOM (eV/MO)" -n m
    sed -i 's/\x0//g' *.tsv
done
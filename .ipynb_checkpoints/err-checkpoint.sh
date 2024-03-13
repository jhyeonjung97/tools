if [[ $1 === '-r' ]]; then
    grep '\-\-\-\-\-\-\-\-\-\-\-\-' */vasp.out
else
    grep '\-\-\-\-\-\-\-\-\-\-\-\-' vasp.out

fi
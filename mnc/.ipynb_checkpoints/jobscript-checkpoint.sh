for dir in /scratch/x2755a09/3_MNC/*d/*_*/*_*S/*_/; do
    cd $dir
    IFS='/' read -r -a path_components <<< $PWD
    metal=$(echo "${path_components[-3]}" | cut -d'_' -f2)
    spin=$(echo "${path_components[-2]}" | cut -d'_' -f2)
    dz=$(echo "${path_components[-1]}" | cut -d'_' -f1)
    sed -i -e "/#PBS -N/c\#PBS -N $metal$spin$dz" submit.sh
    if [[ $spin == 'LS' ]]; then
       sed -i -e 's/mnc-sol-ls.py/mnc-sol-ls.py/g' submit.sh
       sed -i -e 's/mnc-sol-is.py/mnc-sol-ls.py/g' submit.sh
       sed -i -e 's/mnc-sol-hs.py/mnc-sol-ls.py/g' submit.sh
       sed -i -e 's/mnc-sol-ns.py/mnc-sol-ls.py/g' submit.sh
       sed -i -e 's/mnc-sol.py/mnc-sol-ls.py/g' submit.sh
    elif [[ $spin == 'IS' ]]; then
       sed -i -e 's/mnc-sol-ls.py/mnc-sol-is.py/g' submit.sh
       sed -i -e 's/mnc-sol-is.py/mnc-sol-is.py/g' submit.sh
       sed -i -e 's/mnc-sol-hs.py/mnc-sol-is.py/g' submit.sh
       sed -i -e 's/mnc-sol-ns.py/mnc-sol-is.py/g' submit.sh
       sed -i -e 's/mnc-sol.py/mnc-sol-is.py/g' submit.sh
    elif [[ $spin == 'HS' ]]; then
       sed -i -e 's/mnc-sol-ls.py/mnc-sol-hs.py/g' submit.sh
       sed -i -e 's/mnc-sol-is.py/mnc-sol-hs.py/g' submit.sh
       sed -i -e 's/mnc-sol-hs.py/mnc-sol-hs.py/g' submit.sh
       sed -i -e 's/mnc-sol-ns.py/mnc-sol-hs.py/g' submit.sh
       sed -i -e 's/mnc-sol.py/mnc-sol-hs.py/g' submit.sh
    fi
done
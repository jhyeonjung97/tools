if [[ $1 -lt 10 ]]; then
    time='0'$1
else
    time=$1
fi
sed -i "/#SBATCH -t/c\#SBATCH -t $time:00:00" submit.sh
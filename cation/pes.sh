for i in {0..29}; do
    j=$(printf "%02d" $i)
    mkdir ${j}_
    mv restart_$j.json ${j}_/restart.json
done
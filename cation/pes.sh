for i in {1..30}; do
    j=$(printf "%02d" $i)
    mkdir -p ${j}_
    mv restart_$j.json ${j}_/restart.json
done
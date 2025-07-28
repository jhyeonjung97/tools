# basic
mkdir -p no3rr/
mv *.json no3rr/; mv no3rr/h*.json no3rr/oh*.json no3rr/o.json ./
python ~/bin/tools/pourbaix/pourbaix.py \
--figx 4 --figy 4 --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 \
--cmap-2d Greens --cmin-2d 0.1 --cmax-2d 0.6 --cgap-2d 0.0 --suffix HER

mv *.json no3rr/; mv no3rr/h*.json no3rr/oh*.json no3rr/o.json no3rr/vac*.json ./
python ~/bin/tools/pourbaix/pourbaix.py \
--figx 4 --figy 4 --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 \
--cmap-2d Greens --cmin-2d 0.3 --cmax-2d 0.6 --cgap-2d 0.0 --suffix HER_vac

mv no3rr/* ./; mv h*.json oh*.json o.json *mono*.json no3rr/
python ~/bin/tools/pourbaix/pourbaix.py \
--figx 4 --figy 4 --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 \
--cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.6 --cgap-2d 0.0 --suffix NO3RR_bi

mv no3rr/* ./; mv h*.json oh*.json o.json *bi*.json no3rr/
python ~/bin/tools/pourbaix/pourbaix.py \
--figx 4 --figy 4 --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 \
--cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.6 --cgap-2d 0.0 --suffix NO3RR_mono

mv no3rr/* ./; mv nh.json nh2*.json nh3*.json no3rr/
python ~/bin/tools/pourbaix/pourbaix.py --hybrid \
--thermo-data ./thermodynamic_data.jsonc \
--figx 4 --figy 4 --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 \
--cmap Greens --cmin 0.1 --cmax 0.9 --cgap 0.0 \
--cmap-2d Blues --cmin-2d 0.6 --cmax-2d 0.9 --cgap-2d 0.0 --suffix NO3RR_oxi

mv no3rr/* ./
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --png \
--figx 4 --figy 4 --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 --suffix NO3RR

python ~/bin/tools/pourbaix/pourbaix.py --hybrid --no-bulk \
--figx 4 --figy 4 --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 \
--cmap Greens --cmin 0.1 --cmax 0.9 --cgap 0.0 \
--cmap-2d Blues --cmin-2d 0.6 --cmax-2d 0.9 --cgap-2d 0.0 --suffix NO3RR

rm -r no3rr/
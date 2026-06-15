# surface (OxHy adsorbates)
mkdir -p P_pois_sqrt3_sqrt7_v2/
 cp *.json P_pois_sqrt3_sqrt7_v2/
 cp label.csv P_pois_sqrt3_sqrt7_v2/
 cd P_pois_sqrt3_sqrt7_v2/
# mv no3rr/h*.json no3rr/oh*.json no3rr/o.json .
# python ../../pourbaix.py --legend-out --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 \
# --colors-2d '#5c90c2' '#f0f9db' '#fec779' --suffix OxHy

# surface (NOxHy adsorbates)
# mv no3rr/* .
# mv h*.json oh*.json o.json no3rr/
python ../../pourbaix.py --legend-out --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 --png
# --cmap-2d PiYG_r --cmin-2d 0.15 --cmax-2d 0.85 --cgap-2d 0.0 --suffix NOxHy

# hybrid (bulk)
# mv P_pois/* .
python ../../pourbaix.py --legend-out --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 --hybrid --png --OER --HER
# --colors-bulk white darkgray gainsboro white darkgray

# hybrid
python ../../pourbaix.py --legend-out --Umin -2 --Umax 2 --pHmin -2 --pHmax 16 --hybrid --no-bulk --OER --HER \
 --cmap RdYlBu_r --cmin 0.15 --cmax 0.85 --cgap 0.0 --colors-2d '#D14895' '#821616'

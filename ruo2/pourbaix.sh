# basic
rm OH_OO.json
python ~/bin/HybridPB/pourbaix.py --hybrid --legend-out --show-transitions \
--Umin -0.5 --Umax 2.5 --Gmin -15 --Gmax 15 --cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0

python ~/bin/HybridPB/pourbaix.py --hybrid --legend-out --no-bulk \
--Umax 2 --Gmin -15 --Gmax 15 --cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 --suffix RdYlBu
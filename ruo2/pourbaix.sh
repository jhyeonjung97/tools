# basic
rm OH_OO.json

python ~/bin/HybridPB/pourbaix.py --legend-in --show-transitions \
--Umin 0.0 --Umax 2.0 --Gmin -5 --Gmax 5 --cmin-1d 0.1 --cmax-1d 0.9 --cgap-1d 0.2

python ~/bin/HybridPB/pourbaix.py --hybrid --show-transitions --no-bulk \
--Umin 0.0 --Umax 2.0 --Gmin -15 --Gmax 15 --cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0

python ~/bin/HybridPB/pourbaix.py --hybrid --no-bulk \
--Umax 2.5 --Gmin -15 --Gmax 15 --cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 --suffix RdYlBu

python ~/bin/HybridPB/pourbaix.py --hybrid --legend-out --no-bulk \
--Umax 2.5 --Gmin -15 --Gmax 15 --cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 --suffix RdYlBu
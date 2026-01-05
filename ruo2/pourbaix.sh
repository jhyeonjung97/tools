# basic
rm OH_OO.json

python ~/bin/HybridPB/pourbaix.py --legend-in --show-transitions \
--Umin 0.0 --Umax 2.0 --Gmin -5 --Gmax 5 --cmin-1d 0.1 --cmax-1d 0.9 --cgap-1d 0.2

python ~/bin/HybridPB/pourbaix.py --hybrid --show-transitions --no-bulk \
--Umin 0.0 --Umax 2.0 --Gmin -15 --Gmax 15 --cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0

python ~/bin/HybridPB/pourbaix.py --hybrid --no-bulk \
--Umin -0.5 --Umax 2.5 --Gmin -15 --Gmax 15 --cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 --suffix RdYlBu

python ~/bin/HybridPB/pourbaix.py --hybrid --no-bulk --legend-out \
--Umin -0.5 --Umax 2.5 --Gmin -15 --Gmax 15 --cmap-2d RdYlBu \
--cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 --colors-bulk white gainsboro whitesmoke white

python ~/bin/HybridPB/pourbaix.py --hybrid --no-bulk \
--Umin -0.5 --Umax 2.5 --Gmin -15 --Gmax 15 --cmap-2d RdYlBu \
--cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 --colors-bulk white whitesmoke white

python ~/bin/HybridPB/pourbaix.py --hybrid --Umin -2 --Umax 2 --cmin 0.1 --cmax 0.6

--cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 --cmin 0.0 --cmax 0.3

python ~/bin/HybridPB/pourbaix.py --hybrid --no-bulk --colors-bulk white whitesmoke white \
--Umin -0.5 --Umax 2.5 --cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 --legend-out --suffix top

python ~/bin/HybridPB/pourbaix.py --hybrid --no-bulk \
--Umin -0.5 --Umax 2.5 --cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 --legend-out --suffix vojvodic

python ~/bin/HybridPB/pourbaix.py \
--Umin -0.5 --Umax 2.5 --Gmin -15 --Gmax 15 \
--cmap-2d RdYlBu --cmin-2d 0.3 --cmax-2d 0.8 --cgap-2d 0.0 \
--cmin 0.0 --cmax 0.6
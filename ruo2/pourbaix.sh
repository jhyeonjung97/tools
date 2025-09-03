# basic
mv */*.json .
python ~/bin/HybridPB/pourbaix.py --hybrid --legend-out \
--Umax 2 --Gmin -15 --Gmax 15 --cmap-2d Purples --cmin-2d 0.2 --cmax-2d 0.8
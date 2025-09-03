# basic
mv */*.json .; mkdir -p K/; mv *.json K/; mv K/*KK-*.json .
python ~/bin/HybridPB/pourbaix.py --legend-out --hybrid \
--Umax 2 --Gmin -15 --Gmax 15 --cmap-2d Purples --cmin-2d 0.2 --cmax-2d 0.8
# basic
mv */*.json .; mv *OOH*json OOH/
python ~/bin/HybridPB/pourbaix.py --hybrid --legend-out \
--Umax 2 --cmap-2d Purples --cmin-2d 0.2 --cmax-2d 0.8 --Gmin -15 --Gmax 15 --suffix H2O

mv */*.json .; mv *OH2*json OH2/
python ~/bin/HybridPB/pourbaix.py --hybrid --legend-out \
--Umax 2 --cmap-2d Purples --cmin-2d 0.2 --cmax-2d 0.8 --Gmin -15 --Gmax 15 --suffix OOH

mv */*.json .; mv *OOH*json OOH/; mv *OH2*json OH2/
python ~/bin/HybridPB/pourbaix.py --hybrid --legend-out \
--Umax 2 --cmap-2d Purples --cmin-2d 0.2 --cmax-2d 0.8 --Gmin -15 --Gmax 15

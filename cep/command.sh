# basic
mv */*.json .; mkdir -p K/; mv *K-*.json K/
python ~/bin/HybridPB/pourbaix.py --legend-out --hybrid \
--suffix O --Umax 2 --cmin-2d 0.2 --cmax-2d 0.8 --cgap-2d 0.2

mv */*.json .; mkdir -p O/; mv O*.json O/
python ~/bin/HybridPB/pourbaix.py --legend-out --hybrid \
--suffix K --Umax 2 --cmin-2d 0.2 --cmax-2d 0.8 --cgap-2d 0.2

mv */*.json .; mkdir -p V/; mv *Mn*.json V/
python ~/bin/HybridPB/pourbaix.py --legend-out --hybrid \
--suffix V --Umax 2 --cmin-2d 0.2 --cmax-2d 0.8 --cgap-2d 0.2

mv */*.json .
python ~/bin/HybridPB/pourbaix.py --legend-out --hybrid \
--Umax 2 --cmin-2d 0.2 --cmax-2d 0.8 --cgap-2d 0.2
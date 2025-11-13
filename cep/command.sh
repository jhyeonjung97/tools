# basic
# mv */*.json .; mkdir -p K/; mv *K-*.json K/
# python ~/bin/HybridPB/pourbaix.py --legend-out --hybrid \
# --suffix O --Umax 2 --cmin-2d 0.2 --cmax-2d 0.8 --cgap-2d 0.2

# mv */*.json .; mkdir -p O/; mv O*.json O/
# python ~/bin/HybridPB/pourbaix.py --legend-out --hybrid \
# --suffix K --Umax 2 --cmin-2d 0.2 --cmax-2d 0.8 --cgap-2d 0.2

# mv */*.json .; mkdir -p V/; mv *Mn*.json V/
# python ~/bin/HybridPB/pourbaix.py --legend-out --hybrid \
# --suffix V --Umax 2 --cmin-2d 0.2 --cmax-2d 0.8 --cgap-2d 0.2

# mv */*.json .
# python ~/bin/HybridPB/pourbaix.py --legend-out --hybrid --no-bulk --show-ref \
# --Gmin -40 --Gmax 40 --Umin -2 --Umax 3 --cmin-2d 0.3 --cmax-2d 0.7 --cgap-2d 0.2

python ~/bin/HybridPB/pourbaix.py --legend-out --hybrid --no-bulk \
--Gmin -40 --Gmax 40 --Umin -2 --Umax 3 --cmap-2d RdYlBu --suffix RdYlBu38 --cmin-2d 0.3 --cmax-2d 0.8 

python ~/bin/HybridPB/pourbaix.py --legend-out --hybrid --no-bulk \
--Gmin -40 --Gmax 40 --Umin -2 --Umax 3 --cmap-2d RdYlBu --suffix RdYlBu39 --cmin-2d 0.3 --cmax-2d 0.9 


python ~/bin/HybridPB/pourbaix.py --legend-up --hybrid --no-bulk \
--Gmin -40 --Gmax 40 --Umin -2 --Umax 3 --cmap-2d RdYlBu --cmin-2d 0.2 --cmax-2d 0.9

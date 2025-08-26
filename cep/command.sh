# basic
mv */*.json .; mkdir -p K/; mv K*.json K/
python ~/bin/HybridPB/pourbaix.py \
--suffix O --Umax 2 --cmin-2d 0.1 --cmax-2d 0.45 --cgap-2d 0.0
python ~/bin/HybridPB/pourbaix.py --gc \
--suffix O --Umax 2 --cmin-2d 0.1 --cmax-2d 0.45 --cgap-2d 0.0
python ~/bin/HybridPB/pourbaix.py --legend-out \
--suffix O --Umax 2 --cmin-2d 0.1 --cmax-2d 0.45 --cgap-2d 0.0
python ~/bin/HybridPB/pourbaix.py --legend-out --gc \
--suffix O --Umax 2 --cmin-2d 0.1 --cmax-2d 0.45 --cgap-2d 0.0

mv */*.json .; mkdir -p O/; mv O*.json O/
python ~/bin/HybridPB/pourbaix.py \
--suffix K --Umax 2 --cmin-2d 0.55 --cmax-2d 0.9 --cgap-2d 0.0
python ~/bin/HybridPB/pourbaix.py --gc \
--suffix K --Umax 2 --cmin-2d 0.55 --cmax-2d 0.9 --cgap-2d 0.0
python ~/bin/HybridPB/pourbaix.py --legend-out \
--suffix K --Umax 2 --cmin-2d 0.55 --cmax-2d 0.9 --cgap-2d 0.0
python ~/bin/HybridPB/pourbaix.py --legend-out --gc \
--suffix K --Umax 2 --cmin-2d 0.55 --cmax-2d 0.9 --cgap-2d 0.0

mv */*.json .
python ~/bin/HybridPB/pourbaix.py --hybrid --legend-out \
--suffix K --Umax 2 --cmin-2d 0.2 --cmax-2d 0.8 --cgap-2d 0.2
python ~/bin/HybridPB/pourbaix.py --hybrid --legend-out --gc \
--suffix K --Umax 2 --cmin-2d 0.2 --cmax-2d 0.8 --cgap-2d 0.2

python ~/bin/HybridPB/pourbaix.py --hybrid --concentration 1e-99 --legend-out \
--suffix K99 --Umax 2 --cmin-2d 0.2 --cmax-2d 0.8 --cgap-2d 0.2
python ~/bin/HybridPB/pourbaix.py --hybrid --concentration 1e-99 --legend-out --gc \
--suffix K99 --Umax 2 --cmin-2d 0.2 --cmax-2d 0.8 --cgap-2d 0.2
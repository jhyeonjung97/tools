# basic
python ~/bin/tools/pourbaix/co-sac.py --figx 4 --figy 4 --cmap-2d Oranges \
--cmin-2d 0.1 --cmax-2d 0.6
python ~/bin/tools/pourbaix/co-sac.py --figx 4 --figy 4 --cmap-2d Oranges \
--cmin-2d 0.1 --cmax-2d 0.6 --legend-in --suffix legend

# gc
python ~/bin/tools/pourbaix/co-sac.py --figx 4 --figy 4 --cmap-2d Oranges \
--cmin-2d 0.1 --cmax-2d 0.6 --cmin-2d 0.1 --cmax-2d 0.9 --gc
python ~/bin/tools/pourbaix/co-sac.py --figx 4 --figy 4 --cmap-2d Oranges \
--cmin-2d 0.1 --cmax-2d 0.6 --cmin-2d 0.1 --cmax-2d 0.9 --gc --legend-in --suffix legend --Gmin -8 --Gmax 8

# hybrid
python ~/bin/tools/pourbaix/co-sac.py --gc \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/co-sac.py --gc --legend-in \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/co-sac.py --gc --legend-out \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8

# gc + hybrid
python ~/bin/tools/pourbaix/co-sac.py --gc --hybrid \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.3 --cmax-2d 0.5 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/co-sac.py --gc --hybrid --legend-in \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.3 --cmax-2d 0.5 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/co-sac.py --gc --hybrid --legend-out \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.3 --cmax-2d 0.5 --Gmin -8 --Gmax 8
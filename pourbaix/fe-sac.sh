# basic
python ~/bin/tools/pourbaix/pourbaix.py \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/pourbaix.py --legend-in \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/pourbaix.py --legend-out \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8

# gc
python ~/bin/tools/pourbaix/pourbaix.py --gc \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/pourbaix.py --gc --legend-in \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/pourbaix.py --gc --legend-out \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8

# hybrid
python ~/bin/tools/pourbaix/pourbaix.py --hybrid \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --legend-in \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --legend-out \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.2 --cmax-2d 0.9 --Gmin -8 --Gmax 8

# gc + hybrid
python ~/bin/tools/pourbaix/pourbaix.py --gc --hybrid \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.3 --cmax-2d 0.5 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/pourbaix.py --gc --hybrid --legend-in \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.3 --cmax-2d 0.5 --Gmin -8 --Gmax 8
python ~/bin/tools/pourbaix/pourbaix.py --gc --hybrid --legend-out \
--figx 4 --figy 4 --cmin 0.1 --cmax 0.7 --cmap-2d Oranges --cmin-2d 0.3 --cmax-2d 0.5 --Gmin -8 --Gmax 8
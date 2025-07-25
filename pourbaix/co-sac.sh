# basic
python ~/bin/tools/pourbaix/pourbaix.py \
--figx 4 --figy 4 --cmin 0.2 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5
python ~/bin/tools/pourbaix/pourbaix.py --legend-in \
--figx 4 --figy 4 --cmin 0.2 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5
python ~/bin/tools/pourbaix/pourbaix.py --legend-out \
--figx 4 --figy 4 --cmin 0.2 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5

# hybrid (bulk)
python ~/bin/tools/pourbaix/pourbaix.py --hybrid \
--figx 4 --figy 4 --cmax 0.1 --cmax 0.5 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --legend-in \
--figx 4 --figy 4 --cmax 0.1 --cmax 0.5 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --legend-out \
--figx 4 --figy 4 --cmax 0.1 --cmax 0.5 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5

# hybrid (hybrid)
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --no-bulk \
--figx 4 --figy 4 --cmin 0.2 --cmax 0.3 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --no-bulk --legend-in \
--figx 4 --figy 4 --cmin 0.2 --cmax 0.3 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5 --Gmin -8 --Gmax 8 --pH 0
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --no-bulk --legend-in \
--figx 4 --figy 4 --cmin 0.2 --cmax 0.3 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5 --Gmin -8 --Gmax 8 --pH 7
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --no-bulk --legend-in \
--figx 4 --figy 4 --cmin 0.2 --cmax 0.3 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5 --Gmin -8 --Gmax 8 --pH 14
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --no-bulk --legend-out \
--figx 4 --figy 4 --cmin 0.2 --cmax 0.3 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5 --Gmin -8 --Gmax 8 --pH 0
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --no-bulk --legend-out \
--figx 4 --figy 4 --cmin 0.2 --cmax 0.3 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5 --Gmin -8 --Gmax 8 --pH 7
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --no-bulk --legend-out \
--figx 4 --figy 4 --cmin 0.2 --cmax 0.3 --cmap-2d Blues --cmin-2d 0.1 --cmax-2d 0.5 --Gmin -8 --Gmax 8 --pH 14
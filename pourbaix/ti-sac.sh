# basic
python ~/bin/tools/pourbaix/pourbaix.py \
--figx 4 --figy 4 --Umin -2 --Umax 2 --cmin 0.2 --cmax 0.3 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4
python ~/bin/tools/pourbaix/pourbaix.py --legend-in \
--figx 4 --figy 4 --Umin -2 --Umax 2 --cmin 0.2 --cmax 0.3 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4
python ~/bin/tools/pourbaix/pourbaix.py --legend-out \
--figx 4 --figy 4 --Umin -2 --Umax 2 --cmin 0.2 --cmax 0.3 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4

# hybrid (bulk)
python ~/bin/tools/pourbaix/pourbaix.py --hybrid \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc6
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --legend-in \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc6
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --legend-out \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc6

# hybrid (hybrid)
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --no-bulk \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc6
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --no-bulk --legend-in \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc6
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --no-bulk --legend-out \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc6

# concentration (10^-4)
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --conc 0.0001 \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc4
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --conc 0.0001 --legend-in \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc4
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --conc 0.0001 --legend-out \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc4

# concentration (10^-8)
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --conc 0.00000001 \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc8
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --conc 0.00000001 --legend-in \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc8
python ~/bin/tools/pourbaix/pourbaix.py --hybrid --conc 0.00000001 --legend-out \
--figx 4 --figy 4 --pHmin -2 --pHmax 16 --Umin -2 --Umax 2 --cmin 0.0 --cmax 0.6 --cmap-2d Purples --cmin-2d 0.1 --cmax-2d 0.4 --suffix conc8
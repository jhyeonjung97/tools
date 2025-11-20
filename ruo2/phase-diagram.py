import os
import mendeleev as md

for i in list[int](range(21, 31)) + list(range(39, 49)) + [57] + list(range(72, 81)):
    element = md.element(i).symbol
    png_name = f'{i}_{element}_phase_diagram.png'
    os.system(f'python ~/bin/tools/ruo2/Re-phase-diagram.py {element} O -o {png_name} --no-show')
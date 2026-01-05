python ~/bin/HybridPB/pourbaix.py --legend-out --show-transitions

python ~/bin/tools/irfe/energetics.py -V 0.0 --xmin -0.5 --xmax 8.5 --ymin -0.5 --ymax 8.5 --no-show-text
python ~/bin/tools/irfe/energetics.py -V 1.0 --xmin -0.5 --xmax 7.5 --ymin -0.5 --ymax 2.5 --no-show-text
python ~/bin/tools/irfe/energetics.py -V 1.2 --xmin -0.5 --xmax 9.5 --ymin -1.5 --ymax 1.0 --no-show-text
python ~/bin/tools/irfe/energetics.py -V 1.2 --xmin -0.5 --xmax 9.5 --ymin -1.5 --ymax 1.5 --no-show-text
python ~/bin/tools/irfe/energetics.py -V 1.5 --xmin -0.5 --xmax 9.5 --ymin -2.5 --ymax 1.0 --no-show-text

python ~/bin/tools/irfe/energetics_dual.py -V 0.0 --no-show-text --show-legend \
--xmin -0.5 --xmax 6.0 --ymin -0.5 --ymax 9.5 --gap 0.5 --figsize 8 3

python ~/bin/tools/irfe/energetics_dual.py -V 1.23 --no-show-text \
--xmin -0.5 --xmax 6.5 --ymin -1.5 --ymax 2.5 --suffix 1.23 --gap 0.5 --figsize 10 4

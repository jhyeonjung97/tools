from ase import io
from ase.build import surface

bulk = io.read("restart.json")
slab = surface(bulk, (1, 1, 1), layers=4, vacuum=15.0)
slab = slab.repeat((3, 3, 1))
io.write("surface111.json", slab)
print(slab)

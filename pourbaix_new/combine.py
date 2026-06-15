from ase.db import connect
from ase.io import read
import glob
import os

if os.path.exists('combined.db'):
    print('combined.db already exists')
    exit()

# Create or open the database
db = connect('combined.db')

# Loop over all JSON files
for filename in sorted(glob.glob('*.json')):
    atoms = read(filename)
    # Write to database
    db.write(atoms, filename=filename)

print(f"All JSON files have been written to combined.db")
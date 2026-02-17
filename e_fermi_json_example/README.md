# Storing and Reading E-fermi in ASE JSON Files

How to store VASP Fermi energy (E-fermi) in ASE-format JSON and read it back with ASE.

## Summary

| Step | Method |
|------|--------|
| **Write** | `ase.db.connect(file).write(atoms, e_fermi=value)` or `db.update(id, e_fermi=value)` |
| **Read** | `read_json_with_info(file)` → `atoms.info['e_fermi']` |

Using only `ase.io.read('file.json')` does **not** load custom keys like `e_fermi`. You need to write with **ase.db** and read with **read_json_with_info()** as in this example.

---

## 1. Writing e_fermi to a JSON File

### Option A: Add e_fermi to an existing JSON file

When you have already parsed E-fermi from OUTCAR and want to add it to an existing JSON:

```python
import ase.db

json_path = 'my_structure.json'
e_fermi = 19.99   # value parsed from OUTCAR, etc.

with ase.db.connect(json_path, serial=True) as db:
    db.update(1, e_fermi=e_fermi)   # add e_fermi to row id=1
```

Use the appropriate row `id` if your database has multiple rows.

### Option B: Write Atoms to JSON with e_fermi from the start

When creating the JSON file:

```python
from ase.io import read
import ase.db

# e.g. read from traj or structure file
atoms = read('structure.traj')
e_fermi = 19.99

# ase.db allows storing key-value pairs together with the structure
with ase.db.connect('output.json', serial=True) as db:
    db.write(atoms, e_fermi=e_fermi)
```

### Option C: Parse E-fermi from OUTCAR and write to JSON

```python
import ase.db

def get_e_fermi_from_occar(outcar_path='OUTCAR'):
    with open(outcar_path) as f:
        text = f.read()
    return float(text.split('E-fermi :')[-1].split()[0])

json_path = 'final_with_calculator.json'
e_fermi = get_e_fermi_from_occar('OUTCAR')

with ase.db.connect(json_path, serial=True) as db:
    db.update(1, e_fermi=e_fermi)
```

---

## 2. Reading with ASE

### Option 1: Use read_json_with_info() (recommended)

Using `read_json_with_info.py` in this folder gives direct access via `atoms.info['e_fermi']`.

```python
from read_json_with_info import read_json_with_info

atoms = read_json_with_info('my_structure.json')
print(atoms.info['e_fermi'])   # 19.99
```

### Option 2: Use ase.db only

```python
from ase.db import connect

with connect('my_structure.json') as db:
    atoms = db.get_atoms(1, add_additional_information=True)

e_fermi = atoms.info['key_value_pairs']['e_fermi']
```

---

## Running the Example Scripts

1. **Create a test sample.json and add e_fermi**
   ```bash
   python 01_add_e_fermi_to_json.py
   ```
   → Creates `sample.json` with e_fermi=5.123 stored.

2. **Read the JSON and print e_fermi**
   ```bash
   python 02_read_json_with_e_fermi.py
   ```
   → Prints the value of `atoms.info['e_fermi']`.

3. **(Optional) Parse E-fermi from OUTCAR and add to an existing JSON**
   ```bash
   python 03_add_e_fermi_from_occar.py --json my.json --occar OUTCAR --id 1
   ```

Adjust `json_path` or script arguments to match your own file paths as needed.

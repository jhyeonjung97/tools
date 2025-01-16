import re

def parse_zval_from_potcar(potcar_content):
    symbols = []
    zvals = []
    symbol_pattern = re.compile(r"TITEL\s*=\s*PAW_PBE\s*(\w+)")
    zval_pattern = re.compile(r"ZVAL\s*=\s*(\d+\.\d+)")

    lines = potcar_content.splitlines()
    for line in lines:
        match1 = symbol_pattern.search(line)
        if match1:
            symbols.append(match1.group(1))
        match2 = zval_pattern.search(line)
        if match2:
            zvals.append(int(float(match2.group(1))))

    if len(symbols) != len(zvals):
        raise ValueError("Mismatch between number of symbols and ZVALs in POTCAR content.")

    zval_mapping = dict(zip(symbols, zvals))
    return zval_mapping


with open("POTCAR", "r") as potcar_file:
    potcar_content = potcar_file.read()

zval_mapping = parse_zval_from_potcar(potcar_content)
print(zval_mapping)
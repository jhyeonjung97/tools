import re

lobsterout_file = "lobsterout"

# Read the lobsterout file
with open(lobsterout_file, "r") as f:
    lines = f.readlines()

# Extract spilling values using regex
spilling_values = []
for line in lines:
    match = re.search(r"abs\. charge spilling:\s+([\d.]+)%", line)
    if match:
        spilling_values.append(float(match.group(1)))

# Check if any spilling is greater than 5%
print(spilling_values)
if any(spill > 5.0 for spill in spilling_values):
    warning_msg = "\nWARNING: Charge spilling exceeds 5%! Check your basis set and projection parameters.\n"

    # Append warning to lobsterout
    with open(lobsterout_file, "a") as f:
        f.write(warning_msg)
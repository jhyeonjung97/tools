#! /usr/bin/envpython 

### Created by: Md. Delowar Hossain ###

### Modified by: Evan Carlson ###
### Date: 08/25/2024 ###

### This script runs GCDFT calculations ###

### Imports ###
import os
import numpy as np
from ase.io import *
from ase.visualize import *
import matplotlib.pyplot as plt

# Read the OUTCAR file
outcar_path = os.path.join(os.getcwd(), 'OUTCAR')
vaspout_path = os.path.join(os.getcwd(), 'vasp.out')

# Extract the # of electrons
with open(outcar_path, 'r') as f:
    outcar = f.read()
neutral_electrons = float(outcar.split('NELECT =')[1].split()[0]) # Number of electrons in the uncharged slab

# Extract the Fermi energy
with open(vaspout_path, 'r') as f:
    vaspout = f.read()
fermi_shift = float(vaspout.split('FERMI_SHIFT =')[-1].split()[0]) # Fermi shift from VaspSol (eV)

## Constants
WF_SHE = 4.4 # Experimental WF of SHE (eV)
e_chg = 1.602e-19 # Elementary Charge (Coulombs)

## Are all charge states used for the fitting?
FULL = False
# FULL = True

## Charge states to analyze
num_electrons = [+1.0, +0.5, 0.0, -0.5, -1.0, -1.5, -2.0]
# num_electrons = [+0.5, 0.0, -0.5, -1.0, -1.5]

folders =[]

for numElectrons in num_electrons:
    folders+=["nelect_"+str(neutral_electrons+numElectrons)]

# print(folders)

## Empty lists to store the data
nelects= []   # Number of electrons in the slab
charges = []  # Charge of the slab 
energies = [] # Total energy of the slab
fermi_energies = [] # Fermi energy of the slab
gc_free_energies =[] # GC free energy of the slab
applied_potentials = [] # Potential (vs. SHE) of the slab
latticeParams = [] # Lattice parameters of the slab

## Loop over the directories and extract the data
for folder in folders:
   
    # Read the OUTCAR file
    outcar_path = os.path.join(folder, 'OUTCAR')
    with open(outcar_path, 'r') as f:
        outcar = f.read()

    # Extract the # of electrons
    nelect = float(outcar.split('NELECT =')[1].split()[0])
    # print(nelect)
    nelects.append(nelect)

    # Extract the charge
    charge = nelect-neutral_electrons
    # print(charge)
    charges.append(charge)

    # Extract the Fermi energy
    fermi_energy = float(outcar.split('E-fermi :')[-1].split()[0])
    # print(fermi_energy)
    fermi_energies.append(fermi_energy)

    # Calculate the applied potential (vs. SHE)
    applied_potential = (fermi_shift-fermi_energy)-WF_SHE 
    # print(applied_potential)
    applied_potentials.append(applied_potential)

    # Extract the DFT energy
    energy = float(outcar.split('energy  without entropy=')[-1].split()[0])
    # print(energy)
    energies.append(energy)

    # Calculate the GC free energy
    gc_free_energies.append(energy-charge*(fermi_energy)) 

    # Extract the lattice parameters
    aLat = float(outcar.split(' length of vectors')[-1].split()[0])
    bLat = float(outcar.split(' length of vectors')[-1].split()[1])
    latticeParams.append(aLat*bLat)
    # print(aLat)
    # print(bLat)

## Fit the data quadratically
p1, risid1, rank1, sing1, rcon1 = np.polyfit(charges, energies, 2, full=True)
p2, risid2, rank2, sing2, rcon2 = np.polyfit(applied_potentials, gc_free_energies, 2, full=True)
p3, risid3, rank3, sing3, rcon3 = np.polyfit(applied_potentials, charges, 2, full=True)
a1 = p1[0]; b1 = p1[1]; c1 = p1[2]
a2 = p2[0]; b2 = p2[1]; c2 = p2[2]
a3 = p3[0]; b3 = p3[1]; c3 = p3[2]

## Calculate the variance of the fits
ybar1 = np.sum(energies)/len(energies)
ybar2 = np.sum(gc_free_energies)/len(gc_free_energies)
ybar3 = np.sum(charges)/len(charges)
sstot1 = np.sum((energies-ybar1)**2)
sstot2 = np.sum((gc_free_energies-ybar2)**2)
sstot3 = np.sum((charges-ybar3)**2)

## Create a range of x values
x = np.linspace(min(applied_potentials), max(applied_potentials), 100)
x_charge = np.linspace(min(charges)-0.3, max(charges)+0.3, 100)

## Calculate the corresponding y values for each fit
y1_fit_EnVChg = a1*x_charge**2 + b1*x_charge + c1
y2_fit_GCEnVPot = a2*x**2 + b2*x + c2
y3_fit_ChgVPot = a3*x**2 + b3*x + c3
y4_fit_CapVPot = (2*a3*x + b3)*-1*e_chg/(aLat*bLat*(10**-8)**2)*10**6

## Plot DFT energy vs. charge
plt.figure(figsize=(8,5))
plt.plot(charges, energies, 'ro', label='data')
plt.plot(x_charge, y1_fit_EnVChg, label='quadratic fitting')
plt.xlabel("Number of Extra Electrons")
plt.ylabel("DFT_energy, eV")
plt.legend()
plt.text(0.1, 0.1, r'RSS=%.2f, ' % (risid1.item()) + 'R$^{\sf 2}$ = %.2f' % (1-risid1.item()/sstot1.item()), transform=plt.gca().transAxes) # Fit statistics

## Save and display the DFT energy vs. charge plot
if FULL:
    plt.savefig('energy4p4FULL.png')
else:
    plt.savefig('energy4p4.png')
plt.show()
plt.close()

## Plot GC free energy vs. applied potential (vs. SHE)
plt.figure(figsize=(8,5))
plt.plot(applied_potentials, gc_free_energies, 'go', label = 'data')
plt.plot(x, y2_fit_GCEnVPot, label='quadratic fitting')
plt.xlabel("Potential vs. SHE")
plt.ylabel("Free Energy, eV")
plt.legend()
plt.text(0.1, 0.1, r'RSS=%.2f, ' % (risid2.item()) + 'R$^{\sf 2}$ = %.2f' % (1-risid2.item()/sstot2.item()), transform=plt.gca().transAxes) # Fit statistics

## Save and display the GC free energy vs. applied potential plot
if FULL:
    plt.savefig('Free_energy4p4FULL.png')
else:
    plt.savefig('Free_energy4p4.png')
plt.show()
plt.close()

## Plot charge vs. applied potential (vs. SHE)
plt.figure(figsize=(8,5))
plt.plot(applied_potentials, charges, 'ro', label = 'data')
plt.plot(x, y3_fit_ChgVPot, label = 'quadratic fitting')
plt.legend()
plt.xlabel("Potential vs. SHE")
plt.ylabel("Number of Extra Electrons")
plt.text(0.1, 0.1, r'RSS=%.2f, ' % (risid3.item()) + 'R$^{\sf 2}$ = %.2f' % (1-risid3.item()/sstot3.item()), transform=plt.gca().transAxes) # Fit statistics

## Save and display the charge vs. applied potential plot
if FULL:
    plt.savefig('Charge4p4FULL.png')
else:
    plt.savefig('Charge4p4.png')
plt.show()
plt.close()

## Plot Fermi energy vs. applied potential (vs. SHE) 
plt.figure(figsize=(8,5))
plt.plot(applied_potentials,fermi_energies, '-o')
plt.xlabel("Potential vs. SHE")
plt.ylabel("Fermi_Energy, eV")

## Save and display the Fermi energy vs. applied potential plot
if FULL:
    plt.savefig('Fermi_energy4p4FULL.png')
else:
    plt.savefig('Fermi_energy4p4.png')
plt.show()
plt.close()

## Plot capacitance vs. applied potential (vs. SHE)
plt.figure(figsize=(8,5))
plt.plot(x, y4_fit_CapVPot, label='Fit')
plt.xlabel("Potential vs. SHE")
plt.ylabel(r"Capacitance, uF/cm$^2$")
plt.legend()

## Save and Display the Capacitance vs. Applied Potential Plot
if FULL:
    plt.savefig('Capacitance4p4FULL.png')
else:
    plt.savefig('Capacitance4p4.png')
plt.show()
plt.close()

## Save the GC Free Energy Fit and Data Points 
coeffs = 'Ax^2 + Bx + C, where A,B,C = ' + str(a2) + ', ' + str(b2) + ', ' + str(c2) # Fit coefficients of GC Free Energy vs. Applied Potential

if FULL:
    np.savetxt('Applied_Pot_vs_GCFE_Fit_FULL.dat', np.column_stack((x, y2_fit_GCEnVPot)), header=coeffs)
    np.savetxt('GCFE_data_FULL.dat', np.column_stack((applied_potentials, gc_free_energies)))
else:
    np.savetxt('Applied_Pot_vs_GCFE_Fit.dat', np.column_stack((x, y2_fit_GCEnVPot)), header=coeffs)
    np.savetxt('GCFE_data.dat', np.column_stack((applied_potentials, gc_free_energies)))

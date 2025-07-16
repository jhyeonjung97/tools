import os
import json
from collections import defaultdict
from mendeleev import element
from scipy import optimize
from scipy.integrate import cumtrapz
from ase.io import read
from ase.data import covalent_radii as CR
import numpy as np

try:
    with open('./ads_data.json', 'r') as f:
        ads_data = json.load(f)

    with open('./metal_fits.json', 'r') as f:
        cohp_fit_dict = json.load(f)
except:
    print('couldn\'t load data')

def read_COHPCAR(filename):
    with open(filename, 'r') as f:
        txt = f.readlines()

    #print(filename)
    ln_start = 3
    names = []

    for line in txt[2:]:
        if 'No' in line or 'Average' in line:
            ln_start += 1
            cut = line.split('(')[0]
            cut = cut.strip()
            names.extend([cut, cut + ' Integral'])
        else:
            break
    names *= 2

    txt = txt[ln_start:]
    engs = []
    proj_dict = defaultdict(list)

    for line in txt:
        ln = [float(a) for a in line.split()]
        #print(len(ln))

        engs.append(ln[0])

        end = int((len(ln)-1)/2)
        for i, proj in enumerate(ln[1:end+1]):
            proj_dict[names[i]].append(proj)
        start_spin_2 = end+1
        end_spin_2 = len(ln)
        for i, proj in enumerate(ln[start_spin_2:end_spin_2]):
            proj_dict[names[i]+' spin 2'].append(proj)
    return proj_dict, engs

def post_process_projections(proj_dict, engs):
    t2g_keys = []
    eg_keys = []
    p_t2g_keys = []
    p_eg_keys = []
    d_keys = []
    s_keys = []
    O_s_keys = []
    f_keys = []
    met_p = []
    subset_sums=['t2g', 'eg', 'd tot','neg d tot','pos d tot' ]

    for i in proj_dict:
        proj_dict[i] = np.array(proj_dict[i])
        if '2p' in i and 'Integral' not in i:
            if 'd_xy' in i or 'd_yz' in i or 'd_xz' in i:
                t2g_keys.append(i)
            elif 'd_z^2' in i or 'd_x^2' in i:
                eg_keys.append(i)

    for key in subset_sums: 
        proj_dict[key] = np.array([0.]* len(engs))
        proj_dict[key + ' Integral'] = np.array([0.]* len(engs))

    for key in eg_keys + t2g_keys:
        proj_dict['d tot'] += np.array(proj_dict[key]).copy()
        # filter to have just bonding/anit-bonding
        proj_dict['pos d tot'] += np.where(np.array(proj_dict[key]).copy()<0, 0, np.array(proj_dict[key]).copy())
        proj_dict['neg d tot'] += np.where(-1 * np.array(proj_dict[key]).copy()<0, 0, np.array(proj_dict[key]).copy())
        if 'spin' in key:
            splt = key.split()
            n_key = ' '.join([splt[0], 'Integral'] + splt[1:])
        else:
            n_key = key + ' Integral'
        proj_dict['d tot Integral'] +=np.array(proj_dict[n_key]).copy()
        # integrate bonding/anti-bonding
        proj_dict['pos d tot Integral'] += cumtrapz(np.where(np.array(proj_dict[key]).copy()<0, 0, np.array(proj_dict[n_key]).copy()),
                                                    engs, initial=0)
        proj_dict['neg d tot Integral'] += cumtrapz(np.where(-1 * np.array(proj_dict[key]).copy()<0, 0, np.array(proj_dict[n_key]).copy()), 
                                                    engs, initial=0)

        #else:
        #    proj_dict['d tot Integral'] +=np.array(proj_dict[key + ' Integral']).copy()
        #    proj_dict['pos d tot Integral'] += np.where(np.array(proj_dict[key + ' Integral']).copy()<0,
        #                                                0, np.array(proj_dict[key]).copy())
        #    proj_dict['neg d tot Integral'] += np.where(-1 * np.array(proj_dict[key + ' Integral']).copy()<0,
        #                                                 0, np.array(proj_dict[key]).copy())
        #proj_dict['neg d tot Integral'] = simps(proj_dict['neg d tot'], x=engs)
        #proj_dict['pos d tot Integral'] = simps(proj_dict['pos d tot'], x=engs)
    #proj_dict['pos d tot Integral'] += simps()
    #proj_dict['neg d tot Integral'] += np.where(-1 * np.array(proj_dict[key]).copy()<0, 0, np.array(proj_dict[d_]).copy())

    return proj_dict, engs


def parse_list_file(filename, return_avg_bond_len=False, orbital_interaction='d', orbital_removed='s'):
    with open(filename, 'r') as f:
        lines = f.readlines()

    total_integral = 0
    per_bond_integral_dict = defaultdict(float)
    bond_lengths = []
    bl_icohp_pairs = []
    valence = None
    for line in lines:
        if 'spin' in line:
            continue
        splt = line.split()
        int_num = int(splt[0])
        orbitals = [splt[1], splt[2]]
        full_config = [[a.split('_')[0].rstrip('0123456789')] + a.split('_')[1:] for a in orbitals]
        orbitals = [a.split('_')[1:] for a in orbitals]
        if orbitals == [[],[]]: # this indicates a 'total cohp' line
            #grabbing out the bond length
            bond_lengths.append(float(splt[3]))
            if valence is None:
                sym = [a for a in [splt[1].rstrip('0123456789'), splt[2].rstrip('0123456789')] if a != 'O'][0]
                valence = element(sym).ec.get_valence().to_str()
                valence = [a[:2] for a in valence.split()]
                if sym == 'Pd':
                    valence += ['5s']  # Annoying...
            continue
        # sum isn't working for some reason uggghhh
        orbitals = orbitals[0] + orbitals[1]
        orbs = ''
        for o in orbitals:
            orbs += o
        orbitals = orbs
        #print(orbs)
        if orbital_interaction == 'valence':
            not_screened = [a in orbitals for a in valence]
        elif orbital_interaction == 'valence_p':
            p_valence = [a for a in valence if 'p' in a]
            not_screened = [a in orbitals for a in p_valence]
        elif orbital_interaction == 'valence_p_s':
            s_p_valence = [a for a in valence if 'd' not in a]
            not_screened = [a in orbitals for a in s_p_valence]
            #print(orbitals, s_p_valence, valence, not_screened)
        else:
            if type(orbital_interaction) == str:
                not_screened = [orbital_interaction in orbitals and orbital_removed not in orbitals]
            else:
                is_screened = []
                for interaction in orbital_interaction:
                    is_not_screened.append(interaction in orbitals and orbital_removed not in orbitals)
                if True in is_not_screened:
                    not_screened = [True]

        #if True:
        if any(not_screened) == True:
            integral = float(line.split()[-1])
            per_bond_integral_dict[line.split()[0]] += integral
        #else:
        #    print(orbitals, valence)
    for key, item in per_bond_integral_dict.items():
        total_integral += item
        bl_icohp_pairs.append((bond_lengths[int(key)-1], item))
    
    avg_bond_len = np.mean(bond_lengths)
    if return_avg_bond_len:
        return total_integral, bl_icohp_pairs
    return total_integral

def parse_long_list_file(filename, return_avg_bond_len=False, orbital_interaction='d', orbital_removed='s'):
    with open(filename, 'r') as f:
        lines = f.readlines()

    total_integral = 0
    per_bond_integral_dict = defaultdict(float)
    fid = defaultdict(dict)

    bond_lengths = []
    bl_icohp_pairs = []
    valence = None
    for line in lines:
        if 'spin' in line:
            continue
        splt = line.split()
        int_num = int(splt[0])
        orbitals = [splt[1], splt[2]]
        full_config = [[a.split('_')[0].rstrip('0123456789')] + a.split('_')[1:] for a in orbitals]
        orbitals = [a.split('_')[1:] for a in orbitals]
        if orbitals == [[],[]]: # this indicates a 'total cohp' line
            #grabbing out the bond length
            bond_lengths.append(float(splt[3]))
            inds =  [int(''.join([i for i in ind if i.isdigit()]))-1 for ind in [splt[1],splt[2]]]
            sym = [a for a in [splt[1].rstrip('0123456789'), splt[2].rstrip('0123456789')]]
            fid[int_num]['indicies'] = inds
            fid[int_num]['bond_length'] = float(splt[3])
            if fid[int_num].get('integral') is None:
                fid[int_num]['integral'] = 0
            fid[int_num]['symbols'] = sym
            #if valence is None:
            m_sym = [a for a in sym if a != 'O'][0]
            valence = element(m_sym).ec.get_valence().to_str()
            valence = [a[:2] for a in valence.split()]
            if m_sym == 'Pd':
                valence += ['5s']  # Annoying...
            continue
        # sum isn't working for some reason uggghhh
        orbitals = orbitals[0] + orbitals[1]
        orbs = ''
        for o in orbitals:
            orbs += o
        orbitals = orbs
        #print(orbs)
        if orbital_interaction == 'valence':
            not_screened = [a in orbitals for a in valence]
        elif orbital_interaction == 'valence_p':
            p_valence = [a for a in valence if 'p' in a]
            not_screened = [a in orbitals for a in p_valence]
        elif orbital_interaction == 'valence_p_s':
            s_p_valence = [a for a in valence if 'd' not in a]
            not_screened = [a in orbitals for a in s_p_valence]
            #print(orbitals, s_p_valence, valence, not_screened)
        else:
            if type(orbital_interaction) == str:
                not_screened = [orbital_interaction in orbitals and orbital_removed not in orbitals]
            else:
                is_screened = []
                for interaction in orbital_interaction:
                    is_not_screened.append(interaction in orbitals and orbital_removed not in orbitals)
                if True in is_not_screened:
                    not_screened = [True]

        #if True:
        if any(not_screened) == True:
            integral = float(line.split()[-1])
            per_bond_integral_dict[line.split()[0]] += integral
            fid[int_num]['integral'] += integral
        #else:
        #    print(orbitals, valence)
    for key, item in per_bond_integral_dict.items():
        total_integral += item
        bl_icohp_pairs.append((bond_lengths[int(key)-1], item))

    avg_bond_len = np.mean(bond_lengths)
    return dict(fid)
    #if return_avg_bond_len:
    #    return total_integral, bl_icohp_pairs
    #return total_integral


def get_center_sums(ints_dict):
    # find all the metal centers
    metal_indicies = []
    completed_indicies = []
    centers_dict = {}
    for i, inter in ints_dict.items():
        #indicies = [(a,b) for a, b in zip(inter['indicies'], inter['symbols'])]
        #print({a:{'sum':0., 'symbol':b} for a,b in zip(inter['indicies'], inter['symbols'])})
        centers_dict.update({a:{'sum':0., 'symbol':b} for a,b in zip(inter['indicies'], inter['symbols'])})
        
    #metal_indicies = list(set(metal_indicies)) # removed repeated indicies
    #centers_dict = {a:None for a in metal_indicies}
    for i, inter in ints_dict.items():
        #if inter['indicies'] in completed_indicies or inter['indicies'][::-1] in completed_indicies:
        #    continue
        for index in inter['indicies']:
            centers_dict[index]['sum'] += inter['integral']
            centers_dict[index]['bond_length'] = inter['bond_length']
        completed_indicies.append(inter['indicies'])
    return centers_dict


def plot_integrals(int_file):
    integral_dict={}
    for row in [d3, d4, d5]:
        for i, metal in enumerate(row):
            if i < 10:
                dire = '0{}_{}O2'.format(i, metal)
            else:
                dire = '{}_{}O2'.format(i, metal)
            #os.chdir(dire+'/no_sym/')
            if not os.path.isfile('./'+dire+'/no_sym/'+int_file):
                print(dire+'/no_sym/'+int_file)
                continue
            total_int = parse_list_file(dire+'/no_sym/'+int_file)
            integral_dict[metal] = total_int

    for name, row in zip(['3d', '4d', '5d'], [d3, d4, d5]):
        plt.scatter([i for i,a in enumerate(row)], [integral_dict.get(a) for a in row], label=name)
        for ind, icob, txt in zip([i for i,a in enumerate(row)], [integral_dict.get(a) for a in row], row):
            if icob is not None:
                plt.text(ind, icob, txt)
    plt.legend()
    #plt.show()
    plt.close()
    return integral_dict


def get_bond_length(dire_):
    with open(dire_ + '/lobsterin', 'r') as f:
        txt = f.read()
    txt = [a.strip() for a in txt.split('atom')]
    ind1 = int(txt[-1].split()[0]) - 1 # not zero indexing for some reason
    ind2 = int(txt[-2].split()[0]) - 1
    atoms = read(dire_ + '/final_with_calculator.json')
    dist = atoms.get_distance(ind1, ind2, mic=True)
    return dist

def get_ionic_radius(element_name):
    # this function sucks, don't use it
    # but maybe it's a useful reference for making
    # your own
    #hs_elements = ['V', 'Cr', 'Mn', 'Fe', 'Co']
    hs_elements = []
    ionic_radii = element(element_name).ionic_radii
    filtered_list = []
    for radius in ionic_radii:
        if radius.coordination == 'VI':
            if element_name in hs_elements:
                if radius.spin == 'HS':
                    filtered_list.append(radius)
            else:
                if radius.spin != 'HS':
                    filtered_list.append(radius)
    #if element in hs_elements and filtered_list == []:

    relible_filtered = [a for a in filtered_list if a.most_reliable== True]
    # if one of them is the "most relible" return that, otherwise just pick the first one
    # converting from pm to angstrom
    if relible_filtered != []:
        return relible_filtered[0].ionic_radius / 100
    else:
        return filtered_list[0].ionic_radius / 100
    #ionic_radius = [a for a in element(element_name).ionic_radii if a.most_reliable == True][0].ionic_radius

def get_normal_oxidation_states(element_name):
    oxi_states = element(element_name).oxistates
    return oxi_states

def get_adsorption_energy(metal, oxidation_state, adsorbate, facet=None):
    os_facet_match = {2:'100',3:'012',4:'110',5:'18'}
    if facet is None:
        facet = os_facet_match[oxidation_state]
    try:
        #if oxidation_state != 5:
        #    if adsorbate == 'O-OH':
        #        return ads_data[metal][str(oxidation_state)][facet]['O'] - ads_data[metal][str(oxidation_state)][facet]['OH']
        #    else:
        #        return ads_data[metal][str(oxidation_state)][facet][adsorbate]
        #else:
        #    if adsorbate == 'O-OH':
        #        return ads_data[metal][str(oxidation_state)][facet]['O'] - ads_data[metal][str(oxidation_state)][facet]['OH']
        #    else:
        #        return ads_data[metal][str(oxidation_state)][facet][adsorbate]
        if adsorbate == 'O-OH':
            return ads_data[metal][str(oxidation_state)][facet]['O'] - ads_data[metal][str(oxidation_state)][facet]['OH']
        else:
            return ads_data[metal][str(oxidation_state)][facet][adsorbate]
    except:
        return None

def oxidation_state_subtraction(oxys, COHPs, target_ox, diff=2):
    sam_ind = oxys.index(target_ox)
    transformed = COHPs[index+2] - COHPs[index]
    return transformed

def predict_cohp(metal, bond_length):
    if metal not in cohp_fit_dict.keys():
        print(metal)
        return None
    return cohp_fit_dict[metal]['slope'] / bond_length ** 6 + cohp_fit_dict[metal]['intercept']

def predict_cohp_log(metal, bond_length):
    if metal not in cohp_fit_dict.keys():
        print(metal)
        return None
    return np.exp(cohp_fit_dict[metal]['log_intercept']) * bond_length ** cohp_fit_dict[metal]['log_slope']

def predict_cohp_kr(metal, bond_length):
    if metal not in cohp_fit_dict.keys():
        print(metal)
        return None
    a = cohp_fit_dict[metal]['kr_a']
    b = cohp_fit_dict[metal]['kr_b']
    x = bond_length
    return a/(x-b)**6

def predict_cohp_pauling(metal, bond_length):
    D1 = cohp_fit_dict[metal]['pf_D1']
    b = cohp_fit_dict[metal]['pf_b']
    return -1 * np.exp((D1-bond_length)/b)

def predict_bond_order_pauling(metal, bond_length):
    D1 = cohp_fit_dict[metal]['bo_D1']
    b = cohp_fit_dict[metal]['bo_b']
    return np.exp((D1-bond_length)/b)

def get_outer_electrons(atomic_symbol, oxidation_state):
    return int(element(atomic_symbol).group.group_id) - oxidation_state

def check_d_zero(atomic_symbol, oxidation_state):
    if get_outer_electrons(atomic_symbol, oxidation_state) < 1:
        return True
    else:
        return False

def find_zero_eng(engs):
    for energy in engs:
        if energy >= 0:
            zero_eng = engs.index(energy)
            break
    return zero_eng

def approximate_oxygen_position(M1, M2, d):
    a1 = cohp_fit_dict[M1]['slope']
    a2 = cohp_fit_dict[M2]['slope']
    def bl_f(r):
        return bl_derivative(r,a1,a2)
    rt = optimize.root(bl_f, 0.5)
    rt = rt.x[0] * d
    return rt

def bl_derivative(r, a1, a2):
    return a1/r ** 7 - a2/(1-r) ** 7

def bl(r, a1, a2):
    return a1/r ** 6 + a2/(d - r) ** 6

def bs(r, a, b):
    return a/r ** 6 + b

def appromixate_O_3_xy(c1, c2, c3, s1, s2, s3):
    #make the initial guess right in the middle
    initial_guess = [np.mean([a[0] for a in [c1,c2,c3]]), np.mean([a[1] for a in [c1,c2,c3]]) ]
    print(initial_guess)
    def bs_3_xy(xy):
        x = xy[0]
        y = xy[1]
        bs1 = bs(np.sqrt((c1[0]-x)**2 + (c1[1]-y)**2), cohp_fit_dict[s1]['slope'], cohp_fit_dict[s1]['intercept'])
        bs2 = bs(np.sqrt((c2[0]-x)**2 + (c2[1]-y)**2), cohp_fit_dict[s2]['slope'], cohp_fit_dict[s2]['intercept'])
        bs3 = bs(np.sqrt((c3[0]-x)**2 + (c3[1]-y)**2), cohp_fit_dict[s3]['slope'], cohp_fit_dict[s3]['intercept'])
        return -1 * (bs1 + bs2 + bs3)
    opt = optimize.minimize(bs_3_xy, initial_guess)
    x = opt.x[0]
    y = opt.x[1]
    print(np.sqrt((c1[0]-x)**2 + (c1[1]-y)**2), np.sqrt((c2[0]-x)**2 + (c2[1]-y)**2), np.sqrt((c3[0]-x)**2 + (c3[1]-y)**2))
    return opt

def get_potential(filename, metal_symbol):
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if ' {} '.format(metal_symbol) in line:
        #if ' O ' in line:
            metal_potential = float(line.split()[2])
        if ' O ' in line:
            O_potential = float(line.split()[2])

    return {'{}'.format(metal_symbol):metal_potential,
            'O':O_potential}

def get_average_bader_charges(filename):
    atoms = read(filename)
    chrg_dict = {}

    for sym in set(atoms.get_chemical_symbols()):
        n_sp = atoms.get_chemical_symbols().count(sym)
        chrg_dict[sym] = {'n':n_sp, 'average_charge':0}
    for atom in atoms:
        sym = atom.symbol
        chrg = atom.charge
        chrg_dict[sym]['average_charge'] += chrg
    for sym in chrg_dict.keys():
        chrg_dict[sym]['average_charge'] /= chrg_dict[sym]['n']
        del chrg_dict[sym]['n']
    return chrg_dict


def get_dos_center(filename, bound_eng=None):
    with open(filename, 'r') as f:
        dos = json.load(f)
    for key in dos.keys():
        dos[key] = np.array(dos[key])

    if bound_eng is not None:
        for i, eng in enumerate(dos['energy']):
            if eng > bound_eng:
                si = i
                break
        p_center = np.sum(dos['O_p'][:si]*dos['energy'][:si]) / np.sum(dos['O_p'][:si])
    else:
        p_center = np.sum(dos['O_p']*dos['energy']) / np.sum(dos['O_p'])
    return p_center

def get_center(x, engs):
    return np.sum(x*engs) / np.sum(x)

# functions take from: https://gist.github.com/tgmaxson/8b9d8b40dc0ba4395240



def get_coordination_numbers(atoms, covalent_percent=1.25):
    """Returns an array of coordination numbers and an array of existing bonds determined by
    distance and covalent radii.  By default a bond is defined as 120% of the combined radii
    or less. This can be changed by setting 'covalent_percent' to a float representing a
    factor to multiple by (default = 1.2).
    If 'exclude' is set to an array,  these atomic numbers with be unable to form bonds.
    This only excludes them from being counted from other atoms,  the coordination
    numbers for these atoms will still be calculated,  but will be unable to form
    bonds to other excluded atomic numbers.
    """

    # Get all the distances
    distances = np.divide(atoms.get_all_distances(mic=True), covalent_percent)

    # Atomic Numbers
    numbers = atoms.numbers
    # Coordination Numbers for each atom
    cn = []
    cr = np.take(CR, numbers)
    # Array of indices of bonded atoms.  len(bonded[x]) == cn[x]
    bonded = []
    indices = list(range(len(atoms)))
    for i in indices:
        bondedi = []
        for ii in indices:
            # Skip if measuring the same atom
            if i == ii:
                continue
            if (cr[i] + cr[ii]) >= distances[i,ii]:
                bondedi.append(ii)
        # Add this atoms bonds to the bonded list
        bonded.append(bondedi)
    for i in bonded:
        cn.append(len(i))
    return cn, bonded


def get_gcn(site, cn, bond, surface_type="fcc"):
    """Returns the generalized coordination number of the site given.  To define
    A site,  just give a list containing the indices of the atoms forming the site.
    The calculated coordination numbers and bonding needs to be supplied, and the
    surface type also needs to be changed if the surface is not represented by bulk fcc.
    """
    # Define the types of bulk accepted
    gcn_bulk = {"fcc": [12., 18., 22., 26.], "bcc": [14., 22., 28., 32.]}
    gcn = 0.
    if len(site) == 0:
        return 0
    sur = site
    counted = []
    for i in site:
        counted.append(i)
    for i in sur:
        for ii in bond[i]:
            if ii in counted:
                continue
            counted.append(ii)
            gcn += cn[ii]
    return gcn / gcn_bulk[surface_type][len(site) - 1]


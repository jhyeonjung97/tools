#!/usr/bin/env python

from __future__ import print_function

import os, re
import numpy as np
from ase.io import read
from optparse import OptionParser


import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['axes.unicode_minus'] = False
# Set font to Helvetica
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['font.sans-serif'] = ['Helvetica']

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Polygon
from matplotlib.font_manager import FontProperties
import matplotlib.colors as mcolors
plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['font.sans-serif'] = ['Helvetica']

############################################################
__version__ = "1.0"
############################################################


def parseList(string):
    def parseRange(rng):
        # print(rng)
        m = re.match(r'(\d+)(?:[-:](\d+))?(?:[-:](\d+))?$', rng)
        if not m:
            raise ValueError(
                """
The band index should be assigned with combination of the following ways:
-> 10, a single band with index 10
-> 20:30, or '20-30', a continuous range from 20 to 30, 30 included
-> 30:50:2, or '30-50:2', a continues range from 30 to 50, with step size 2
-> '1 2 3', all the patterns above in one string separated by spaces.
For example: '1 4:6 8:16:2' will be converted to '1 4 5 6 8 10 12 14 16'
"""
            )
        ii = m.group(1)
        jj = m.group(2) or ii
        ss = m.group(3) or 1
        return [x-1 for x in range(int(ii), int(jj)+1, int(ss))]

    ret = []
    for rng in string.split():
        ret += parseRange(rng)
    return list(set(ret))


def getElemIdx(poscar='POSCAR'):
    '''
    Get the atom indices corresponding to elements in POSCAR
    - poscar: POSCAR file name, default is 'POSCAR'

    returns a dictionary with keys are symbols, vals are list of
    corresponding indices, starting from 0

    e.g. {  'Mo': [0, 1],
            'S': [2, 3, 4, 5],
            'Se': [6, 7],
            'W': [8]}

    Ionizing
    '''
    pos = read(poscar)
    symb_num = {}
    for idx, elem in enumerate(pos.get_chemical_symbols()):
        if elem not in symb_num:
            symb_num[elem] = [idx]
        else:
            symb_num[elem].append(idx)

    return symb_num


def parseSpdProjection(spd):
    '''
    Parse spdProjections string.  str -> [int]

    # Ionizing
    '''
    spd_dict = {
            's'   : [0],
            'p'   : [1, 2, 3],
            'd'   : [4, 5, 6, 7, 8],
            'f'   : [9, 10, 11, 12, 13, 14, 15],
            'py'  : [1],
            'pz'  : [2],
            'px'  : [3],
            'dxy' : [4],
            'dyz' : [5],
            'dz2' : [6],
            'dxz' : [7],
            'dx2' : [8],
    "fy(3x2-y2)"  : [9],
    "fxyz  "      : [10],
    "fyz2  "      : [11],
    "fz3   "      : [12],
    "fxz2  "      : [13],
    "fz(x2-y2)"   : [14],
    "fx(x2-3y2) " : [15],
    }

    ret = []
    for l in spd.split():
        try:
            assert int(l) <= 15, "Maximum spd index should be <= 15."
            ret += [int(l)]
        except:
            if l.lower() not in spd_dict:
                raise ValueError(
                   "Spd-projected wavefunction character of each KS orbital.\n"
                   "    s orbital: 0\n"
                   "    py, pz, px orbital: 1 2 3\n"
                   "    dxy, dyz, dz2, dxz, dx2 orbital: 4 5 6 7 8 \n"
                   "    fy(3x2-y2), fxyz, fyz2, fz3, fxz2, fz(x2-y2), fx(x2-3y2) orbital: 9 10 11 12 13 14 15\n"
                   "\nFor example, --spd 's dxy 10' specifies the s/dxy/fxyz components\n"
                )
            ret += spd_dict[l]

    return list(set(ret))


def WeightFromPro(infile='PROCAR', lsorbit=False):
    """
    Contribution of selected atoms to the each KS orbital
    """

    assert os.path.isfile(infile), '%s cannot be found!' % infile
    FileContents = [line for line in open(infile) if line.strip()]

    # when the band number is too large, there will be no space between ";" and
    # the actual band number. A bug found by Homlee Guo.
    # Here, #kpts, #bands and #ions are all integers
    nkpts, nbands, nions = [int(xx) for xx in re.sub(
        '[^0-9]', ' ', FileContents[1]).split()]

    # Weights = np.asarray([line.split()[-1] for line in FileContents
    #                       if not re.search('[a-zA-Z]', line)], dtype=float)
    Weights = np.asarray([line.split()[1:-1] for line in FileContents
                          if not re.search('[a-zA-Z]', line)], dtype=float)

    kpt_weight = np.asarray(
        [line.split()[-1] for line in FileContents if 'weight' in line], dtype=float)

    energies = np.asarray([line.split()[-4] for line in FileContents
                           if 'occ.' in line], dtype=float)

    nlmax = Weights.shape[-1]
    nspin = Weights.shape[0] // (nkpts * nbands * nions)
    nspin //= 4 if lsorbit else 1

    if lsorbit:
        Weights.resize(nspin, nkpts, nbands, 4, nions, nlmax)
        Weights = Weights[:, :, :, 0, :, :]
    else:
        Weights.resize(nspin, nkpts, nbands, nions, nlmax)

    kpt_weight.resize(nspin, nkpts)
    energies.resize(nspin, nkpts, nbands)

    return energies, kpt_weight, Weights

############################################################


def gradient_fill(x, y, fill_color=None, ax=None, direction=1, **kwargs):
    """
    Plot a line with a linear alpha gradient filled beneath it.

    Parameters
    ----------
    x, y : array-like
        The data values of the line.
    fill_color : a matplotlib color specifier (string, tuple) or None
        The color for the fill. If None, the color of the line will be used.
    ax : a matplotlib Axes instance
        The axes to plot on. If None, the current pyplot axes will be used.
    Additional arguments are passed on to matplotlib's ``plot`` function.

    Returns
    -------
    line : a Line2D instance
        The line plotted.
    im : an AxesImage instance
        The transparent gradient clipped to just the area beneath the curve.
    """

    line, = ax.plot(x, y, **kwargs)
    if fill_color is None:
        fill_color = line.get_color()

    # print fill_color
    zorder = line.get_zorder()
    alpha = line.get_alpha()
    alpha = 1.0 if alpha is None else alpha

    z = np.empty((100, 1, 4), dtype=float)
    rgb = mcolors.colorConverter.to_rgb(fill_color)
    z[:, :, :3] = rgb
    if direction == 1:
        z[:, :, -1] = np.linspace(0, alpha, 100)[:, None]
    else:
        z[:, :, -1] = np.linspace(alpha, 0, 100)[:, None]

    xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()
    im = ax.imshow(z, aspect='auto', extent=[xmin, xmax, ymin, ymax],
                   origin='lower', zorder=zorder)

    xy = np.column_stack([x, y])
    if direction == 1:
        xy = np.vstack([[xmin, ymin], xy, [xmax, ymin], [xmin, ymin]])
    else:
        xy = np.vstack([[xmin, ymax], xy, [xmax, ymax], [xmin, ymax]])
    clip_path = Polygon(xy, lw=0.0, facecolor='none',
                        edgecolor='none', closed=True)
    ax.add_patch(clip_path)
    im.set_clip_path(clip_path)

    ax.autoscale(True)

    return line, im
############################################################


def lorentz_smearing(x, x0, sigma=0.03):
    '''
    Lorentz smearing of a Delta function.
    '''
    return sigma / np.pi / ((x-x0)**2 + sigma**2)


def gaussian_smearing(x, x0, sigma=0.05):
    '''
    Gaussian smearing of a Delta function.
    '''
    smear = np.zeros(x.size)
    condition = np.abs(x - x0) < (sigma * 5.)
    smear[condition] = 1. / (np.sqrt(2*np.pi) * sigma) * \
        np.exp(-(x[condition] - x0)**2 / (2*sigma**2))

    return smear


def gaussian_smearing_org(x, x0, sigma=0.05):
    '''
    Gaussian smearing of a Delta function.
    '''

    return 1. / (np.sqrt(2*np.pi) * sigma) * np.exp(-(x - x0)**2 / (2*sigma**2))


def generateDos(opts):
    '''
    generate dos
    '''

    ens, kptw, whts = WeightFromPro(opts.procar, opts.lsorbit)
    nspin, nkpts, nbands, nions, nlmax = whts.shape

    emin = ens.min()
    emax = ens.max()
    eran = emax - emin
    emin = emin - eran * opts.extra
    emax = emax + eran * opts.extra

    xen = np.linspace(emin, emax, opts.nedos)
    # tDOS = np.empty((opts.nedos, nspin))
    pDOS = []

    # make all the k-points weight equal
    if opts.homoKpts:
        kptw[...] = 1.0

    tdos_smear = np.empty((nspin, nkpts, nbands, opts.nedos))
    for IS in range(nspin):
        sign = 1 if IS == 0 else -1
        for Ik in range(nkpts):
            for Ib in range(nbands):
                x0 = ens[IS, Ik, Ib]
                # tdos_smear[IS, Ik, Ib] = sign * gaussian_smearing(xen, x0, opts.sigma) * kptw[IS, Ik]
                tdos_smear[IS, Ik, Ib] = sign * \
                    gaussian_smearing_org(xen, x0, opts.sigma) * kptw[IS, Ik]
                # tdos_smear[IS, Ik, Ib] = sign * lorentz_smearing(xen, x0, opts.sigma) * kptw[IS, Ik]
    tDOS = np.sum(tdos_smear, axis=(1, 2)).T

    # tdos_smear = np.empty((nspin, nkpts, nbands, opts.nedos))
    # for IS in range(nspin):
    #     sign = 1 if IS == 0 else -1
    #     for Ik in range(nkpts):
    #         tdos_smear[IS, Ik, :] = sign * gaussian_smearing_org(xen[np.newaxis, ...],
    #                                         ens[IS,Ik,:,np.newaxis], opts.sigma) * kptw[IS, Ik]
    # tDOS = np.sum(tdos_smear, axis=(1, 2)).T

    # tdos_smear = gaussian_smearing_org(xen[np.newaxis,...], ens[...,np.newaxis],
    #         opts.sigma) * kptw[..., np.newaxis, np.newaxis]
    # tDOS = np.sum(tdos_smear, axis=(1, 2)).T
    # if nspin == 2: tDOS[:,1] *= -1

    if len(opts.elem_list) != 0:
        elem_idx = getElemIdx(opts.posfile)
        for elem in opts.elem_list:
            if elem.upper() == 'XX':
                opts.pdosAtom = [
                    ' '.join([str(i+1) for i in elem_idx[k]]) for k in elem_idx]
                # print(opts.pdosAtom)
                opts.pdosLabel = list(elem_idx.keys())
                break
            opts.pdosAtom.append(" ".join([str(i+1) for i in elem_idx[elem]]))
            opts.pdosLabel.append(elem)

    if len(opts.pdosAtom) != 0:

        # add a factor for each PDOS
        factor = np.ones(len(opts.pdosAtom))
        for ii in range(min(len(opts.pdosAtom), len(opts.pdosFactor))):
            factor[ii] = opts.pdosFactor[ii]

        # only plot DOS from selected kpoints, default from all kpoints
        selected_kpts_index = [slice(None, None)
                               for ii in range(len(opts.pdosAtom))]
        for ii in range(min(len(opts.pdosAtom), len(opts.pdosKpts))):
            alist = np.array(opts.pdosKpts[ii].split(), dtype=int)
            # index starting from 1
            assert np.max(
                alist) <= nkpts, "Maximum kpoint index must <{:d}".format(nkpts)
            assert np.min(alist) >= 1, "Minimum kpoint index must >=1"

            nlist = [x for x in alist if not x == -1]
            cmark, = np.where(alist == -1)
            for ii in cmark:
                nlist += range(alist[ii + 1], alist[ii + 2] + 1)
            # index starting from 0
            nlist = [x - 1 for x in set(nlist)]

            selected_kpts_index[ii] = nlist

        for ia, atoms in enumerate(opts.pdosAtom):
            # p = np.zeros((opts.nedos, nspin))

            #  alist = np.array(atoms.split(), dtype=int)
            if '0' in atoms.split():
                nlist = range(nions)
            else:
                nlist = parseList(atoms)
                # nlist = [x for x in alist if not x == -1]
                # cmark, = np.where(alist == -1)
                # for ii in cmark:
                #     nlist += range(alist[ii + 1], alist[ii + 2] + 1)
                # nlist = [x - 1 for x in set(nlist)]

            if ia <= len(opts.spdProjections) - 1:
                # spdList = [int(x) for x in opts.spdProjections[ia].split()]
                spdList = parseSpdProjection(opts.spdProjections[ia])   # Ionizing
                pwhts = np.sum(whts[..., spdList], axis=-1)
            else:
                pwhts = np.sum(whts, axis=-1)

            pwhts = np.sum(pwhts[:, :, :, nlist], axis=-1)

            # p = np.sum(pwhts[..., np.newaxis] * tdos_smear, axis=(1, 2)).T
            p = np.sum(pwhts[:, selected_kpts_index[ia], :, np.newaxis] *
                       tdos_smear[:, selected_kpts_index[ia], ...], axis=(1, 2)).T

            for IS in range(nspin):
                sign = 1 if IS == 0 else -1

                p[:, IS] = opts.pdosOffset * ia * sign + p[:, IS] * factor[ia]
                if ia == 0:
                    tDOS[:, IS] += sign * opts.pdosOffset * len(opts.pdosAtom)

            pDOS += [p]

    return xen, tDOS, pDOS

############################################################


def readDOSFromFile(opts):
    '''
    Read DOS info from file.

    the format of the DOS file:
    first line: ISPIN, NEDOS
    second line: labels of the dos
    next lines:
        if ISPIN = 1:
            Energy pDOS1 PDOS2 ... TotalDOS
        else:
            Energy pDOS1_up PDOS2_up ... TotalDOS_up pDOS1_down PDOS2_down ... TotalDOS_down
    '''

    inp = open(opts.dosFromFile).readlines()

    # the dos basic info
    nspin, nedos = [int(x) for x in inp[0].split()[1:]]
    labels = inp[1].split()[1:]
    # data
    DOS = np.array([line.split() for line in inp[2:] if line.strip()],
                   dtype=float)
    NoPdos = (DOS.shape[1] - 1) // nspin - 1

    tDOS = np.empty((nedos, nspin))
    pDOS = []
    xen = DOS[:, 0]
    for ii in range(nspin):
        tDOS[:, ii] = DOS[:, (ii + 1) * (NoPdos + 1)]
    for pp in range(NoPdos):
        tmp = []
        for ii in range(nspin):
            tmp += [DOS[:, (pp + 1) + ii * (NoPdos + 1)]]
        pDOS += [np.array(tmp).T]

    opts.nedos = nedos
    opts.pdosAtom = ['' for x in range(NoPdos)]
    opts.pdosLabel = labels

    return xen, tDOS, pDOS


def saveDOSToFile(opts, xen, tDOS, pDOS):
    '''
    save DOS info to file.

    the format of the DOS file:
    first line: ISPIN, NEDOS
    second line: labels of the dos
    next lines:
        if ISPIN = 1:
            Energy pDOS1 PDOS2 ... TotalDOS
        else:
            Energy pDOS1_up PDOS2_up ... TotalDOS_up pDOS1_down PDOS2_down ... TotalDOS_down
    '''

    nspin = tDOS.shape[1]
    nedos = tDOS.shape[0]
    NoPdos = len(pDOS)

    out = open(opts.dosToFile, 'w')
    out.write('# %5d %8d\n' % (nspin, nedos))
    labels = '# ' + ' '.join(opts.pdosLabel) + '\n'
    out.write(labels)

    for nn in range(nedos):
        line = '%8.4f ' % xen[nn]
        for ii in range(nspin):
            for p in pDOS:
                line += '%8.4f ' % p[nn, ii]
            line += '%8.4f ' % tDOS[nn, ii]
        line += '\n'
        out.write(line)

############################################################


def plot_single_dos(ax, xen, tdos, pdos, opts, vlines1, vlines2, labels1, labels2):
    '''
    Plot DOS on a single axes
    '''
    LINES = []
    nspin = tdos.shape[1]

    plabels = []
    LWs = []
    LCs = []
    if opts.pdosAtom:
        plabels = ['p_%d' % ii for ii in range(len(opts.pdosAtom))]
        LWs = [0.5 for ii in range(len(opts.pdosAtom))]
        LCs = [None for ii in range(len(opts.pdosAtom))]

        for ii in range(min(len(opts.pdosAtom), len(opts.pdosLabel))):
            plabels[ii] = opts.pdosLabel[ii]
        for ii in range(min(len(opts.pdosAtom), len(opts.linewidth))):
            LWs[ii] = opts.linewidth[ii]
        for ii in range(min(len(opts.pdosAtom), len(opts.linecolors))):
            LCs[ii] = opts.linecolors[ii]

    plabels += ['total']

    xen_plot = xen - opts.zero
    line = None
    for ip, p in enumerate(pdos):
        for ii in range(nspin):
            fill_direction = 1 if ii == 0 else -1
            if ii == 0:
                lc = LCs[ip] if LCs[ip] is not None else None
            else:
                lc = line.get_color() if line is not None else None

            if opts.fill:
                line, im = gradient_fill(xen_plot, p[:, ii], ax=ax, lw=LWs[ip],
                                         color=lc,
                                         direction=fill_direction)
            else:
                line, = ax.plot(xen_plot, p[:, ii], lw=LWs[ip], alpha=0.6,
                                color=lc)
            if ii == 0:
                LINES += [line]

    if opts.showtotal:
        for ii in range(nspin):
            fill_direction = 1 if ii == 0 else -1
            if ii == 0:
                lc = 'k'
            else:
                lc = line.get_color() if line is not None else None

            if opts.fill:
                line, im = gradient_fill(xen_plot, tdos[:, ii], ax=ax,
                                         color=lc,
                                         lw=0.5,
                                         # zorder=-1,
                                         direction=fill_direction,
                                         )
            else:
                line, = ax.plot(xen_plot, tdos[:, ii], color=lc,
                                lw=0.5, alpha=0.6)
            if ii == 0:
                LINES += [line]

    # Create semibold font for labels
    # Arial doesn't support intermediate weights, so we'll use 'normal' with slight adjustment
    # or try 'semibold' if the font supports it
    try:
        semibold_font = FontProperties(weight='semibold')
    except:
        # Fallback to normal if semibold not supported
        semibold_font = FontProperties(weight='normal')
    ax.set_ylabel('DOS [a.u.]',  # fontsize='small',
                  labelpad=10, fontproperties=semibold_font)
    ax.tick_params(which='both', labelsize='small', direction='in')

    xmin, xmax = opts.xlim
    ax.set_xlim(xmin, xmax)
    if opts.ylim is not None:
        ymin, ymax = opts.ylim
        ax.set_ylim(ymin, ymax)

    ax.set_yticks([])
    ax.axhline(0.0, color='black', lw=0.7)
    ax.axvline(0.0, color='black', lw=0.7)
    ax.axvline(vlines1, color='C2', lw=0.7)
    ax.axvline(vlines2, color='C1', lw=0.7)
    ax.text(vlines1+0.1, ax.get_ylim()[1]*0.5, labels1, color='C2', fontsize=8, rotation=90)
    ax.text(vlines2+0.1, ax.get_ylim()[1]*0.5, labels2, color='C1', fontsize=8, rotation=90)

    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    opts.pdosLabel = plabels
    leg = ax.legend(LINES, plabels,
                     loc=opts.legendloc,
                     fontsize='small',
                     frameon=True,
                     framealpha=1.0,
                     labelspacing=0.3)  # Adjust vertical spacing between legend entries
    
    # Make lines in legend thicker
    for line in leg.get_lines():
        line.set_linewidth(1.0)

    return LINES


def dosplot_bulk(xen1, tdos1, pdos1, opts1, xen2, tdos2, pdos2, opts2):
    '''
    Use matplotlib to plot two DOS plots in one figure, stacked vertically with shared x-axis
    '''

    width, height = opts1.figsize
    # Double the height for two plots
    height = height * 2
    dpi = opts1.dpi

    plt.style.use(opts1.mpl_style)
    # DO NOT use unicode minus regardless of the style
    mpl.rcParams['axes.unicode_minus'] = False

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(width, height))
    
    # Plot first DOS (top, pydos-ruo2 style)
    plot_single_dos(ax1, xen1, tdos1, pdos1, opts1, 
                    -2.4773, -3.7085, '-2.48', '-3.71')

    # Plot second DOS (bottom, pydos-reruo2 style)
    plot_single_dos(ax2, xen2, tdos2, pdos2, opts2,
                    -2.5551, -3.3965, '-2.56', '-3.40')

    # Set xlabel only on bottom plot
    try:
        semibold_font = FontProperties(weight='semibold')
    except:
        # Fallback to normal if semibold not supported
        semibold_font = FontProperties(weight='normal')
    ax2.set_xlabel('Energy [eV]',  # fontsize='small',
                  labelpad=5, fontproperties=semibold_font)

    # Remove x-axis ticks and labels from top plot
    ax1.tick_params(labelbottom=False, bottom=True, top=False, direction='in')
    # Remove top ticks from bottom plot
    ax2.tick_params(top=False, direction='in')
    
    # Adjust subplots to be tightly packed with no space between them
    fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1, hspace=0.0, wspace=0.0)
    
    # Manually set positions to ensure they are exactly touching
    # Get current positions to preserve x coordinates
    pos1 = ax1.get_position()
    pos2 = ax2.get_position()
    
    # Make them exactly touch with equal heights
    # Total available height: 0.95 - 0.10 = 0.85
    # Each plot gets half: 0.425
    # Top plot: from 0.525 to 0.95 (height = 0.425)
    # Bottom plot: from 0.10 to 0.525 (height = 0.425)
    plot_height = 0.30
    ax1.set_position([pos1.x0, 0.10 + plot_height, pos1.width, plot_height])  # top plot
    ax2.set_position([pos2.x0, 0.10, pos2.width, plot_height])  # bottom plot
    
    plt.savefig(opts1.dosimage, dpi=dpi, transparent=True)

############################################################


def command_line_arg():
    usage = "usage: %prog [options] arg1 arg2"
    par = OptionParser(usage=usage, version=__version__)

    par.add_option("-i", '--input',
                   action='store', type="string", dest='procar',
                   default='PROCAR',
                   help='location of the first PROCAR')

    par.add_option('--input2',
                   action='store', type="string", dest='procar2',
                   default=None,
                   help='location of the second PROCAR')

    par.add_option("-d", '--dir1',
                   action='store', type="string", dest='dir1',
                   default=None,
                   help='first directory path (PROCAR will be read from this directory)')

    par.add_option('--dir2',
                   action='store', type="string", dest='dir2',
                   default=None,
                   help='second directory path (PROCAR will be read from this directory)')

    par.add_option("-p", '--pdos',
                   action='append', type="string", dest='pdosAtom',
                   default=[],
                   help='specify which atoms to plot the pdos, 0 for all atoms (for first plot)')

    par.add_option('--pdos2',
                   action='append', type="string", dest='pdosAtom2',
                   default=[],
                   help='specify which atoms to plot the pdos, 0 for all atoms (for second plot)')

    par.add_option("-k", '--kpts',
                   action='append', type="string", dest='pdosKpts',
                   default=[],
                   help='specify which k-points to plot the pdos')

    par.add_option("--homokpts",
                   action='store_true', dest='homoKpts',
                   default=False,
                   help='Homogeneous k-points weights')

    par.add_option('--pdosoffset',
                   action='store', type="float", dest='pdosOffset',
                   default=0.0,
                   help='offset in pdos plot')

    par.add_option("-l", '--label',
                   action='append', type="string", dest='pdosLabel',
                   default=[],
                   help='label of the pdos (for first plot)')

    par.add_option('--label2',
                   action='append', type="string", dest='pdosLabel2',
                   default=[],
                   help='label of the pdos (for second plot)')

    par.add_option('--lloc',
                   action='store', type="string", dest='legendloc',
                   default='upper right',
                   help='legend location of dos plot')

    par.add_option('--fac',
                   action='append', type="float", dest='pdosFactor',
                   default=[],
                   help='scale factor of the pdos')

    par.add_option('-z', '--zero',
                   action='store', type="float",
                   dest='zero', default=0.0,
                   help='energy reference of the band plot (for first plot)')

    par.add_option('--zero2',
                   action='store', type="float",
                   dest='zero2', default=None,
                   help='energy reference of the band plot (for second plot, defaults to first plot value)')

    par.add_option('--sigma',
                   action='store', type="float",
                   dest='sigma', default=0.05,
                   help='smearing parameter, default 0.05')

    par.add_option('-n', '--nedos',
                   action='store', type="int",
                   dest='nedos', default=5000,
                   help='number of point in DOS plot')

    par.add_option('-o', '--output',
                   action='store', type="string", dest='dosimage',
                   default='dos-bulk.pdf',
                   help='output image name, "dos-bulk.pdf" by default')

    par.add_option('-s', '--size', nargs=2,
                   action='store', type="float", dest='figsize',
                   default=(5.0, 2.0),
                   help='figure size of the output plot')

    par.add_option('-x', nargs=2,
                   action='store', type="float", dest='xlim',
                   default=(-6, 6),
                   help='x limit of the dos plot')

    par.add_option('-y', nargs=2,
                   action='store', type="float", dest='ylim',
                   default=None,
                   help='energy range of the band plot')

    par.add_option('-e',
                   action='store', type="float", dest='extra',
                   default=0.05,
                   help='extra energy range of the band plot')

    par.add_option('--lw',
                   action='append', type="float", dest='linewidth',
                   default=[],
                   help='linewidth of the band plot')

    par.add_option('--lc',
                   action='append', type="string", dest='linecolors',
                   default=[],
                   help='linecolors of the band plot (for first plot)')

    par.add_option('--lc2',
                   action='append', type="string", dest='linecolors2',
                   default=[],
                   help='linecolors of the band plot (for second plot)')

    par.add_option('--fill',
                   action='store_true', dest='fill',
                   default=True,
                   help='fill under the DOS')

    par.add_option('--nofill',
                   action='store_false', dest='fill',
                   help='no fill under the DOS')

    par.add_option('--dpi',
                   action='store', type="int", dest='dpi',
                   default=360,
                   help='resolution of the output image')

    par.add_option('--tot',
                   action='store_true', dest='showtotal',
                   default=True,
                   help='show total dos')

    par.add_option('--notot',
                   action='store_false', dest='showtotal',
                   help='not show total dos')

    par.add_option('--style',
                   action='store', type='string', dest='mpl_style',
                   default='default',
                   help='plot style of matplotlib. See "plt.style.available" for list of available styles.')

    par.add_option('--fromfile',
                   action='store', type='string', dest='dosFromFile',
                   default=None,
                   help='plot the dos contained in the file')

    par.add_option('--tofile',
                   action='store', type='string', dest='dosToFile',
                   default=None,
                   help='save DOS to file.')

    par.add_option('--spd',
                   action='append', type="string", dest='spdProjections',
                   default=[],
                   help="Spd-projected wavefunction character of each KS orbital (for first plot).\n"
                   "    s orbital: 0\n"
                   "    py, pz, px orbital: 1 2 3\n"
                   "    dxy, dyz, dz2, dxz, dx2 orbital: 4 5 6 7 8 \n"
                   "    fy(3x2-y2), fxyz, fyz2, fz3, fxz2, fz(x2-y2), fx(x2-3y2) orbital: 9 10 11 12 13 14 15\n"
                   "\nFor example, --spd 's dxy 9\n"
                   )

    par.add_option('--spd2',
                   action='append', type="string", dest='spdProjections2',
                   default=[],
                   help="Spd-projected wavefunction character of each KS orbital (for second plot).\n"
                   "    s orbital: 0\n"
                   "    py, pz, px orbital: 1 2 3\n"
                   "    dxy, dyz, dz2, dxz, dx2 orbital: 4 5 6 7 8 \n"
                   "    fy(3x2-y2), fxyz, fyz2, fz3, fxz2, fz(x2-y2), fx(x2-3y2) orbital: 9 10 11 12 13 14 15\n"
                   "\nFor example, --spd2 's dxy 9\n"
                   )

    par.add_option('--lsorbit',
                   action='store_true', dest='lsorbit',
                   help='Spin orbit coupling on, special treament of PROCAR')

    par.add_option('-q', '--quiet',
                   action='store_true', dest='quiet',
                   help='not show the resulting image')

    par.add_option('--elem',
                   action='append', type='string', dest='elem_list',
                   default=[],
                   help='display the PDOS of given element(s), if --elem == "XX", show all of them respectively. (for first plot)')  # Ionizing

    par.add_option('--elem2',
                   action='append', type='string', dest='elem_list2',
                   default=[],
                   help='display the PDOS of given element(s), if --elem2 == "XX", show all of them respectively. (for second plot)')  # Ionizing

    par.add_option('--poscar',
                   action='store', type='string', dest='posfile',
                   default='POSCAR',
                   help='specify which poscar to read, if "--elem" is not empty (for first plot)')  # Ionizing

    par.add_option('--poscar2',
                   action='store', type='string', dest='posfile2',
                   default='POSCAR',
                   help='specify which poscar to read, if "--elem2" is not empty (for second plot)')  # Ionizing

    return par.parse_args()


############################################################
if __name__ == '__main__':
    from time import time
    import copy
    import os
    opts, args = command_line_arg()

    # Handle directory paths
    if opts.dir1:
        if not os.path.isabs(opts.procar):
            opts.procar = os.path.join(opts.dir1, opts.procar)
        if not os.path.isabs(opts.posfile):
            opts.posfile = os.path.join(opts.dir1, opts.posfile)
    
    if opts.dir2:
        if opts.procar2 is None:
            opts.procar2 = 'PROCAR'
        if not os.path.isabs(opts.procar2):
            opts.procar2 = os.path.join(opts.dir2, opts.procar2)
        if not os.path.isabs(opts.posfile2):
            opts.posfile2 = os.path.join(opts.dir2, opts.posfile2)
    elif opts.procar2 is None:
        opts.procar2 = 'PROCAR2'

    # Create options for second PROCAR (copy from first)
    opts2 = copy.deepcopy(opts)
    opts2.procar = opts.procar2
    opts2.posfile = opts.posfile2
    
    # Override with second plot specific options if provided
    if opts.pdosAtom2:
        opts2.pdosAtom = opts.pdosAtom2
    if opts.pdosLabel2:
        opts2.pdosLabel = opts.pdosLabel2
    if opts.spdProjections2:
        opts2.spdProjections = opts.spdProjections2
    if opts.linecolors2:
        opts2.linecolors = opts.linecolors2
    if opts.elem_list2:
        opts2.elem_list = opts.elem_list2
    if opts.zero2 is not None:
        opts2.zero = opts.zero2

    # Check if PROCAR files exist
    if not os.path.isfile(opts.procar):
        print("Error: First PROCAR file not found: %s" % opts.procar)
        exit(1)
    if not os.path.isfile(opts2.procar):
        print("Error: Second PROCAR file not found: %s" % opts2.procar)
        print("Please check if the file exists or specify the correct path with --input2 or --dir2")
        exit(1)

    if opts.dosFromFile:
        xen1, tdos1, pdos1 = readDOSFromFile(opts)
        # For second file, use same function but with opts2
        if hasattr(opts2, 'dosFromFile') and opts2.dosFromFile:
            xen2, tdos2, pdos2 = readDOSFromFile(opts2)
        else:
            # If no second file specified, use same data (for testing)
            xen2, tdos2, pdos2 = xen1, tdos1, pdos1
    else:
        t0 = time()
        xen1, tdos1, pdos1 = generateDos(opts)
        t1 = time()
        print('First DOS calc completed! Time Used: %.2f [sec]' % (t1 - t0))

        t0 = time()
        xen2, tdos2, pdos2 = generateDos(opts2)
        t1 = time()
        print('Second DOS calc completed! Time Used: %.2f [sec]' % (t1 - t0))

    t0 = time()
    dosplot_bulk(xen1, tdos1, pdos1, opts, xen2, tdos2, pdos2, opts2)
    t1 = time()
    print('DOS plot completed! Time Used: %.2f [sec]' % (t1 - t0))

    # save dos to file
    if opts.dosToFile:
        saveDOSToFile(opts, xen1, tdos1, pdos1)

    if not opts.quiet:
        try:
            from subprocess import call
            call(['feh', '-xdF', opts.dosimage])
        except:
            # do nothing if image view fails
            pass

#! /user/bin/env python
"""
Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from math import pow
from pylab import *
import matplotlib
#matplotlib.use('TkAgg')
#import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.ticker import *
from scipy import *
import subprocess

#Import plot_settings as ps
from matplotlib import rc

# settings size and font for revtex stylesheet

fig_width_pt = 1.8*246.0  # Get this from LaTeX using \showthe\columnwidth
#fig_width_pt *= 300./72 # convert to 300 dpi
inches_per_pt = 1.0/72.27               # Convert pt to inches
#inches_per_pt = 1.0/300               # Convert pt to inches
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
fig_size = [fig_width,fig_height]
fig = plt.figure(figsize=fig_size,dpi=300)


font_size = 9 
tick_font_size = 8
xlabel_pad = 8
ylabel_pad = 18
matplotlib.rcParams['ps.usedistiller'] = 'xpdf'

matplotlib.rcParams['font.size'] = 10
#matplotlib.rcParams['axes.labelsize'] = 2*font_size
matplotlib.rcParams['axes.labelsize'] = font_size
matplotlib.rcParams['legend.fontsize'] = font_size
matplotlib.rcParams['xtick.labelsize'] = tick_font_size
matplotlib.rcParams['ytick.labelsize'] = tick_font_size

#font_default='helvetica'
font_default='cmss'
#font_default='times'

def setfont(font=font_default,unicode=True):
   r"""
   Set Matplotlibs rcParams to use LaTeX for font rendering.
   Revert all changes by calling rcdefault() from matplotlib.

   Parameters:
   -----------
   font: string
       "Helvetica"
       "Times"
       "Computer Modern"

   usetex: Boolean
       Use unicode. Default: False.    
   """

   # Use TeX for all figure text!
   #plt.rc('text', usetex=True)

   font = font.lower().replace(" ","")
   if font == 'times':
       # Times
       font = {'family':'serif', 'serif':['Times']}
       preamble  = r"""
                      \usepackage{color}
                      \usepackage{mathptmx}
                   """
   elif font == 'helvetica':
       # Helvetica
       # set serif, too. Otherwise setting to times and then
       # Helvetica causes an error.
       font = {'family':'sans-serif','sans-serif':['Helvetica'],
               'serif':['cm10']}
       preamble  = r"""
                      \usepackage{color}
                      \usepackage[tx]{sfmath}
                      \usepackage{helvet}
                      \usepackage{sansmath}
                   """
   else:
       # Computer modern serif
       font = {'family':'serif', 'serif':['cm10']}
       preamble  = r"""
                      \usepackage{color}
                   """
   if font == 'cmss':
       # Computer modern sans serif
       font = {'family':'sans-serif', 'serif':['cmss']}
       preamble  = r"""
                      \usepackage{color}
                      \usepackage[tx]{sfmath}
                   """

   if unicode:
       # Unicode for Tex
       #preamble =  r"""\usepackage[utf8]{inputenc}""" + preamble
       # inputenc should be set automatically
       plt.rcParams['text.latex.unicode']=True

   #print font, preamble
   plt.rc('font',**font)
   plt.rcParams['text.latex.preamble'] = preamble



#setfont(font_default,unicode=True)


matplotlib.rcParams['lines.linewidth'] = 1.
#matplotlib.rcParams['ytick.direction'] = 'out'
#matplotlib.rcParams['xtick.direction'] = 'out'


ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

zoom=0.5
d1=3*zoom
d2=4*zoom
xcenter=1.45 #0.65
#ycenter=1.23#2.4
ycenter=0.73#2.4

x1=xcenter-d1#-0.6
x2=xcenter+d1#2.2
y1=ycenter-d2#1#0.5
y2=ycenter+d2#5
ax.axis([x1, x2, y1, y2])
ax.set_xlabel(r'$\Delta$G$_{\sf O}$ - $\Delta$G$_{\sf OH}$ (eV)')
#ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$ -$\Delta$G$_{\sf O}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf OH}$ (eV)')




#delta = 0.025
delta = 0.01
x = np.arange(x1,x2+delta, delta)
y = np.arange(y1,y2+delta, delta)
X, Y = np.meshgrid(x, y)
#Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# difference of Gaussians
#Z = 10.0 * (Z2 - Z1)


#fit=[0.84527288, 3.38026638]
def ooh_oh_scaling(doh):
    #like ambars
    #dooh=0.5*doh  + 3.0		 #O
    #normal one
    #dooh=doh*1.0 + 3.2
    #dooh=0.68*doh + 3.55
    dooh=0.84*doh + 3.14
    #dooh=0.66*doh + 3.41
    return dooh


def overpotential(doh,do):
    dooh=ooh_oh_scaling(doh)
    dg14=[doh, do-doh, dooh-do, -dooh+4.92]
    m=max(dg14)
    return m-1.23
    #return doh*do

def overpotential2(x,doh):
    dooh=ooh_oh_scaling(doh)
    dg14=[doh, x, -x+2.46, -dooh+4.92]
    m=max(dg14)
    return m-1.23
    #return doh*do

def overpotential3(x,doh):
    dooh=ooh_oh_scaling(doh)
    dg14=[doh, x, dooh-(x+doh), -dooh+4.92]
    m=max(dg14)
    return m-1.23
    #return doh*do

def overpotential_label(doh,do):
    dooh=ooh_oh_scaling(doh)
    dg14=[doh, do-doh, dooh-do, -dooh+4.92]
    m=max(dg14)
    for i in range(len(dg14)):
       if(m==dg14[0]):
          return r'OH lim.'
       if(m==dg14[1]):
          return r'OH-O lim.'
       if(m==dg14[2]):
          return r'O-OOH lim.'
       if(m==dg14[3]):
          return r'OOH-O$_{\sf 2}$ lim.'
    #return doh*do    
def oer_step (i):
        steps=['H2O->OH*','OH*->O*', 'O*->OOH*', 'OOH*->H2O']
        return str(steps[i]) 


def overpotential_oer_full(doh,do,dooh):
    #do=1.556*doh + 0.9951 #from our scaling
    dg23=[doh,do-doh,dooh-do,4.92-dooh] #add solvation soon
    m=max(dg23)
    ##print dg14.index(m)
    #print dg14,m
    print (dg23)
    return [round(m-1.23,2),round(m,2),oer_step(dg23.index(m))]



#Z=overpotential(X,Y)

Z=[]
for j in y:
    tmp=[]
    for i in x:
        tmp.append(overpotential3(i,j))
    Z.append(tmp)


#print overpotential(0.8,2.4)

Z = np.array(Z)


#im = plt.imshow(Z, origin='lower',interpolation='bilinear',
#                cmap=cm.jet_r, extent=(x1,x2,y1,y2), vmin=0, vmax=2)

import numpy as np
import seaborn as sns
#sns.set()

from matplotlib.colors import ListedColormap

origin='lower'
levels = np.arange(0.2, 1.6 ,0.1)
#levels = np.arange(0.2, 2, 0.1)
CS = plt.contourf(X, Y, Z, levels, #20, # [-1, -0.1, 0, 0.1],
                        #alpha=0.8,
                        #cmap=plt.cm.bone,
                        #cmap=plt.cm.coolwarm_r,
                        #cmap=ListedColormap([u'#f07894', u'#f1859e', u'#f392a8', u'#f4a0b3', u'#f6adbe', u'#f7bbc9', u'#f9c8d3', u'#fad6de', u'#fce3e9', u'#f2f2f2', u'#f2f2f2', u'#e6f4f4', u'#d3ebec', u'#bfe1e4', u'#acd9dc', u'#98cfd3', u'#85c6cb', u'#71bdc2', u'#5eb4ba', u'#4bacb2'])  ,
                        #cmap=ListedColormap([u'#800026', u'#bd0026', u'#e31a1c', u'#fc4e2a', u'#fd8d3c', u'#feb24c', u'#fed976', u'#ffeda0', u'#fce3e9', u'#f2f2f2', u'#f2f2f2', u'#e6f4f4', u'#d3ebec', u'#bfe1e4', u'#acd9dc', u'#98cfd3', u'#85c6cb', u'#71bdc2', u'#5eb4ba', u'#fff7fb'])  ,
                        cmap=ListedColormap([
			'#a50026',
			'#d73027',
			'#f46d43',
			'#fdae61',
			'#fee090',
			'#ffffbf',
			#'#ffffff',
			'#ffffe5',
			'#e0f3f8',
			'#abd9e9',
			'#74add1',
			'#4575b4',
			'#313695',
			'#d0d1e6',
			])  ,
                                
                        #cmap=ListedColormap(sns.color_palette("BuGn_r",256))  ,
                        #cmap=ListedColormap(sns.color_palette("cubehelix", 20))  ,
                        #cmap=sns.diverging_palette(0, 200, sep=30, n=20, s=80, l=65,  as_cmap=True)  ,
                        #cmap=sns.color_palette("BuGn_r"),
                        #cmap=plt.cm.jet_r,
                        #extend='both',
                        extend='max',
                        origin=origin)

pal=sns.diverging_palette(0, 200, sep=30, n=20, s=80, l=65).as_hex()
pal_hls = sns.hls_palette(12, l=.3, s=.8).as_hex(); 
'''
colorname= []
colorid= []
import matplotlib
for name, hex in matplotlib.colors.cnames.items():
    colorname.append(name)
    colorid.append(hex)
zippedcolors = list(zip(colorname, colorid))
zippedcolors = sorted(zippedcolors, key=lambda x: x[1])

print zippedcolors
'''
# Note that in the following, we explicitly pass in a subset of
# the contour levels used for the filled contours.  Alternatively,
# We could pass in additional levels to provide extra resolution,
# or leave out the levels kwarg to use all of the original levels.

'''
CS2 = plt.contour(CS, levels=CS.levels,
                        colors = 'white',
                        linewidths=0.05,
                        alpha=0.3,
                        origin=origin,
                        hold='on')
'''

#levels = np.arange(0, 2, 0.05)
#CS = plt.contourf(X,Y,Z, levels, cmap=cm.jet_r, origin='lower')
#CS = plt.contourf(X,Y,Z, levels, origin='lower')

#im = plt.imshow(Z, interpolation='bilinear', origin='lower',
#                cmap=cm.jet, extent=(x1,x2,y1,y2))


#levels2 = [2.0]
#CS2 = plt.contour(CS, levels2,
#                        colors = 'r',
#                        origin='lower',
#                        hold='on')

#CS = plt.contour(Z, levels,
#                 origin='lower',
#                 linewidths=0.5,
#                 extent=(x1,x2,y1,y2))




##Thicken the zero contour.
#zc = CS.collections[6]
#plt.setp(zc, linewidth=2)


cbar = plt.colorbar(CS, ticks=list(np.arange(0.2, 1.6, step=0.1)))
#cbar = plt.colorbar(CS)

#cbar.ax.set_ylabel('Overpotential [V]')
#cbar.ax.set_ylabel(r'$\eta_{\sf calc.}$')
cbar.ax.set_ylabel(r'$\eta_{\sf OER}$ (V)')

cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')


#cbar.add_lines(CS2)


#plt.clabel(CS, levels[1::2],  # label every second level
#           inline=1,
#           fmt='%1.1f',
#           fontsize='x-small')


#plt.title('Lines with colorbar')

# We can still add a colorbar for the image, too.


# This makes the original colorbar look a bit out of place,
# so let's improve its position.

ax.tick_params(axis='both', direction='out')
ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
ax.get_yaxis().tick_left()


#plot(x,ooh_oh_scaling(x),'--',color='orange',lw=1, dashes=(3,1),label='$\Delta$G$_{\sf OOH}$=0.82G$_{\sf OH}$+3.18 eV')

#ax.text(x1+0.02,y2-0.3,'$\Delta$G$_{\sf OOH}$=%.2fG$_{\sf OH}$+%.2f eV' %(fit[0],fit[1]),color='orange',fontsize='x-small',zorder=10,horizontalalignment='left')



#ax.text(x1+0.02,y2-0.3,'$\Delta$G$_{\sf OOH}$=%.2fG$_{\sf OH}$+%.2f eV' %(0.82,3.18),color='orange',fontsize='x-small',zorder=10,horizontalalignment='left')


#plt.show()

offset=[0.0,0.08]




calc_systems=[] 
#energies are dGs on dEs!
#to input dEs
def make_dGs(x,i):
    dgeoh=0.391254
    dgeo=0.0216138
    dgeooh=0.451538
    if(i==0):
       return x+dgeoh
    if(i==1):
       return x+dgeo
    if(i==2):
       return x+dgeooh

#colors
#colors=['#e5d8bd','#ffffcc','#fed9a6','#decbe4','#ccebc5','#fbb4ae','#b3cde3']
colors=['#a6cee3','#fdbf6f','#b2df8a','#1f78b4','#e31a1c','#fb9a99','#33a02c']
#tcolor='#737373'
tcolor='#d9d9d9'
bcolor='#252525'
#bcolor='#99000d'
symbols=['D','^','o','s','h','H','x','*','+']
tcolor2='white'
tcolor3='black'
tcolor4='gray'

#doped stipes
#From Zhenghang Zhao
#OER 
# Ti
#[dgOH, dgO, dgOOH, overpotential, blablabla....]

# Ni

# solvation correction for NiO2 and LNO-NiO2 interfaces ---Jiang Li Feb-18,2020
solv_corr_OH = -0.21
solv_corr_O = -0.003
solv_corr_OOH = -0.06

'''
from DE conversion 

 LNO-Ni-S4-O1c             ,    1.123   ,  3.478  ,   4.307   ,     1.124
 LNO-Ni-S4-O2c             ,    0.824   ,  2.950  ,   4.798   ,     0.895
 LNO-La-S4-O2c             ,    1.355   ,  3.421  ,   4.485   ,     0.835
 NiO2-LNO-La-O3c1          ,    0.484   ,  2.204  ,   3.886   ,     0.490
 NiO2-LNO-La-O3c2          ,    0.615   ,  2.204  ,   3.977   ,     0.543
 NiO2-LNO-La-O3c3          ,    0.450   ,  2.127  ,   3.968   ,     0.611
 NiO2-O3c                  ,    0.527   ,  1.648  ,   3.823   ,     0.945
 LaO2-LNO-La-typeII-O3c    ,    1.015   ,  3.151  ,   4.340   ,     0.906
 LaO2-LNO-La-typeII-O2c    ,    0.751   ,  3.026  ,   3.884   ,     1.045
'''

calc_systems.append([-0.316   , -0.034  ,   3.238 , 0.0  ,r'HT-WO$_3$(001)-6f',tcolor3,0.0,0.08,1.5,'8','lime','o','orangered'])
calc_systems.append([-0.302   , -0.054  ,   3.259 , 0.0  ,r'$\alpha$-WO$_3$(002)-6f',tcolor3,0.0,0.08,1.5,'8','lime','s','orangered'])
calc_systems.append([-0.408   , -0.248  ,   3.013 , 0.0  ,r'$\gamma$-WO$_3$(002)-6f',tcolor3,0.0,0.08,1.5,'8','lime','D','orangered'])
calc_systems.append([-0.588   , -0.354  ,   2.945 , 0.0  ,r'$\gamma$-WO$_3$(020)-6f',tcolor3,0.0,0.08,1.5,'8','lime','d','orangered'])

calc_systems.append([ 2.423   ,  4.577  ,   4.855 , 0.0  ,r'HT-WO$_3$(001)-5f',tcolor3,0.0,0.08,1.5,'8','lime','o','gold'])
calc_systems.append([ 2.473   ,  4.630  ,   4.862 , 0.0  ,r'$\alpha$-WO$_3$(002)-5f',tcolor3,0.0,0.08,1.5,'8','lime','s','gold'])
calc_systems.append([ 2.413   ,  4.542  ,   4.840 , 0.0  ,r'$\gamma$-WO$_3$(002)-5f',tcolor3,0.0,0.08,1.5,'8','lime','D','gold'])
calc_systems.append([ 1.904   ,  4.427  ,   5.012 , 0.0  ,r'$\gamma$-WO$_3$(020)-5f',tcolor3,0.0,0.08,1.5,'8','lime','d','gold'])

calc_systems.append([ 0.099   ,  2.198  ,   3.710 , 0.0  ,r'HT-WO$_3$(110)',tcolor3,0.0,0.08,1.5,'8','lime','o','darkorange'])
calc_systems.append([-0.060   ,  2.142  ,   3.575 , 0.0  ,r'$\alpha$-WO$_3$(110)',tcolor3,0.0,0.08,1.5,'8','lime','s','darkorange'])

calc_systems.append([ 0.033   ,  1.418  ,   3.307 , 0.0  ,r'HT-WO$_3$(001)-6f-Ir',tcolor3,0.0,0.08,1.5,'8','lime','o','steelblue'])
calc_systems.append([ 0.007   ,  1.379  ,   3.288 , 0.0  ,r'$\alpha$-WO$_3$(002)-6f-Ir',tcolor3,0.0,0.08,1.5,'8','lime','s','steelblue'])
calc_systems.append([-0.059   ,  1.338  ,   3.237 , 0.0  ,r'$\gamma$-WO$_3$(002)-6f-Ir',tcolor3,0.0,0.08,1.5,'8','lime','D','steelblue'])
calc_systems.append([ 0.058   ,  1.503  ,   3.213 , 0.0  ,r'$\gamma$-WO$_3$(020)-6f-Ir',tcolor3,0.0,0.08,1.5,'8','lime','d','steelblue'])

calc_systems.append([ 1.058   ,  2.499  ,   4.084 , 0.0  ,r'HT-WO$_3$(001)-5f-Ir',tcolor3,0.0,0.08,1.5,'8','lime','o','skyblue'])
calc_systems.append([ 1.068   ,  2.515  ,   4.111 , 0.0  ,r'$\alpha$-WO$_3$(002)-5f-Ir',tcolor3,0.0,0.08,1.5,'8','lime','s','skyblue'])
calc_systems.append([ 1.058   ,  2.529  ,   4.068 , 0.0  ,r'$\gamma$-WO$_3$(002)-5f-Ir',tcolor3,0.0,0.08,1.5,'8','lime','D','skyblue'])
calc_systems.append([ 0.797   ,  2.341  ,   3.959 , 0.0  ,r'$\gamma$-WO$_3$(020)-5f-Ir',tcolor3,0.0,0.08,1.5,'8','lime','d','skyblue'])

calc_systems.append([ 0.482   ,  1.939  ,   3.725 , 0.0  ,r'HT-WO$_3$(110)-Ir',tcolor3,0.0,0.08,1.5,'8','lime','o','deepskyblue'])
calc_systems.append([ 0.438   ,  1.880  ,   3.755 , 0.0  ,r'$\alpha$-WO$_3$(110)-6f-Ir',tcolor3,0.0,0.08,1.5,'8','lime','s','deepskyblue'])
calc_systems.append([ 0.440   ,  1.956  ,   3.799 , 0.0  ,r'$\alpha$-WO$_3$(110)-5f-Ir',tcolor3,0.0,0.08,1.5,'8','lime','s','deepskyblue'])

calc_systems.append([ 1.045   ,  2.561  ,   4.066 , 0.0  ,r'W-doped $\alpha$-IrO$_3$(110)-Ir','limegreen',0.0,0.08,1.5,'8','lime','^','none'])
calc_systems.append([ 0.869   ,  2.350  ,   3.890 , 0.0  ,r'W-doped $\alpha$-IrO$_3$(110)-Ir','limegreen',0.0,0.08,1.5,'8','lime','^','none'])
calc_systems.append([ 1.501   ,  2.911  ,   4.805 , 0.0  ,r'W-doped $\alpha$-IrO$_3$(110)-W',tcolor4,0.0,0.08,1.5,'8','lime','^','lightgreen'])
calc_systems.append([-0.002   ,  1.365  ,   3.275 , 0.0  ,r'W-doped R-IrO$_2$(110)-Ir',tcolor4,0.0,0.08,1.5,'8','lime','v','hotpink'])
calc_systems.append([-0.502   ,  0.197  ,   3.252 , 0.0  ,r'W-doped R-IrO$_2$(110)-W',tcolor4,0.0,0.08,1.5,'8','lime','v','pink'])



'''
# solvation correction applied for NiO2 and LNO-NiO2 interfaces ----Jiang Li Feb-18,2020
#calc_systems.append([0.458+solv_corr_OH   ,  1.962+solv_corr_O    , 3.801+solv_corr_OOH   ,0.0  ,r'$LNO-NiO2-inter-1$','black',0.0,0.08,1.5,'8','lime',symbols[4],colors[5]])
calc_systems.append([0.484 +solv_corr_OH   ,  2.204+solv_corr_O    , 3.886+solv_corr_OOH   ,0.0  ,r'$LNO-NiO2-inter-2-O_S-1$','turquoise',0.0,0.08,1.5,'8','lime',symbols[4],colors[3]])
calc_systems.append([0.615 +solv_corr_OH   ,  2.204+solv_corr_O    , 3.977+solv_corr_OOH   ,0.0  ,r'$LNO-NiO2-inter-2-O_S-2$','turquoise',0.0,0.08,1.5,'8','lime',symbols[4],colors[4]])
calc_systems.append([0.450 +solv_corr_OH   ,  2.127+solv_corr_O    , 3.968+solv_corr_OOH   ,0.0  ,r'$LNO-NiO2-inter-2-O_S-4$','turquoise',0.0,0.08,1.5,'8','lime',symbols[4],colors[5]])
calc_systems.append([0.527 +solv_corr_OH   ,  1.648+solv_corr_O    , 3.823+solv_corr_OOH   ,0.0  ,r'$NiO2-001$','turquoise',0.0,0.08,1.5,'8','lime',symbols[4],'y'])
calc_systems.append([0.484   ,  2.204  ,   3.886 ,0.0  ,r'$LNO-NiO2-inter-2-O_S-1$','w',0.0,0.08,1.5,'8','lime',symbols[4],colors[3]])
calc_systems.append([0.615   ,  2.204  ,   3.977 ,0.0  ,r'$LNO-NiO2-inter-2-O_S-2$','w',0.0,0.08,1.5,'8','lime',symbols[4],colors[4]])
calc_systems.append([0.450   ,  2.127  ,   3.968 ,0.0  ,r'$LNO-NiO2-inter-2-O_S-4$','w',0.0,0.08,1.5,'8','lime',symbols[4],colors[5]])
calc_systems.append([0.527   ,  1.648  ,   3.823 ,0.0  ,r'$NiO2-001$',tcolor,0.0,0.08,1.5,'8','lime',symbols[4],'y'])
'''
import csv
with open('Final_dGs_for_MOOHx_system_OER.csv', 'w') as myfile:
    fieldnames = ['Surface name','dOH','dO','dOOH', 'overpotential', 'onset potential', 'PLS']
    wr = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL, delimiter=',')
    writer=csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()

    print ('Final dGs for MOOHx_system  system')

    #solv_corr=0.0
    for i in range(len(calc_systems)):
        calc_systems[i][0]+=0 #-0.188  # solvation correction for ad_OH
        calc_systems[i][1]+=0 #-0.003  # solvation correction for ad_O
        calc_systems[i][2]+=0 #-0.057  # solvation correction for ad_OOH
        recalculated_over=overpotential_oer_full(calc_systems[i][0],calc_systems[i][1],calc_systems[i][2])
        #print  calc_systems[i][4], calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], recalculated_over[0],recalculated_over[1], ' limiting step is ', recalculated_over[2]
        calc_systems[i][3]=recalculated_over[0]
        wr.writerow([calc_systems[i][4], calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], recalculated_over[0], recalculated_over[1], recalculated_over[2]])

import csv
with open('Final_dEs_for_MOOHx_system_OeR.csv', 'w') as myfile:
    fieldnames = ['Surface name','dOH','dO','dOOH', 'overpotential', 'onset potential', 'PLS']
    wr = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL, delimiter=',')
    writer=csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()

    print ('Final dEs for MOOHx_system  system')

    for i in range(len(calc_systems)):
        recalculated_over=overpotential_oer_full(calc_systems[i][0],calc_systems[i][1],calc_systems[i][2])
        #print  calc_systems[i][4], calc_systems[i][0]-0.30225 , calc_systems[i][1]-(-0.0145), calc_systems[i][2]-0.34475, recalculated_over[0],recalculated_over[1], ' limiting step is ', recalculated_over[2]
        calc_systems[i][3]=recalculated_over[0]
        wr.writerow([calc_systems[i][4],  calc_systems[i][0]-0.30225  , calc_systems[i][1]-(-0.0145), calc_systems[i][2]-0.34475 , recalculated_over[0], recalculated_over[1], recalculated_over[2]])
#foo=r': %f' % (calc_systems[i][3]) 




for i in range(len(calc_systems)):
   #ax.plot(calc_systems[i][1]-calc_systems[i][0], calc_systems[i][0],'or',color=calc_systems[i][5])
   ax.plot(calc_systems[i][1]-calc_systems[i][0], calc_systems[i][0],calc_systems[i][9],mec=calc_systems[i][5], mfc=calc_systems[i][12],mew=0.8,zorder=4,marker=calc_systems[i][11],label=calc_systems[i][4]+' : %.2f V' % (calc_systems[i][3]))
   #if i!=0 and 1:
#	ax.text(calc_systems[i][1]-calc_systems[i][0]+calc_systems[i][6], calc_systems[i][0]+calc_systems[i][7],calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]), color='black',fontsize=6,horizontalalignment='center',rotation=0,zorder=1)
 #  else:
  #      ax.text(calc_systems[i][1]-calc_systems[i][0]+calc_systems[i][6], calc_systems[i][0]+calc_systems[i][7],calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]), color='white',fontsize=6,horizontalalignment='center',rotation=0,zorder=1)
   #ax.text(calc_systems[i][0],calc_systems[i][1],'%i' %(i+1), color='black',fontsize=4,
   #        horizontalalignment='center',
   #        verticalalignment='center',
   #        rotation=0,zorder=2)

corners=[[1.3,1.0],[x1+(x2-x2)*0.2,y1+(y2-y1)*0.9],[x1+(x2-x2)*0.8,y1+(y2-y1)*0.1],[-2,0]]

#for i in range(len(corners)):
#   ax.text(corners[i][0],corners[i][1], overpotential_label(corners[i][0],corners[i][1]), color='white',fontsize='x-small',horizontalalignment='center',rotation=0,zorder=3) 

ax.legend(bbox_to_anchor=(-0.15, 1.7), loc=2, borderaxespad=1, ncol=3, fancybox=True, shadow=False, fontsize='x-small',handlelength=2)


fig.savefig('OER_contour_v1.pdf', bbox_inches='tight')
# clear the plot
fig.clf()

fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
x1=-1
x2=2.5
ax.axis([x1, x2, x1, ooh_oh_scaling(x2)])


ax.set_xlabel(r'$\Delta$G$_{\sf OH}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$,$\Delta$G$_{\sf O}$ (eV)')

xdata=[]
ydata=[]
y2data=[]

for i in range(len(calc_systems)):
#for i in range(3):
   xdata.append(calc_systems[i][0])
   ydata.append(calc_systems[i][2])
   y2data.append(calc_systems[i][1])


#print xdata
#print ydata

from pylab import *
fit = polyfit(xdata,ydata,1)
fit_fn = poly1d(fit)
print (fit_fn)
aa=fit_fn[1]
bb=fit_fn[0]

fit1 = polyfit(xdata,y2data,1)
fit_fn1 = poly1d(fit1)
print (fit_fn1)
#print fit_fn[0], fit_fn[1]
#how bad is scaling
for i in range(len(calc_systems)):
        error=calc_systems[i][2]-(fit_fn[1]*calc_systems[i][0]+fit_fn[0])
        #print error, calc_systems[i]

xx = np.arange(x1,x2, delta)
ax.plot(xx,fit_fn[1]*xx+fit_fn[0],'--',lw=1, dashes=(3,1),c='black',label='OOH scaling')
ax.plot(xx,xx+3.2,'--',lw=1, dashes=(3,1),c='grey')
ax.plot(xx,xx,'--',lw=1, dashes=(3,1),c='grey')
ax.plot(xx,fit_fn1[1]*xx+fit_fn1[0],'--',lw=1, dashes=(3,1),c='orangered',label='O scaling')


for i in range(len(calc_systems)):
	ax.plot(calc_systems[i][0],calc_systems[i][2],'ro',
           ms=3,marker=calc_systems[i][11],
           #alpha=0.2,
           color='black')
	ax.plot(calc_systems[i][0],calc_systems[i][1],'ro',
           ms=3,marker=calc_systems[i][11],
           #alpha=0.2,
           color='orangered')
	ax.plot(calc_systems[i][0],calc_systems[i][0],calc_systems[i][9],mec=calc_systems[i][5], mfc=calc_systems[i][12],mew=0.8,zorder=4,marker=calc_systems[i][11],label=calc_systems[i][4]+' : %.2f V' % (calc_systems[i][3]))
#   	ax.text(calc_systems[i][0], calc_systems[i][0]+calc_systems[i][7]+0.08,
#           calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]), color='black',fontsize=6,horizontalalignment='center',
 #          rotation=0,zorder=1)

#ax.legend(bbox_to_anchor=(-0.15, 1.55), loc=2, borderaxespad=0.5, ncol=3, fancybox=True, shadow=True, fontsize='x-small',handlelength=2)
eq1=r'$\Delta$G$_{\sf OOH}$=%.2f$\Delta$G$_{\sf OH}$+%.2f' %(fit_fn[1],fit_fn[0])
ax.text(-0.8,5.1,eq1,fontsize='small',color='black')
eq2=r'$\Delta$G$_{\sf OOH}$=$\Delta$G$_{\sf OH}$+3.2'
ax.text(-0.8,4.7,eq2,fontsize='small',color='grey')
eq3=r'$\Delta$G$_{\sf O}$=%.2f$\Delta$G$_{\sf OH}$+%.2f' %(fit_fn1[1],fit_fn1[0])
ax.text(-0.8,4.3,eq3,fontsize='small',color='orangered')

fig.savefig('OER_scaling_v1.pdf', bbox_inches='tight')
# clear the plot
fig.clf()
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

#x1=1.23-1
#x2=1.23+1
#y2=1
#y1=0

x1=0.0
x2=2.8
y2=3.43
y1=1.23

ax.axis([x1, x2, y1, y2])
delta = 0.01
x = np.arange(x1,x2, delta)

ax.set_xlabel(r'$\Delta$G$_{\sf O}-\Delta$G$_{\sf OH}$ (eV)')
#ax.set_ylabel(r'$\Delta$G$_{\sf O}$ (eV)')

ax.set_ylabel(r'U$_{\sf OER}$ (V)')

ax.set_ylim(ax.get_ylim()[::-1])


plot(x,np.maximum(x,3.2-x),'--',color='black',lw=0.67, dashes=(3,1),zorder=2)
#plot(x,1.23,'--',color='black',lw=0.67, dashes=(3,1),zorder=2)


#xy=np.array([xp for xp in x if 1.55<xp<1.66])
#ax.fill_between(xy, y2, np.maximum(xy,3.2-xy)-1.23, zorder=1, color='red', alpha=0.3, edgecolor="none")

for b in x:
	if(np.maximum(b,3.2-b)<0.44):
		print (b)

#plot(x,np.maximum(x,bb-x*(aa)-0.65)-1.23,'--',color='grey',lw=0.67, dashes=(3,1))
#slope not one
plot(x,np.maximum(x,3.18-0.82*x-0.35)-1.23,'--',color='pink',lw=0.67,dashes=(3,1))
#plot(x,np.maximum(x,2.46-x)-1.23,'-',color='black',lw=0.67)

#import matplotlib.patches as patches
#ax.add_patch(
#    patches.Rectangle(
#        (calc_systems[1][1]-calc_systems[1][0],calc_systems[1][3]-0.04),   # (x,y)
#        0.25,          # width
#        0.05,          # height
#        fill=True,
#        edgecolor="none", 
#        facecolor='red',
#    )
#)


for i in range(len(calc_systems)):
   ax.plot(calc_systems[i][1]-calc_systems[i][0],calc_systems[i][3]+1.23,calc_systems[i][9],mec=calc_systems[i][5], mfc=calc_systems[i][12],mew=0.8,zorder=4,marker=calc_systems[i][11],label=calc_systems[i][4]+' : %.2f V' %(calc_systems[i][3]+1.23))
  # if(i!=1):
#   ax.text(calc_systems[i][1]-calc_systems[i][0],calc_systems[i][3]-0.02,calc_systems[i][3])
#		color='black',fontsize=6,horizontalalignment='left',rotation=0,zorder=4)
 #  else:
  # 	ax.text(calc_systems[i][1]-calc_systems[i][0],calc_systems[i][3]-0.02,calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]),
#		color='black',fontsize=6,horizontalalignment='right',rotation=0,zorder=4)
ax.legend(bbox_to_anchor=(-0.15, 1.7), loc=2, borderaxespad=0.5, ncol=3, fancybox=True, shadow=False, fontsize='x-small',handlelength=2)
fig.savefig('OER_1D_v1.pdf', bbox_inches='tight')
#fig.show('OER_1D_v1.pdf')

for i in range(len(calc_systems)):
	print ('%s	Edge	OH	%f	[]	ZhenghangZhao_calculated' %(calc_systems[i][4],calc_systems[i][0]))
	print ('%s	Edge	O	%f	[]	ZhenghangZhao_calculated' %(calc_systems[i][4],calc_systems[i][1]))
	print ('%s	Edge	OOH	%f	[]	ZhenghangZhao_calculated' %(calc_systems[i][4],calc_systems[i][2]))
for i in range(len(calc_systems)):	
	print ("'%s'," %(calc_systems[i][4]))



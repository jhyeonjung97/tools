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

zoom=0.3
d1=3*zoom
d2=4*zoom
xcenter=0.5 #0.65
#ycenter=1.23#2.4
ycenter=-1.0#2.4

x1=xcenter-d1#-0.6
x2=xcenter+d1#2.2
y1=ycenter-d2#1#0.5
y2=ycenter+d2#5
ax.axis([x1, x2, y1, y2])
ax.set_xlabel(r'$\Delta$G$_{\sf COH*}$ - $\Delta$G$_{\sf CO*}$ (eV)')
#ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$ -$\Delta$G$_{\sf O}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf CO*}$ (eV)')




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
def dgc_dgco_scaling(dgco):
    #dgc=3.5*dgco + 2.18 
    dgc=2.2*dgco + 2.34 
    return dgc

corr_onset=0.26 # to CH4 Nitoppi 2019 
ch4=-2.16

#to make the countour 
def overpotential_from_dgc_dgco_scaling(x,dgco): # x=dgcoh-dgco 
    dgc_from_scaling=dgc_dgco_scaling(dgco)
    dg14=[dgco, x, dgc_from_scaling -(x+dgco), (ch4-dgc_from_scaling)/4.0 ]
    #dg14=[dgco, x, dgc_from_scaling -(x+dgco) ]
    m=max(dg14)
    return m+corr_onset
    #return doh*do


def overpotential_corr_full(dgco,dgcoh,dgc,dgch):
    #do=1.556*doh + 0.9951 #from our scaling
    dg23=[dgco,dgcoh-dgco,dgc-dgcoh, dgch-dgc, (ch4-dgch)/3.0] #add solvation soon
    m=max(dg23)
    ##print dg14.index(m)
    #print dg14,m
    print (dg23)
    #           overpotentail, onset potential (is negative)
    return [round(m+corr_onset,2),-round(m,2),corr_step(dg23.index(m))]


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

steps=['CO(g)+*->CO*','CO*->COH*', 'COH*->C*', 'C*->CH*', 'CH*->CH4(g)']
def corr_step (i):
        return str(steps[i]) 



#Z=overpotential(X,Y)

Z=[]
for j in y:
    tmp=[]
    for i in x:
        tmp.append((overpotential_from_dgc_dgco_scaling(i,j)))
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
levels = np.arange(0.0, 1.5 ,0.1)
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


cbar = plt.colorbar(CS, ticks=list(np.arange(0.0, 1.4, step=0.1)))
#cbar = plt.colorbar(CS)

#cbar.ax.set_ylabel('Overpotential [V]')
#cbar.ax.set_ylabel(r'$\eta_{\sf calc.}$')
cbar.ax.set_ylabel(r'$\eta_{\sf CORR}$ (V)')

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

ddgco=0.09+0.45  # COg -> CO*
ddgcoh=0.37+0.485 # COg+Hg -> COH* 
ddgc=0.08
ddgch=0.33  

co_solv=-0.06
#coh_solv=-0.4 #from Carter
coh_solv=-0.1 #from Carter
coh_solv_bonus=-0.3 #from Carter, only for weak COH
c_solv=0.0
ch_solv=0.0

hc_dest=0.65 #high coverage of CO* penatly (2014 PCCP)

calc_systems.append([-1.76+hc_dest ,-1.47+hc_dest ,-1.71 ,-1.94 , 0.0  ,r'Rh(100) (0.66 ML CO*)',tcolor3,0.0,0.08,1.5,'8','grey','o','purple'])
calc_systems.append([-1.75+hc_dest ,-1.20+hc_dest ,-1.26 , -1.26,  0.0  ,r'Rh(211) (0.66 ML CO*)',tcolor3,0.0,0.08,1.5,'8','grey','D','purple'])
                                  #-1.32               

calc_systems.append([-1.63+hc_dest ,-0.94+hc_dest,-1.28,-1.27 , 0.0  ,r'Pd(100) (0.66 ML CO*)',tcolor3,0.0,0.08,1.5,'8','brown','o','brown'])
calc_systems.append([-1.65+hc_dest ,-0.96+hc_dest,-1.16, -1.16 , 0.0  ,r'Pd(211) (0.66 ML CO*)',tcolor3,0.0,0.08,1.5,'8','brown','D','brown'])
                                   #-0.86
                                                 
calc_systems.append([-1.65+hc_dest ,-1.08+hc_dest,-1.06,-1.54, 0.0 ,r'Pt(100) (0.66 ML CO*)',tcolor3,0.0,0.08,1.5,'8','green','o','darkblue'])
calc_systems.append([-1.71+hc_dest ,-1.14+hc_dest,-0.52, -0.52 , 0.0  ,r'Pt(211) (0.66 ML CO*)',tcolor3,0.0,0.08,1.5,'8','green','D','darkblue'])
                                   #-1.27
                                                 
calc_systems.append([-1.64+hc_dest ,-1.34+hc_dest,-1.66,-2.24 , 0.0  ,r'Ni(100) (0.66 ML CO*)',tcolor3,0.0,0.08,1.5,'8','green','o','green'])
calc_systems.append([-1.70+hc_dest ,-1.04+hc_dest,-1.30, -1.30 , 0.0  ,r'Ni(210) (0.66 ML CO*)',tcolor3,0.0,0.08,1.5,'8','green','D','green'])
                                   #-0.84

calc_systems.append([-0.65, 0.34+coh_solv_bonus, 0.38,-0.39, 0.0  ,r'Cu(100) (Beef, ++solv.)',tcolor3,0.0,0.08,1.5,'8','lime','o','red'])

calc_systems.append([-0.65, 0.3+0.4, 0.38,-0.39, 0.0  ,r'Cu(100) (Beef, ++solv.)',tcolor3,0.0,0.08,1.5,'8','lime','o','orange'])

calc_systems.append([-0.49, 0.62+coh_solv_bonus, 1.51, 0.14,0.0  ,r'Cu(111) (++solv.)',tcolor3,0.0,0.08,1.5,'8','lime','^','red'])
calc_systems.append([-0.6 , 1.10+coh_solv_bonus, 1.51, 0.14,0.0  ,r'Cu(111) (DMC, ++solv.)',tcolor3,0.0,0.08,1.5,'8','lime','d','orangered'])
calc_systems.append([-1.21+hc_dest, 0.86+coh_solv_bonus, 0.90,-0.04,0.0  ,r'Cu(211) (0.66 ML CO ++solv.)',tcolor3,0.0,0.08,1.5,'8','lime','D','red'])

#calc_systems.append([-0.05+ddgco,1.37+ddgcoh,2.35+ddgc, 0.0  ,r'Au(111)',tcolor3,0.0,0.08,1.5,'8','gold','^','gold'])
calc_systems.append([-0.05,1.37+coh_solv_bonus,2.35,0.99, 0.0  ,r'Au(111) (++solv.)',tcolor3,0.0,0.08,1.5,'8','gold','^','gold'])



import csv
with open('Final_dGs_for_MOOHx_system_CORR.csv', 'w') as myfile:
    fieldnames = ['Surface name','dGCO','dGCOH','dGC','dGCH', 'overpotential', 'onset potential', 'PLS']
    wr = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL, delimiter=',')
    writer=csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()

    print ('Final dGs for MOOHx_system  system')

    solv_corr= [co_solv,coh_solv,c_solv,ch_solv]
    ddg_corr=[ddgco, ddgcoh, ddgc, ddgch]
    for i in range(len(calc_systems)):
        for k in range(len(steps)-1):
           calc_systems[i][k]+=solv_corr[k]+ddg_corr[k]  # solvation correction for ad_OH

        recalculated_over=overpotential_corr_full(calc_systems[i][0],calc_systems[i][1],calc_systems[i][2],calc_systems[i][3])
        #print  calc_systems[i][4], calc_systems[i][0], calc_systems[i][1], calc_systems[i][2], recalculated_over[0],recalculated_over[1], ' limiting step is ', recalculated_over[2]
        calc_systems[i][4]=recalculated_over[0]
        wr.writerow([calc_systems[i][5], calc_systems[i][0], calc_systems[i][1], calc_systems[i][2],calc_systems[i][3], recalculated_over[0], recalculated_over[1], recalculated_over[2]])

import csv
with open('Final_dEs_for_MOOHx_system_CORR.csv', 'w') as myfile:
    fieldnames = ['Surface name','dECO','dECOH','dEC', 'overpotential', 'onset potential', 'PLS']
    wr = csv.writer(myfile, quoting=csv.QUOTE_MINIMAL, delimiter=',')
    writer=csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()

    print ('Final dEs for MOOHx_system  system')

    for i in range(len(calc_systems)):
        recalculated_over=overpotential_corr_full(calc_systems[i][0],calc_systems[i][1],calc_systems[i][2],calc_systems[i][3])
        #print  calc_systems[i][4], calc_systems[i][0]-0.30225 , calc_systems[i][1]-(-0.0145), calc_systems[i][2]-0.34475, recalculated_over[0],recalculated_over[1], ' limiting step is ', recalculated_over[2]
        calc_systems[i][4]=recalculated_over[0]
        wr.writerow([calc_systems[i][5],  calc_systems[i][0]-ddgco  , calc_systems[i][1]-ddgcoh, calc_systems[i][2]-ddgc ,calc_systems[i][3]-ddgch, recalculated_over[0], recalculated_over[1], recalculated_over[2]])
#foo=r': %f' % (calc_systems[i][3]) 




for i in range(len(calc_systems)):
   #ax.plot(calc_systems[i][1]-calc_systems[i][0], calc_systems[i][0],'or',color=calc_systems[i][5])
   ax.plot(calc_systems[i][1]-calc_systems[i][0], calc_systems[i][0],calc_systems[i][10],mec=calc_systems[i][6], mfc=calc_systems[i][13],mew=0.8,zorder=4,marker=calc_systems[i][12],label=calc_systems[i][5]+' : %.2f V' % (calc_systems[i][4]))
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


fig.savefig('CORR_contour_v1.pdf', bbox_inches='tight')
# clear the plot
fig.clf()

fig = plt.figure(figsize=fig_size,dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
x1=-1.4
x2=0.1
ax.axis([x1, x2, x1, dgc_dgco_scaling(x2)])


ax.set_xlabel(r'$\Delta$G$_{\sf CO}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf COH}$,$\Delta$G$_{\sf C}$ (eV)')

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
ax.plot(xx,fit_fn[1]*xx+fit_fn[0],'--',lw=1, dashes=(3,1),c='black',label='C* scaling')
#ax.plot(xx,3.5*xx+2.18,'--',lw=1, dashes=(3,1),c='grey')
ax.plot(xx,xx,'--',lw=1, dashes=(3,1),c='grey')
ax.plot(xx,fit_fn1[1]*xx+fit_fn1[0],'--',lw=1, dashes=(3,1),c='orangered',label='COH* scaling')


for i in range(len(calc_systems)):
	ax.plot(calc_systems[i][0],calc_systems[i][2],'ro',
           ms=3,marker=calc_systems[i][12],
           #alpha=0.2,
           color='black')
	ax.plot(calc_systems[i][0],calc_systems[i][1],'ro',
           ms=3,marker=calc_systems[i][12],
           #alpha=0.2,
           color='orangered')
	ax.plot(calc_systems[i][0],calc_systems[i][3],'ro',
           ms=3,marker=calc_systems[i][12],
           #alpha=0.2,
           color='gold')
	ax.plot(calc_systems[i][0],calc_systems[i][0],calc_systems[i][10],mec=calc_systems[i][6], mfc=calc_systems[i][13],mew=0.8,zorder=4,marker=calc_systems[i][12],label=calc_systems[i][5]+' : %.2f V' %(-(calc_systems[i][4]+corr_onset)))

ax.legend(bbox_to_anchor=(-0.15, 1.55), loc=2, borderaxespad=0.5, ncol=3, fancybox=True, shadow=True, fontsize='x-small',handlelength=2)
eq1=r'$\Delta$G$_{\sf C}$=%.2f$\Delta$G$_{\sf CO}$+%.2f' %(fit_fn[1],fit_fn[0])
ax.text(-1.1,2.5,eq1,fontsize='small',color='black')
eq3=r'$\Delta$G$_{\sf COH}$=%.2f$\Delta$G$_{\sf CO}$+%.2f' %(fit_fn1[1],fit_fn1[0])
ax.text(-1.1,2.2,eq3,fontsize='small',color='orangered')
eq2=r'$\Delta$G$_{\sf CH}$=xx$\Delta$G$_{\sf CO}$+xx'
ax.text(-1.1,1.9,eq2,fontsize='small',color='gold')

fig.savefig('CORR_scaling_v1.pdf', bbox_inches='tight')
# clear the plot
fig.clf()
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])

#x1=1.23-1
#x2=1.23+1
#y2=1
#y1=0

x1=-1.0
x2=2.0
y2=-2.0
y1=-0.0

ax.axis([x1, x2, y1, y2])
delta = 0.01
x = np.arange(x1,x2, delta)

ax.set_xlabel(r'$\Delta$G$_{\sf COH}-\Delta$G$_{\sf CO}$ (eV)')
#ax.set_ylabel(r'$\Delta$G$_{\sf O}$ (eV)')

ax.set_ylabel(r'U$_{\sf CORR}$ (V)')

ax.set_ylim(ax.get_ylim()[::-1])


co_fixed=-0.5
plot(x,-np.maximum(x, -1.5 +2.34 -x),'--',color='black',lw=0.67, dashes=(3,1),zorder=2)
plot(x,corr_onset+0.0*x,'--',color='black',lw=0.67, dashes=(3,1),zorder=2)
#plot(x,1.23,'--',color='black',lw=0.67, dashes=(3,1),zorder=2)


#xy=np.array([xp for xp in x if 1.55<xp<1.66])
#ax.fill_between(xy, y2, np.maximum(xy,3.2-xy)-1.23, zorder=1, color='red', alpha=0.3, edgecolor="none")

for b in x:
	if(np.maximum(b,3.2-b)<0.44):
		print (b)


for i in range(len(calc_systems)):
   ax.plot(calc_systems[i][1]-calc_systems[i][0],-(calc_systems[i][4]-corr_onset),calc_systems[i][10],mec=calc_systems[i][6], mfc=calc_systems[i][13],mew=0.8,zorder=4,marker=calc_systems[i][12],label=calc_systems[i][5]+' : %.2f V' %(-(calc_systems[i][4]-corr_onset)))
  # if(i!=1):
#   ax.text(calc_systems[i][1]-calc_systems[i][0],calc_systems[i][3]-0.02,calc_systems[i][3])
#		color='black',fontsize=6,horizontalalignment='left',rotation=0,zorder=4)
 #  else:
  # 	ax.text(calc_systems[i][1]-calc_systems[i][0],calc_systems[i][3]-0.02,calc_systems[i][4]+'(%.2f)' %(calc_systems[i][3]),
#		color='black',fontsize=6,horizontalalignment='right',rotation=0,zorder=4)
ax.legend(bbox_to_anchor=(-0.15, 1.7), loc=2, borderaxespad=0.5, ncol=3, fancybox=True, shadow=False, fontsize='x-small',handlelength=2)
fig.savefig('CORR_1D_v1.pdf', bbox_inches='tight')
#fig.show('OER_1D_v1.pdf')



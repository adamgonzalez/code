"""
Author: Adam Gonzalez
Descr.: Use result from Gonzalez+2017 to determine coronal parameters.
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import time
import os
from matplotlib import rcParams
rcParams['xtick.top'] = 'True'
rcParams['ytick.right'] = 'True'
rcParams['xtick.direction'] = 'in'#'inout'
rcParams['ytick.direction'] = 'in'#'inout'
rcParams['xtick.major.size'] = 3.5#7
rcParams['xtick.major.width'] = 1.2
rcParams['ytick.major.size'] = 3.5#7
rcParams['ytick.major.width'] = 1.2
rcParams['xtick.minor.size'] = 2#4
rcParams['ytick.minor.size'] = 2#4
rcParams['axes.linewidth'] = 1.5
rcParams['font.family'] = 'serif'
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = ['\\boldmath']

#matplotlib.rcParams.update({'font.size': 18})
#matplotlib.rcParams['axes.linewidth'] = 1 #set the value globally
#plt.rc('font',family='serif')
#plt.rc('text',usetex=True)

####################################################################################################
# Compute the approximation for R given z and beta
def R_calc(h,v):
    a = 0.986

    mu_in = (2.0 - h**2.0)/(h**2.0 + a**2.0)
    mu_out = (2.0*h)/(h**2.0 + a**2.0)
    mu_in_p = (mu_in - v)/(1.0 - v*mu_in)
    mu_out_p = (mu_out - v)/(1.0 - v*mu_out)

    value = (mu_out_p - mu_in_p)/(1.0 - mu_out_p)

    return value
####################################################################################################

os.chdir("/Users/agonzalez/Documents/Research/Data/Mrk335/lowhi")

minh, maxh = 1.0, 10.0
minv, maxv = 0.0, 1.0

for k in range (0, 1):
    RUN = k+1
    print "RUN: ", RUN

    print "Setting up the uniform grid..."
    t0 = time.time()
    z = np.linspace(minh, maxh, num=(maxh-minh)/9e-3)
    beta = np.linspace(minv, maxv, num=(maxv-minv)/1e-3)
    t1 = time.time()
    print t1-t0, "\n"
    resh, resv = len(z), len(beta)

    Rvs = np.zeros([resh+1,resv+1])

    print "Computing R..."
    t0 = time.time()
    # compute R as function of source height and source velocity
    for i in range (0, resh):
        for j in range (0, resv):
            Rvs[0,j+1] = beta[j]
            Rvs[i+1,j+1] = R_calc(z[i],beta[j])
        Rvs[i+1,0] = z[i]
    t1 = time.time()
    print t1-t0, "\n"

    # plotting up the escape velocity curves
    plt.figure()
    ax = plt.subplot(111)
    ##---------------------------------------------------------------------------------------
    # Compute the escape velocity for a black hole of mass M at a height R above the black hole
    def vesc_calc(G,M,R,c):
        v = np.sqrt((2.0*G*M)/R)/c

        return v

    G = 6.674e-11
    c = 2.998e8
    M_sun = 1.989e30

    # plt.figure()
    # ax = plt.subplot(111)
    col = ['k','k','k']

    res = 100 
    Vesc = np.zeros([5,res])
    R = np.zeros([5,res])

    for j in range (0,3):
        ### I Zw 1
        #if (j==0):
        #    M_bh = pow(10.0, 7.30)*M_sun ; name = 'Negrete et al. (2012)'
        #    r_g0 = (G*M_bh)/(c**2.0)
        #if (j==1):
        #    M_bh = pow(10.0, 7.30+0.23)*M_sun ; name = 'Mass + error'
        #if (j==2):
        #    M_bh = pow(10.0, 7.30-0.19)*M_sun ; name = 'Mass -- error'

        ### III Zw 2
        # if (j==0):
        #     M_bh = 184000000.*M_sun ; name = 'van den Bosch (2016) <- Grier et al. (2012)'
        #     r_g0 = (G*M_bh)/(c**2.0)
        # if (j==1):
        #     M_bh = (184000000.+27000000.)*M_sun ; name = '+'
        # if (j==2):
        #     M_bh = (184000000.-27000000.)*M_sun ; name = '--'

        ### Mrk 335
        if (j==0):
            M_bh = (2.5e7)*M_sun ; name = 'Grier et al. (2012)'
            r_g0 = (G*M_bh)/(c**2.0)
        if (j==1):
            M_bh = (2.5e7 + 3e6)*M_sun ; name = 'Mass + error'
        if (j==2):
            M_bh = (2.5e7 - 3e6)*M_sun ; name = 'Mass -- error'

        R_s = (2.0*G*M_bh)/(c**2.0)
        r_g = (G*M_bh)/(c**2.0)

        R[j][:] = np.logspace(start=np.log10(1.01*R_s), stop=np.log10(1000.0*r_g), num=res)

        for i in range (0,res):
            Vesc[j][i] = vesc_calc(G,M_bh,R[j][i],c)

        # print "Mass of I Zw 1 BH [kg]   = ", M_bh
        # print "Schwarzschild radius [m] = ", R_s
        # print "Gravitationl radius [m]  = ", r_g
        R[j][:] = R[j][:]/r_g0

        if (j!=0):
            ax.plot(R[j][:],Vesc[j][:], color=col[j], dashes=[5,3], alpha=0.75, label=name)
        elif (j==0):
            ax.plot(R[j][:],Vesc[j][:], color=col[j], alpha=1.0, label=name)

    for i in range (0,res):
        R[3][i] = abs(R[0][i]-R[1][i])
        R[4][i] = abs(R[0][i]-R[2][i])

    # ax.fill_betweenx(y=Vesc[0][:], x1=R[0][:]-R[4][:], x2=R[0][:]+R[3][:], facecolor='red', alpha=0.05)
    ax.tick_params(axis='both', which='both', direction='in', top='on', right='on')
    ##---------------------------------------------------------------------------------------


    # Compute and plot the pairs (z,b) that match the reflection fraction desired
    c = 0
    pairs = [[0,0,0]]
    # minR, maxR = 0.54-0.04, 0.54+0.04       ### IZw1
    # minR, maxR = 0.204-0.033, 0.204+0.017   ### IIIZw2
    # minR, maxR = 0.23-0.03, 0.23+0.02   ### IIIZw2 XMM
    # minR, maxR = 0.12-0.02, 0.12+0.01   ### IIIZw2 Suzaku
    minR, maxR = 3.2-0.6, 3.2+0.6 	### Mrk 335
    print "Finding the pairs..."
    t0 = time.time()
    for i in range (0, resh):
        for j in range (0, resv):
            if (Rvs[i+1,j+1]<=maxR) and (Rvs[i+1,j+1]>=minR):
                c += 1
                pairs = np.append(pairs,[[Rvs[i+1,0],Rvs[0,j+1],Rvs[i+1,j+1]]], axis=0)
    t1 = time.time()

    print t1-t0, "\n"
    print 'Number of sources within R = ', minR, ' to ', maxR, ' is ', c
    print ''

    cfset = plt.scatter(pairs[1:,0], pairs[1:,1], c=pairs[1:,2], s=5.0, cmap='coolwarm', vmin=minR, vmax=maxR, alpha=1.0)
    cbar = plt.colorbar(cfset, pad=0.05)#, ticks=[-0.02, 0.0, 0.02, 0.04, 0.06, 0.08, 0.10])
    cbar.ax.set_ylabel('\\textbf{Reflection Fraction, $R$}', rotation='270', labelpad=15.0)
    plt.xlabel('\\textbf{Source Height} $[r_g]$')
    plt.ylabel('\\textbf{Source Velocity} $[c]$')
    plt.xlim(minh,maxh)
    plt.ylim(minv,maxv)

    plt.savefig('Mrk335_RZB.png', bbox_inches='tight', dpi=600)
    # plt.show()

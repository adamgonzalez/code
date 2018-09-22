'''
Author: Adam Gonzalez
Descr.: This is a script to plot the output of the MC simulations testing UFO feature significance.
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import stats
from matplotlib import rcParams
rcParams['xtick.top'] = 'True'
rcParams['ytick.right'] = 'True'
rcParams['xtick.direction'] = 'inout'
rcParams['ytick.direction'] = 'inout'
rcParams['xtick.major.size'] = '7'
rcParams['ytick.major.size'] = '7'
rcParams['xtick.minor.size'] = '4'
rcParams['ytick.minor.size'] = '4'
# rcParams['axes.facecolor'] = 'k'


### Source Redshift
# srcname = 'Mrk279' ; srcz = 0.0304051
# srcname = 'PG1211+143' ; srcz = 0.0809
# srcname = 'Mrk79' ; srcz = 0.022189
# srcname = 'Mrk841' ; srcz = 0.0336422

srcname = 'IZw1' ; srcz = 0.0589

Nthreads = 45


######### FAKE ##################################
### Read in the fake spectra results
fake_data = np.zeros((0,7))
for i in range(0, Nthreads+1):
    tmp = np.genfromtxt("results/mos1/fake_10eVresult_T{0:1d}.dat".format(i))
    fake_data = np.concatenate((fake_data,tmp), axis=0)

### Break down the data into named arrays
fake_spec = fake_data[:,0]
fake_E = fake_data[:,1]
fake_N = fake_data[:,2]
fake_modfit = fake_data[:,3]
fake_basefit = fake_data[:,5]
fake_deltaC = fake_modfit-fake_basefit
nf = len(fake_spec)

### Determine the energies sampled and number of iterations
ebins = np.arange(min(fake_E), max(fake_E)+0.1, 0.1)
trials = int(nf/len(ebins))
######### FAKE ##################################


######### REAL ##################################
### Using the steppar normalisation grid
Nsamples = 62
real_data = np.zeros((0,Nsamples+2))

if (Nthreads > 35):
    for i in range(36, Nthreads+1):
        if (os.stat("results/mos1/real_10eVresult_T{0:1d}.dat".format(i)).st_size < 1):
            print("Ignoring T{0:1d} as file is empty!".format(i))
            continue
        if (os.stat("results/mos1/real_10eVresult_T{0:1d}.dat".format(i)).st_size < 1000):
            tmp = np.genfromtxt("results/mos1/real_10eVresult_T{0:1d}.dat".format(i))
            tmp = np.reshape(tmp, (1,Nsamples+2))
        else:
            tmp = np.genfromtxt("results/mos1/real_10eVresult_T{0:1d}.dat".format(i))
        # tmp = np.genfromtxt("results/mos1/real_10eVresult_T{0:1d}.dat".format(i))
        real_data = np.concatenate((real_data,tmp), axis=0)
    Nthreads = 35

for i in range(0, Nthreads+1):
    if (os.stat("results/mos1/real_10eVresult_T{0:1d}.dat".format(i)).st_size < 1):
        print("Ignoring T{0:1d} as file is empty!".format(i))
        continue
    if (os.stat("results/mos1/real_10eVresult_T{0:1d}.dat".format(i)).st_size < 1000):
        tmp = np.genfromtxt("results/mos1/real_10eVresult_T{0:1d}.dat".format(i))
        tmp = np.reshape(tmp, (1,Nsamples+2))
    else:
        tmp = np.genfromtxt("results/mos1/real_10eVresult_T{0:1d}.dat".format(i))
    # print(tmp.shape)
    # if (i == 35):
    #     # continue
    #     tmp = np.reshape(tmp, (1,Nsamples+2))
    print(i,tmp.shape)
    real_data = np.concatenate((real_data,tmp), axis=0)

normE = real_data[:,1]
normDC = real_data[:,2::]
normN = np.linspace(-3e-5, 3e-5, Nsamples)

normsiggrid = np.zeros((len(normE),len(normN)))
for i in range(0, len(normE)):
    print(normE[i], end='\r')
    for j in range(0, len(normN)):
        sig = 0.0
        if (normDC[i][j] < 0.0):
            # Cycle through all fake results ***LONG RUNTIME***
            num = 0
            for k in range(0, nf):
                if (fake_deltaC[k] <= normDC[i][j]):
                    num += 1
            sig = 1.0 - num/nf

            # # Shortcut for quick plots
            # if (normDC[i][j] >= -1):
            #     sig = 0.1
            # elif (-1 > normDC[i][j] >= -4):
            #     sig = 0.6828
            # elif (-4 > normDC[i][j] >= -9):
            #     sig = 0.9546
            # elif (-9 > normDC[i][j]):
            #     sig = 0.9974

        normsiggrid[i][j] = sig
######### REAL ##################################


######### RESIDUALS #############################
residuals = np.genfromtxt('../fakes/mos1/residuals_1-10_mos1.qdp', skip_header=3)
energy = residuals[:,0]
energy_err = residuals[:,1]
ratio = residuals[:,2]
ratio_err = residuals[:,3]
######### RESIDUALS #############################


######### PLOT ##################################
### Three subplots sharing both x/y axes
f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)

### Plot the residuals first
ax1.axhline(1, color='k', dashes=[5,3], lw=1, zorder=0)
ax1.axvline(6.4/(1+srcz), color='k', dashes=[5,3], lw=1, zorder=0)
ax1.errorbar(energy, ratio, xerr=energy_err, yerr=ratio_err, fmt='+', ms=1, color='royalblue', ecolor='royalblue', lw=1, alpha=1.0, zorder=1)
ax1.set_ylabel('Data / Model')
# ax1.set_ylim(min(ratio)-0.125*min(ratio), max(ratio)+0.125*max(ratio))
ax1.set_ylim(min(ratio)-max(ratio_err), max(ratio)+max(ratio_err))
ax1.set_xlim(min(normE),max(normE))
ax1.set_title(srcname)
f.subplots_adjust(hspace=0)

### Plot up the contours
X,Y = np.meshgrid(normE,normN*10**5)
# for i in range(0,len(normE)):
#     for j in range(0, len(normN)):
#             if (normsiggrid[i][j] <= 0.6827):
#                 pc = 'grey'
#                 ps = 10
#             elif (0.6827 < normsiggrid[i][j] <= 0.9545):
#                 pc = 'b'
#                 ps = 20
#             elif (0.9545 < normsiggrid[i][j] <= 0.9973):
#                 pc = 'g'
#                 ps = 30
#             elif (0.9973 < normsiggrid[i][j]):
#                 pc = 'r'
#                 ps = 40
#             else:
#                 break
#             plt.scatter(normE[i]/(1+srcz), normN[j]*10**5, s=5, c=pc)
ax2.axhline(0, color='k', dashes=[5,3], lw=1)
ax2.axvline(6.4/(1+srcz), color='k', dashes=[5,3], lw=1, zorder=0)
ax2.contourf(X,Y,normsiggrid.T, colors=['#DCEDC8','#42B3D5','#1A237E'], levels=[0.6827,0.9545,0.9973,1.0], alpha=0.9)
# ax2.contourf(X,Y,normsiggrid.T, colors=['#FEEB65','#E4521B','#4D342F'], levels=[0.6827,0.9545,0.9973,1.0], alpha=0.9)
# ax2.contourf(X,Y,normsiggrid.T, colors=['#FFECB3','#E85285','#6A1B9A'], levels=[0.6827,0.9545,0.9973,1.0], alpha=0.9)
# ax2.scatter(X,Y, c=normsiggrid.T, s=10, marker='s', cmap='binary', norm=colors.LogNorm(vmin=0.1,vmax=1.0))
# ax2.contourf(X,Y,normsiggrid.T, cmap='viridis_r', norm=colors.LogNorm(vmin=0.6827, vmax=1.0), levels=[0.6827,0.9545,0.9973,1.0], alpha=0.9)


ax2.set_xlabel('Observed Energy [keV]')
ax2.set_ylabel('Line Normalisation\n[$\\times 10^{-5}$ photons cm$^{-2}$ s$^{-1}$]')
ax2.set_xlim(min(normE),max(normE))
ax2.set_ylim(bottom=-0.75, top=0.75)

plt.savefig("{0:s}_1-10_contours_10eV_mos1.png".format(srcname), bbox_inches='tight', dpi=600)
# plt.show()
######### PLOT ##################################


### Print out some useful / important information
print('\n****************************************')
print('Energy range: {0:.1f} - {1:.1f} keV'.format(min(normE),max(normE)))
# print('Number of iterations: {0:2d}'.format(trials))
print('****************************************\n')

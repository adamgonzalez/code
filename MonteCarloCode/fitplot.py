'''
Author: Adam Gonzalez
Descr.: This is a script to plot the output of the MC simulations testing UFO feature significance.
'''

import numpy as np
import matplotlib.pyplot as plt

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

### Redshift of the source
z = 0.0304051

### Read in the fake spectra results
fake_data = np.genfromtxt("fake_result100.dat")
fake = np.zeros((0,7))
for i in range(len(fake_data)):
    if ((fake_data[i,1] - 7.0) > -1e-3):
        fake = np.append(fake, [[fake_data[i,0], fake_data[i,1], fake_data[i,2], fake_data[i,3], fake_data[i,4], fake_data[i,5], fake_data[i,6]]], axis=0)

fake_spec = fake[:,0]
fake_E = fake[:,1] #; fake_E /= (1+z)
fake_N = fake[:,2]
fake_modfit = fake[:,3]
fake_basefit = fake[:,5]
fake_deltaC = fake_modfit-fake_basefit

### Determine mean of norms for each energy
ebins = np.arange(min(fake_E), max(fake_E)+0.1, 0.1)
# ebins = np.arange(7.0, max(fake_E)+0.1, 0.1)
trials = int(len(fake_spec)/len(ebins))
normgrid = np.zeros((trials,len(ebins)))
avggrid = np.zeros(len(ebins))
stdgrid = np.zeros(len(ebins))
for i in range(0,len(ebins)):
    c = 0
    for j in range(0,len(fake_N)):
        if (abs(fake_E[j] - ebins[i]) < 1e-3):
            normgrid[c][i] = fake_N[j]
            c += 1

for i in range(0,len(ebins)):
    avggrid[i] = np.average(normgrid[:,i])
    for j in range(0,trials):
        stdgrid[i] += (normgrid[j][i]-avggrid[i])**2.0
    stdgrid[i] /= (trials)
    stdgrid[i] = np.sqrt(stdgrid[i])


### Read in the actual data results
real = np.genfromtxt("real_result.dat")
real_E = real[:,1]
real_N = real[:,2]
real_modfit = real[:,3]
real_basefit = real[:,5]
real_deltaC = real_modfit-real_basefit

### Compute significance of detection(s) --- 1.00 -> 68.3 , 4.00 -> 95.4 , 9.00 -> 99.7
# num = 0
# sig = np.zeros((0,4))
# for i in range(0,len(fake_N)):
#     if (abs(fake_deltaC[i]) >= max(real_deltaC)):
#         sig = np.append(sig, [[real_deltaC, fake_deltaC[i], fake_E[i], fake_N[i]]], axis=0)
#         num += 1
# sig = 1.0 - num/trials
# print(sig)
# print(sig)

tot = len(fake_N)
siggrid = np.zeros(len(real_E))
for i in range(0, len(real_deltaC)):
    num = 0
    for j in range(0, len(fake_deltaC)):
        if (abs(fake_deltaC[j]) >= abs(real_deltaC[i])):
            num += 1
    sig = 1.0 - num/tot
    siggrid[i] = sig
# print(siggrid)

### Scale things up
avggrid *= 10**6
stdgrid *= 10**6
real_N *= 10**6


### Plot things up
plt.figure()

### Plot the MC results
plt.axhline(0, color='k', dashes=[5,3], lw=1)
# plt.fill_between(ebins, y1=avggrid-3*stdgrid, y2=avggrid+3*stdgrid, facecolor='k', alpha=0.2)
# plt.fill_between(ebins, y1=avggrid-2*stdgrid, y2=avggrid+2*stdgrid, facecolor='k', alpha=0.4)
# plt.fill_between(ebins, y1=avggrid-1*stdgrid, y2=avggrid+1*stdgrid, facecolor='k', alpha=0.5)
# # plt.fill_between(ebins, y1=avggrid-1*stdgrid, y2=avggrid-2*stdgrid, facecolor='k', alpha=0.5)
# # plt.fill_between(ebins, y1=avggrid+1*stdgrid, y2=avggrid+2*stdgrid, facecolor='k', alpha=0.5)
# # plt.fill_between(ebins, y1=avggrid-2*stdgrid, y2=avggrid-3*stdgrid, facecolor='k', alpha=0.3)
# # plt.fill_between(ebins, y1=avggrid+2*stdgrid, y2=avggrid+3*stdgrid, facecolor='k', alpha=0.3)
# # plt.scatter(ebins, avggrid, c='m', s=5, marker='s')
# # plt.errorbar(ebins, avggrid, yerr=3*stdgrid, fmt='s', ms=3, color='m', ecolor='r', lw=1)
# # plt.errorbar(ebins, avggrid, yerr=2*stdgrid, fmt='s', ms=3, color='m', ecolor='g', lw=1.5)
# # plt.errorbar(ebins, avggrid, yerr=1*stdgrid, fmt='s', ms=3, color='m', ecolor='b', lw=2)

### Plot the real data fit results
# # plt.plot(real_E, real_N, '-ok', ms=3, lw=1)
# plt.scatter(real_E, real_N, s=20, marker='o', c=real_deltaC, cmap='plasma', vmax=0.0)
# plt.colorbar(label=r'$\Delta C$')

plt.plot(real_E, real_N, color='k', lw=1, zorder=1)
for i in range(0, len(real_E)):
    if (siggrid[i] < 0.6827):
        pc = 'k'
        ps = 10
    elif (0.6827 <= siggrid[i] < 0.9545):
        pc = 'b'
        ps = 20
    elif (0.9545 <= siggrid[i] < 0.9973):
        pc = 'g'
        ps = 30
    elif (0.9973 <= siggrid[i]):
        pc = 'r'
        ps = 40
    else:
        break
    plt.scatter(real_E[i], real_N[i], s=20, marker='o', c=pc, zorder=2)
plt.xlabel(r'Energy [keV]')
plt.ylabel(r'Line Normalisation [$\times 10^{-6}$ photons cm$^{-2}$ s$^{-1}$]')

# plt.xlim(min(fake_E),max(fake_E))
plt.xlim(7.0,10.0)

# ### Plot significance?
# plt.axhline(y=0.6827, color='b', dashes=[5,3], lw=1)
# plt.axhline(y=0.9545, color='g', dashes=[5,3], lw=1)
# plt.axhline(y=0.9973, color='r', dashes=[5,3], lw=1)
# plt.plot(real_E, siggrid, '-ok', lw=1)
# plt.xlabel(r'Energy [keV]')
# plt.ylabel(r'Detection Significance')
# plt.ylim(bottom=0.5)
# plt.xlim(left=7.0, right=10.0)

# ### Histogram check
# # plt.axvline(0, color='k', dashes=[5,3], lw=1)
# # plt.hist(normgrid[:,0]*10**6)
# # plt.hist(fake_modfit)
# n, bins, patches = plt.hist(fake_deltaC, bins=int(abs(np.floor(min(fake_deltaC)))), range=[int(np.floor(min(fake_deltaC))),0])
# # print(bins)
# plt.yscale('log')
# plt.xlabel(r'$\Delta C$')
# plt.ylabel(r'Number of Instances')

plt.savefig("NormvsEnergy_Sig.png", bbox_inches='tight', dpi=300)
# plt.show()

### Print out some useful / important information
print('\n****************************************')
print('Energy range: {0:.1f} - {1:.1f} keV'.format(min(ebins),max(ebins)))
print('Number of iterations: {0:2d}'.format(trials))
print('****************************************\n')

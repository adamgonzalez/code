'''
Author: Adam Gonzalez
Descr.: This script is designed to provide an estimate of the expected C-statistic value and its variance given a model
        spectrum fit to a given observed spectrum. The approximations throughout can be found in Kaastra (2017), A&A, 605,
        A51. As described in the paper, you must use these functions to approximate the C-statistic in each energy bin i,
        then sum them up over the entire spectrum and compare to the C-statistic output of your model fit in order to
        determine the goodness of the model.
'''

import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
import sys
from scipy.special import erfinv
from astropy.io import fits

def exp_cstat(mu):
    # Function to compute the expected C-statistic value given mu the number of counts
    c_e = 0
    if (0 <= mu <= 0.5):
        c_e = -0.25*mu**3 + 1.38*mu**2 - 2*mu*math.log(mu)
    elif (0.5 < mu <= 2):
        c_e = -0.00335*mu**5 + 0.04259*mu**4 - 0.27331*mu**3 + 1.381*mu**2 - 2*mu*math.log(mu)
    elif (2 < mu <= 5):
        c_e = 1.019275 + 0.1345*mu**(0.461-0.9*math.log(mu))
    elif (5 < mu <= 10):
        c_e = 1.00624 + 0.604/(mu**1.68)
    elif(10 < mu):
        c_e = 1 + 0.1649/mu + 0.226/(mu**2)
    else:
        print('Error: no matching approximation')

    return c_e

def var_cstat(mu):
    # Function to compute the expected C-statistic variance given mu the number of counts
    c_v = 0
    if (0 <= mu <= 0.1):
        print('Error: 0 <= mu < 1 --- I have not included this...yet')
    elif (0.1 < mu <= 0.2):
        c_v = -262*mu**4 + 195*mu**3 - 51.24*mu**2 + 4.34*mu + 0.77005
    elif (0.2 < mu <= 0.3):
        c_v = 4.23*mu**2 - 2.8254*mu + 1.12522
    elif (0.3 < mu <= 0.5):
        c_v = -3.7*mu**3 + 7.328*mu**2 - 3.6926*mu + 1.20641
    elif (0.5 < mu <= 1.0):
        c_v = 1.28*mu**4 - 5.191*mu**3 + 7.666*mu**2 - 3.5446*mu + 1.15431
    elif (1 < mu <= 2):
        c_v = 0.1125*mu**4 - 0.641*mu**3 + 0.859*mu**2 + 1.0914*mu - 0.05748
    elif (2 < mu <= 3):
        c_v = 0.089*mu**3 - 0.872*mu**2 + 2.8422*mu - 0.67539
    elif (3 < mu <= 5):
        c_v = 2.12336 + 0.012202*mu**(5.717-2.6*math.log(mu))
    elif (5 < mu <= 10):
        c_v = 2.05159 + 0.331*mu**(1.343-math.log(mu))
    elif (10 < mu):
        c_v = 12/(mu**3) + 0.79/(mu**2) + 0.6747/mu + 2
    else:
        print('Error: no matching approximation')

    return c_v

# Open up the spectrum file to retrieve exact exposure time
parser = argparse.ArgumentParser(description='This script will read in the desired spectrum file and report the expected C-statistic and its variance.',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-s", "--specfilename", help="Name of spectrum", type=str, default='none')
parser.add_argument("-r", "--rmffilename", help="Name of rmf", type=str, default='none')
# parser.add_argument("-a", "--arffilename", help="Name of arf", type=str, default='none')
parser.add_argument("-low", "--lowelim", help="Low energy limit (keV)", type=float, default=0)
parser.add_argument("-hi", "--highelim", help="High energy limit (keV)", type=float, default=0)
if len(sys.argv) == 1:
    parser.print_help()
    print(help)
    sys.exit(1)
args = parser.parse_args()
specname = args.specfilename
rmfname = args.rmffilename
# arfname = args.arffilename
low_e = args.lowelim
high_e = args.highelim

### DO THINGS WITH RAW FILES ###############
# Read in spectrum file
specfitsfile = fits.open(specname)
# specfitsfile.info()
exp = specfitsfile['SPECTRUM'].header['EXPOSURE']
spec = specfitsfile['SPECTRUM'].data
specfitsfile.close()
spec_chan = spec.field('CHANNEL')
counts = spec.field('COUNTS')
flag = spec.field('GROUPING')

# Read in rmf
rmffitsfile = fits.open(rmfname)
# rmffitsfile.info()
rmf = rmffitsfile['EBOUNDS'].data
rmffitsfile.close()
rmf_chan = rmf.field('CHANNEL')
emin = rmf.field('E_MIN')
emax = rmf.field('E_MAX')

# # Read in arf
# arffitsfile = fits.open(arfname)
# # arffitsfile.info()
# arf = arffitsfile[1].data
# arffitsfile.close()
# elo = arf.field('ENERG_LO')
# ehi = arf.field('ENERG_HI')
# eff = arf.field('SPECRESP')

# Use the group flags to bin up the counts
bin_counts = np.array([])
bin_chan = np.array([])
for i in range(0,len(counts)):
    s = 0
    t = 0
    c = 1
    if (flag[i] == 1):
        s += counts[i]
        t += spec_chan[i]
        # print(spec_chan[i], counts[i])
        while (flag[i+c] == -1):
            s += counts[i+c]
            t += spec_chan[i]
            c += 1
            if ((i+c) >= len(counts)):
                break
        bin_counts = np.append(bin_counts,[s])
        bin_chan = np.append(bin_chan, [t/(c-1)])

# Compute the rmf energy vs chan relationship
eavg = np.zeros(len(rmf_chan))
for i in range(0,len(rmf_chan)):
    eavg[i] = (emax[i]+emin[i])/2
# coeffs = np.polyfit(rmf_chan,eavg,1)
# chan2energy = np.poly1d(coeffs)
coeffs = np.polyfit(eavg,rmf_chan,1)
energy2chan = np.poly1d(coeffs)

# Compute the corresponding bin limits for desired range
# bin_energy = chan2energy(bin_chan)
low_chan = math.floor(energy2chan(low_e))
high_chan = math.ceil(energy2chan(high_e))

# Call the C-stat functions to compute relevant quantities
cstat_exp = 0
cstat_var = 0
for i in range(0,len(bin_chan)):
    if (bin_chan[i] >= low_chan) and (bin_chan[i] <= high_chan):
        # print(bin_chan[i],bin_counts[i])
        cstat_exp += exp_cstat(bin_counts[i])
        cstat_var += var_cstat(bin_counts[i])
### DO THINGS WITH RAW FILES ###############


### DO THINGS WITH XSPEC ##############
# # Open up the temp.qdp file to obtain count information
# data = np.genfromtxt('temp.qdp', skip_header=3, comments='NO')
# counts = np.zeros(len(data))
#
# # Call the functions above to compute relevant quantities
# cstat_exp = 0
# cstat_var = 0
# for i in range(0,len(data)):
#     # counts[i] = data[i,-1]*exp*data[i,0]      # use model spectrum
#     counts[i] = data[i,2]*exp*data[i,0]       # use observed spectrum
#     cstat_exp += exp_cstat(counts[i])
#     cstat_var += var_cstat(counts[i])
### DO THINGS WITH XSPEC ##############


### FINAL CALCULATIONS AND OUTPUT ###############
# Compute the desired confidence interval
conf = np.array([0.683,0.90,0.954,0.997])
#lowconf = cstat_exp-np.sqrt(2)*erfinv(conf)*np.sqrt(cstat_var)
#uppconf = cstat_exp+np.sqrt(2)*erfinv(conf)*np.sqrt(cstat_var)

# Report results in terminal
print('\n****************************************')
print('Spectrum:  {}'.format(specname))
print('Exposure:  {0:.2f} sec'.format(exp))
print('Total counts:  {0:.1f}'.format(np.sum(bin_counts)))
print('Energy range:  {0:.1f} - {1:.1f} keV'.format(low_e,high_e))
print('')
print('C-stat expected:  {0:.2f}'.format(cstat_exp))
print('C-stat variance:  {0:.2f}'.format(np.sqrt(cstat_var)))
for i in range(0,len(conf)):
    lowconf = cstat_exp-np.sqrt(2)*erfinv(conf[i])*np.sqrt(cstat_var)
    uppconf = cstat_exp+np.sqrt(2)*erfinv(conf[i])*np.sqrt(cstat_var)
    print('{0:2d}% confidence:  [{1:.2f} , {2:.2f}]'.format(int(conf[i]*100),lowconf,uppconf))
print('****************************************\n')
### FINAL CALCULATIONS AND OUTPUT ###############

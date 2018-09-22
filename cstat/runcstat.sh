#!/bin/bash

# Author: Adam Gonzalez
# Descr.: This script calls the cstatcalc.py Python 3 script in order to estimate the expected value and variance of the
#         C-statistic of the loaded model spectrum. It should be sufficient to alter this file only.

# input spectrum information
datafile=PN_IZW1_opt.pha	# individual source spectrum
rmffile=PN_IZW1.rmf        	# rmf
low_e=0.4       # low energy limit (format 0.0 keV)
high_e=10.0     # high energy limit (format 0.0 keV)

# Run the Python 3 script
python ~/Dropbox/Graduate/PhD/cstat/cstatcalc.py -s ${datafile} -r ${rmffile} -low ${low_e} -hi ${high_e}

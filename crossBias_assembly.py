#########################################################
# This code will go into a directory, find the reuslts  #
# .txt files from meanstd.py and then plot the P(k) as  #
# well as the bias from the resolution.                 #
#
#########################################################

import numpy as np
import matplotlib.pyplot as plt
import sys, os
from matplotlib import rc
rc('text',usetex=True)
############ INPUT ######################################

root = '/mnt/xfs1/home/tcourt/test/'
folder1 = 'CDM_bias_0.6eV_64/'
folder2 = 'CDM_bias_0.6eV_128/'
folder3 = 'CDM_bias_0.6eV_256/'
folder4 = 'CDM_bias_0.6eV_32/'


####### Locate files #################

plt.title('Neutrino Mass 0.6eV CDM Cross Bias Split',fontsize=20)
plt.ylabel(r'$b_h(k)$',fontsize=16)
plt.xlabel(r'$k [h^{-1}Mpc]$',fontsize=16)
plt.xlim([6e-3,0.5])#,0.5,3.5])
plt.ylim([-3,6])

folders = [folder1, folder2, folder3, folder4]
colors  = ['r','g','b','magenta']
grids   = [64, 128, 256, 32]

b_h_total = root+folder1+"results_bhtt_c.txt"
k, b_ht, b_hterr = np.loadtxt(b_h_total,unpack=True)
plt.errorbar(k, b_ht, yerr=b_hterr, fmt='.-',c='black', label = 'Total Halos')

for folder,color,grid in zip(folders,colors,grids):

    b_h_above = root+folder+"results_bhaa_c.txt"
    b_h_below = root+folder+"results_bhbb_c.txt"

    k, b_ha, b_haerr = np.loadtxt(b_h_above,unpack=True)
    k, b_hb, b_hberr = np.loadtxt(b_h_below,unpack=True)

    plt.errorbar(k, b_ha, yerr=b_haerr, fmt='^-',c=color, label = 'Above %d'%grid)
    plt.errorbar(k, b_hb, yerr=b_hberr, fmt='v--',c=color, label = 'Below %d'%grid)
    

plt.xscale('log')
plt.legend(loc=1)

#plt.show()
plt.savefig('CDM_bias_0.6eV_assemblybiastotal.png')

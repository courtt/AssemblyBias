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
folder1 = 'random_0.0eV_64/'
folder2 = 'random_0.0eV_128/'
folder3 = 'random_0.0eV_256/'
folder4 = 'random_0.0eV_512/'


####### Locate files #################

plt.title('Neutrino Mass 0.0eV Total Random Split',fontsize=20)
plt.ylabel(r'$b_h(k)$',fontsize=16)
plt.xlabel(r'$k [h^{-1}Mpc]$',fontsize=16)
plt.xlim([6e-3,0.5])#,0.5,3.5])
plt.ylim([0,7])

folders = [folder1, folder2, folder3, folder4]
colors  = ['r','g','b','magenta']
grids   = [64, 128, 256, 512]

b_h_total = root+folder1+"results_bhtt.txt"
k, b_ht, b_hterr = np.loadtxt(b_h_total,unpack=True)
plt.errorbar(k, b_ht, yerr=b_hterr, fmt='.-',c='black', label = 'Total Halos')

for folder,color,grid in zip(folders,colors,grids):

    b_h_above = root+folder+"results_bhaa.txt"
    b_h_below = root+folder+"results_bhbb.txt"
    

    k, b_ha, b_haerr = np.loadtxt(b_h_above,unpack=True)
    k, b_hb, b_hberr = np.loadtxt(b_h_below,unpack=True)
    

    plt.errorbar(k, b_ha, yerr=b_haerr, fmt='^-',c=color, label = 'Above %d'%grid)
    plt.errorbar(k, b_hb, yerr=b_hberr, fmt='v--',c=color, label = 'Below %d'%grid)
    

plt.xscale('log')
plt.legend(loc=1)

#plt.show()
plt.savefig('r_0.0eV_assemblybiastotal.png')

"""
#64 
b_h_above1  = root+folder1+"results_bias_halos_above.txt"
b_h_below1  = root+folder1+"results_bias_halos_below.txt"
b_h_total1  = root+folder1+"results_bias_halos_total.txt"

#128
b_h_above2  = root+folder2+"results_bias_halos_above.txt"
b_h_below2  = root+folder2+"results_bias_halos_below.txt"
b_h_total2  = root+folder2+"results_bias_halos_total.txt"

#256
b_h_above3  = root+folder3+"results_bias_halos_above.txt"
b_h_below3  = root+folder3+"results_bias_halos_below.txt"
b_h_total3  = root+folder3+"results_bias_halos_total.txt"

#512
b_h_above4  = root+folder4+"results_bias_halos_above.txt"
b_h_below4  = root+folder4+"results_bias_halos_below.txt"
b_h_total4  = root+folder4+"results_bias_halos_total.txt"



###### Read files ####################
#64
k, b_h_a1, b_a_err1 = np.loadtxt(b_h_above1,unpack=True)
k, b_h_b1, b_b_err1 = np.loadtxt(b_h_below1,unpack=True)
k, b_h_t1, b_t_err1 = np.loadtxt(b_h_total1,unpack=True)

#128
k, b_h_a2, b_a_err2 = np.loadtxt(b_h_above2,unpack=True)
k, b_h_b2, b_b_err2 = np.loadtxt(b_h_below2,unpack=True)
k, b_h_t2, b_t_err2 = np.loadtxt(b_h_total2,unpack=True)

#256
k, b_h_a3, b_a_err3 = np.loadtxt(b_h_above3,unpack=True)
k, b_h_b3, b_b_err3 = np.loadtxt(b_h_below3,unpack=True)
k, b_h_t3, b_t_err3 = np.loadtxt(b_h_total3,unpack=True)

#512
k, b_h_a4, b_a_err4 = np.loadtxt(b_h_above4,unpack=True)
k, b_h_b4, b_b_err4 = np.loadtxt(b_h_below4,unpack=True)
k, b_h_t4, b_t_err4 = np.loadtxt(b_h_total4,unpack=True)


######### Plotting the data and error ###############

plt.title('Assembly Bias of Multiple Resolutions',fontsize=20)
plt.ylabel(r'$b_h(k)$',fontsize=16)
plt.xlabel(r'$k [h^{-1}Mpc]$',fontsize=16)
plt.xlim([6e-3,0.5])#,0.5,3.5])
plt.ylim([0,7])
#64
plt.errorbar(k, b_h_a1, yerr=b_a_err1, fmt='b^-', label = 'Above 64')
plt.errorbar(k, b_h_b1, yerr=b_b_err1, fmt='b^--', label = 'Below 64')
plt.errorbar(k, b_h_t1, yerr=b_t_err1, fmt='r.-', label = 'Total Halos')

#128
plt.errorbar(k, b_h_a2, yerr=b_a_err2, fmt='gv-', label = 'Above 128')
plt.errorbar(k, b_h_b2, yerr=b_b_err2, fmt='gv--', label = 'Below 128')
#plt.errorbar(k, b_h_t2, yerr=b_t_err2, fmt='r.-', label = 'Total Halos')

#256
plt.errorbar(k, b_h_a3, yerr=b_a_err3, fmt='cd-', label = 'Above 256')
plt.errorbar(k, b_h_b3, yerr=b_b_err3, fmt='cd--', label = 'Below 256')
#plt.errorbar(k, b_h_t3, yerr=b_t_err3, fmt='r.-', label = 'Total Halos')

#512
plt.errorbar(k, b_h_a4, yerr=b_a_err4, fmt='mD-', label = 'Above 512')
plt.errorbar(k, b_h_b4, yerr=b_b_err4, fmt='mD--', label = 'Below 512')
#plt.errorbar(k, b_h_t4, yerr=b_t_err4, fmt='r.-', label = 'Total Halos')



plt.xscale('log')
plt.legend(loc=1)

#plt.show()
plt.savefig('nassemblybiastotal.png')
"""





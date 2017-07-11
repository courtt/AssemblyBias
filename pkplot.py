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

for root_folder in ['0.15eV/','0.6eV/']:
    print root_folder

    for dims1 in [64,128,256,512,1024]:
        print dims1
        folder = 'results_cross_nu_%s_%i/'%(root_folder[:-1],dims1)

        ####### Locate files #################
        pkhaa   = root+folder+"results_haa.txt"
        pkhab   = root+folder+"results_hab.txt"
        pkhaabb = root+folder+"results_haabb.txt"
        pkhbb   = root+folder+"results_hbb.txt"
        pkhtt   = root+folder+"results_htt.txt"
        pkm     = root+folder+"results_m.txt"
        bhaa    = root+folder+"results_bhaa.txt"
        bhab    = root+folder+"results_bhab.txt"
        bhaabb  = root+folder+"results_bhaabb.txt"
        bhbb    = root+folder+"results_bhbb.txt"
        bhtt    = root+folder+"results_bhtt.txt"


        ###### Read files ####################
        k, pk_h_aa, aa_err     = np.loadtxt(pkhaa,unpack=True)
        k, pk_h_ab, ab_err     = np.loadtxt(pkhab,unpack=True)
        k, pk_h_aabb, aabb_err = np.loadtxt(pkhaabb,unpack=True)
        k, pk_h_bb, bb_err     = np.loadtxt(pkhbb,unpack=True)
        k, pk_h_tt, tt_err     = np.loadtxt(pkhtt,unpack=True)
        k, pk_m, m_err         = np.loadtxt(pkm,unpack=True)  
        k, b_h_aa, b_aa_err    = np.loadtxt(bhaa,unpack=True)
        k, b_h_ab, b_ab_err    = np.loadtxt(bhab,unpack=True)
        k, b_h_aabb, b_aabb_err= np.loadtxt(bhaabb,unpack=True)
        k, b_h_bb, b_bb_err    = np.loadtxt(bhbb,unpack=True)
        k, b_h_tt, b_tt_err    = np.loadtxt(bhtt,unpack=True)

        ######### Plotting the data and error ###############
        plt.title('Power Spectrum of Neutrino Mass %s at Resolution %i'%(root_folder[:-1],dims1),fontsize=20)
        plt.ylabel(r'$P(k) [(h^{-1}Mpc)^3]$',fontsize=16)
        plt.xlabel(r'$k [h^{-1}Mpc]$',fontsize=16)
        plt.xlim([6e-3,3.5])#,100,1e7])
        plt.errorbar(k, pk_h_aa, yerr=aa_err, fmt='.-', label = 'Halos Above Median')
        plt.errorbar(k, pk_h_bb, yerr=bb_err, fmt='.-', label = 'Halos Below Median')
        plt.errorbar(k, pk_h_ab, yerr=ab_err, fmt='.-', label = 'Halos Cross Median')
        plt.errorbar(k, pk_h_aabb, yerr=aabb_err, fmt='.-', label = 'Halos Above+Below+Cross') 
        plt.errorbar(k, pk_h_tt, yerr=tt_err, fmt='.-', label = 'Total Halos')
        plt.errorbar(k, pk_m, yerr=m_err, fmt='.-', label = 'Dark Matter')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc=0)
        
        #plt.show()
        plt.savefig('Nu_%s_%i.png'%(root_folder[:-1],dims1))
        plt.clf()
        
        plt.title('Assembly Bias of Neutrino Mass %s at Resolution %i'%(root_folder[:-1],dims1),fontsize=20)
        plt.ylabel(r'$b_h(k)$',fontsize=16)
        plt.xlabel(r'$k [h^{-1}Mpc]$',fontsize=16)
        plt.xlim([6e-3,0.5])#,0.5,3.5])
        plt.ylim([0.5,3.5])
        plt.errorbar(k, b_h_aa, yerr=b_aa_err, fmt='.-', label = 'Halos Above Median')
        plt.errorbar(k, b_h_bb, yerr=b_bb_err, fmt='.-', label = 'Halos Below Median')
        plt.errorbar(k, b_h_ab, yerr=b_ab_err, fmt='.-', label = 'Halos Cross Median')
        plt.errorbar(k, b_h_aabb, yerr=b_aabb_err, fmt='.-', label = 'Halos Above+Below+Cross')
        plt.errorbar(k, b_h_tt, yerr=b_tt_err, fmt='.-', label = 'Total Halos')
        plt.xscale('log')
        plt.legend(loc=0)

        #plt.show()
        plt.savefig('Nu_%s_%i.png'%(root_folder[:-1],dims1))






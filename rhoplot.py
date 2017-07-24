######################## 100split.py ##################################
# This code will loop through all of the simulations, calculate the   #
# density of each bin and then calculate the local median for a number#
# of mass bins. It then sorts the values to be above or below the med #
# and calculates P(k) above and below the med and puts that data in a #
# .txt file. It also computes the combined P(k), unsorted.            #
#######################################################################

import numpy as np
import readsnap, readfof
import MAS_library as MASL
import Pk_library as PKL
import sys,os
import matplotlib.pyplot as plt

#root_sims = '/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2/0.0eV/'
root_sims = '/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2/0.6eV/'
root_output = '/mnt/xfs1/home/tcourt/test/'
############################## INPUT ##################################
dims1 = 1024
dims = 512
bins = 40
realization=1
#######################################################################
# name of output folder
folder = root_output+'results_cross_nu_0.15eV__%i/'%dims1

# create output folder if it doesnt exist
if not(os.path.exists(folder)):
    os.system('mkdir '+folder)
    
    
#while count < 34:
print 'REALIZATION NUMBER', realization
snapshot_fname = root_sims+'%i/snapdir_004/snap_004'%realization
snapdir = root_sims+'%i/'%realization
snapnum = 4

# read snapshot head and obtain BoxSize, Omega_m and Omega_L    
print '\nREADING SNAPSHOTS PROPERTIES'
head     = readsnap.snapshot_header(snapshot_fname)
BoxSize  = head.boxsize/1e3  #Mpc/h                  
Nall     = head.nall
Masses   = head.massarr*1e10 #Msun/h                                     
Omega_m  = head.omega_m
Omega_l  = head.omega_l
redshift = head.redshift
Hubble   = 100.0*np.sqrt(Omega_m*(1.0+redshift)**3+Omega_l)#km/s/(Mpc/h)       
h        = head.hubble
z        = '%.3f'%redshift

# read positions and velocites
pos_nu = readsnap.read_block(snapshot_fname,"POS ",parttype=2)/1e3 #Mpc/h
pos_c = readsnap.read_block(snapshot_fname,"POS ",parttype=1)/1e3 #Mpc/h
#vel = readsnap.read_block(snapshot_fname,"VEL ",parttype=1)     #km/s         
print 'positions calculated' 
# compute density
rhoFunc_nu = np.zeros((dims1,dims1,dims1), dtype=np.float32)
MASL.CIC(pos_nu,rhoFunc_nu,BoxSize)

rhoFunc_c = np.zeros((dims1,dims1,dims1), dtype=np.float32)
MASL.CIC(pos_c,rhoFunc_c,BoxSize)

print ' densities found'

# read positions and velocities of halos
FoF   = readfof.FoF_catalog(snapdir,snapnum,long_ids=False,
                            swap=False,SFR=False)
pos_h = FoF.GroupPos/1e3            #Mpc/h
mass  = FoF.GroupMass*1e10          #Msun/h
vel_h = FoF.GroupVel*(1.0+redshift) #km/s

# compute the density in the positions of the halos
rho_nu = np.zeros(len(pos_h), dtype=np.float32)
MASL.CIC_interp(rhoFunc_nu,BoxSize,pos_h,rho_nu)

rho_c = np.zeros(len(pos_h), dtype=np.float32)
MASL.CIC_interp(rhoFunc_c,BoxSize,pos_h,rho_c)

plt.title('Density plot of Resolution %i'%dims1,fontsize=20)
plt.ylabel(r'$\rho_\nu$',fontsize=16)
plt.xlabel(r'$\rho_c$',fontsize=16)
plt.xscale('log')
plt.yscale('log')
plt.plot(rho_c,rho_nu,'bo')
#plt.show()
plt.savefig('rho_%i.png'%dims1)

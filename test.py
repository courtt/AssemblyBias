import numpy as np
import readsnap
import sys,os

snapshot_fname = '/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2/0.0eV/1/snapdir_004/snap_004'

print('SNAPSHOT #' + snapshot_fname[-1])

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

pos = readsnap.read_block(snapshot_fname,"POS ",parttype=1)/1e3 #Mpc/h
vel = readsnap.read_block(snapshot_fname,"VEL ",parttype=1)     #km/s


# to know the position of the first particle type
print pos[0]

#to know the position of the last particle type
print pos[-1]

# to know the velocity of the first particle
print vel[0]
#print Omega_m
#print Omega_l
#print BoxSize

maxv=np.amax(vel[:])
maxvindex=np.argmax(vel[:])
print "maxv", maxv
print "maxvindex", maxvindex
print "Hubble", Hubble
#print h
print 'the redshift is ', z

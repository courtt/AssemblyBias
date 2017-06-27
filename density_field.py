#################### density_field.py #################################
# This code calculates the density field for a grid cell of halos by  #
# determining the cell in which a halo is and calculating the desnsity#
# in that cell. It then plots a density vs. mass and shows the best   #
# fit line of the data to first and second order                      #
#                                                                     #
#######################################################################

import numpy as np
import readsnap, readfof
import MAS_library as MASL
import sys,os
import matplotlib.pyplot as plt

############################## INPUT ##################################
# for dark matter
snapshot_fname = '/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2/0.0eV/1/snapdir_004/snap_004'

# for halos
snapdir = '/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2/0.0eV/1/'
snapnum = 4

dims = 512
#######################################################################

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


rho = np.zeros((dims,dims,dims), dtype=np.float32)
MASL.CIC(pos,rho,BoxSize)

print len(pos)
print np.sum(rho, dtype=np.float64)

print rho[0,1,5]  #for any box, the density can be calculated
print BoxSize     # how large the box is as a float




# read positions and velocities of halos
FoF   = readfof.FoF_catalog(snapdir,snapnum,long_ids=False,
                            swap=False,SFR=False)
pos_h = FoF.GroupPos/1e3            #Mpc/h
mass  = FoF.GroupMass*1e10          #Msun/h
vel_h = FoF.GroupVel*(1.0+redshift) #km/s

print pos_h[0]
print len(pos_h)
print len(mass)
x=[]
for i in mass:
    x.append(i)

print len(x)
y=[]
boxdiv=BoxSize/dims     #width of each box 
counter=0
for i in pos_h:
    boxnum= np.trunc(pos_h[counter]/boxdiv)  #places each positional vector into respective box
    #print boxnum
    y1= boxnum.astype(int)
    #print x1
    #print boxnum.astype(int)
    #print rho[130, 172,40]
    y.append(rho[y1[0],y1[1],y1[2]])
    counter += 1

print len(y)

#################### Plot Input#################################
x_min=min(x)
x_max=max(x)
y_min=min(y)
y_max=max(y)




coefficients = np.polyfit(x, y, 1)
co_square = np.polyfit(x, y, 2)
polynomial = np.poly1d(coefficients)
poly_square= np.poly1d(co_square)
ys = polynomial(x)
ys2= poly_square(x)
print coefficients
print polynomial
print poly_square
plt.plot(x, y, 'o')
plt.plot(x, ys,linestyle='-',lw=3)
plt.plot(x, ys2,linestyle='--',lw=3,color='magenta')
ax = plt.gca()
#ax.scatter(x,y)
plt.axis([x_min,x_max,y_min,y_max])
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('log(Mass of Halo)')
plt.ylabel('log(Density Field)')
plt.show()



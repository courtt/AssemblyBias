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

root_sims = '/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2/0.0eV/'
root_output = '/mnt/xfs1/home/tcourt/test/'
############################## INPUT ##################################
dims1 = 1024
dims = 512
bins = 40
#######################################################################

# name of output folder
folder = root_output+'results_%i/'%dims1

# create output folder if it doesnt exist
if not(os.path.exists(folder)):
    os.system('mkdir '+folder)

count=49
while count < 50:
    print 'REALIZATION NUMBER', count
    snapshot_fname = root_sims+'%i/snapdir_004/snap_004'%count
    snapdir = root_sims+'%i/'%count
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
    pos = readsnap.read_block(snapshot_fname,"POS ",parttype=1)/1e3 #Mpc/h
    vel = readsnap.read_block(snapshot_fname,"VEL ",parttype=1)     #km/s         
    
    # compute density
    rhoFunc = np.zeros((dims1,dims1,dims1), dtype=np.float32)
    MASL.CIC(pos,rhoFunc,BoxSize)
    
    
    # read positions and velocities of halos
    FoF   = readfof.FoF_catalog(snapdir,snapnum,long_ids=False,
                                swap=False,SFR=False)
    pos_h = FoF.GroupPos/1e3            #Mpc/h
    mass  = FoF.GroupMass*1e10          #Msun/h
    vel_h = FoF.GroupVel*(1.0+redshift) #km/s
    
    # compute the density in the positions of the halos
    rho = np.zeros(len(pos_h), dtype=np.float32)
    MASL.CIC_interp(rhoFunc,BoxSize,pos_h,rho)



############### SORT HALOS BY DENSITY #######################
    # create bins for local median     
    intervals= np.logspace(np.log10(np.min(mass)),np.log10(np.max(mass)),bins) # creates 40 bins between Mmin and Mmax 
    delta_m = np.log10(intervals[1])-np.log10(intervals[0])          # Calculates the width of each bin using 1st and 0th spot
    index = np.int_((np.log10(mass)-np.log10(intervals[0]))/delta_m)  # Creates an index of bin numbers corresponding to masses
        
    rho_median= np.zeros(bins,dtype=np.float64)
        
    
        
    # create empty arrays since they vary in length
    massAbove = []
    rhoAbove = []
    indexAbove = [] 
    massBelow = []
    rhoBelow = []
    indexBelow = []
        

    for i in xrange(bins):
        index2 = np.where(index==i)[0]       
        rho_median[i] = np.median(rho[index2]) # index2 stores the indeces of what's in each bin
        rho_diff = rho[index2]- rho_median[i]
    
        for j in xrange(len(rho_diff)):
            if rho_diff[j] >= 0:
                massAbove.append(mass[index2[j]])
                rhoAbove.append(rho[index2[j]])
                indexAbove.append(index2[j])
            else:
                massBelow.append(mass[index2[j]])
                rhoBelow.append(rho[index2[j]])    
                indexBelow.append(index2[j])

    massAbove = np.array(massAbove)
    rhoAbove = np.array(rhoAbove)
    indexAbove = np.array(indexAbove)
    massBelow = np.array(massBelow)
    rhoBelow = np.array(rhoBelow)
    indexBelow = np.array(indexBelow)


    print 'Total above', len(indexAbove)
    print 'Total below', len(indexBelow)

################ COMPUTE DENSITY AND P(k) #############################

    # compute density field of halos above
    delta_h = np.zeros((dims,dims,dims), dtype=np.float32)
    MASL.CIC(pos_h[indexAbove],delta_h,BoxSize)
    delta_h /= np.mean(delta_h,dtype=np.float64);  delta_h -= 1.0
   
    print '%.3e < delta_h < %.3e'%(np.min(delta_h),np.max(delta_h))
    
    # compute power spectrum of halos above
    Pk_h = PKL.Pk(delta_h, BoxSize, axis=0, MAS='CIC', threads=1)
    Pk_h.Pk[:,0] = Pk_h.Pk[:,0] - BoxSize**3*1.0/(len(pos_h[indexAbove]))
    fout = folder+'Pk_h_above_local%i.txt'%count
    np.savetxt(fout, np.transpose([Pk_h.k3D,Pk_h.Pk[:,0]]))


    # compute density field of halos below
    delta_h = np.zeros((dims,dims,dims), dtype=np.float32)
    MASL.CIC(pos_h[indexBelow],delta_h,BoxSize)
    delta_h /= np.mean(delta_h,dtype=np.float64);  delta_h -= 1.0
    
    print '%.3e < delta_h < %.3e'%(np.min(delta_h),np.max(delta_h))
    
    # compute power spectrum of halos below
    Pk_h = PKL.Pk(delta_h, BoxSize, axis=0, MAS='CIC', threads=1)
    Pk_h.Pk[:,0] = Pk_h.Pk[:,0] - BoxSize**3*1.0/(len(pos_h[indexBelow]))
    fout = folder+'Pk_h_below_local%i.txt'%count
    np.savetxt(fout, np.transpose([Pk_h.k3D,Pk_h.Pk[:,0]]))
    

    # compute density field of all halos 
    index_total = np.hstack([indexAbove,indexBelow])
    delta_h = np.zeros((dims,dims,dims), dtype=np.float32)
    MASL.CIC(pos_h[index_total],delta_h,BoxSize)
    delta_h /= np.mean(delta_h,dtype=np.float64);  delta_h -= 1.0
    
    print '%.3e < delta_h < %.3e'%(np.min(delta_h),np.max(delta_h))
    
    # compute power spectrum of all halos 
    Pk_h = PKL.Pk(delta_h, BoxSize, axis=0, MAS='CIC', threads=1)
    Pk_h.Pk[:,0] = Pk_h.Pk[:,0] - BoxSize**3*1.0/(len(pos_h[index_total]))
    fout = folder+'Pk_h_total_local%i.txt'%count
    np.savetxt(fout, np.transpose([Pk_h.k3D,Pk_h.Pk[:,0]]))
    
    count+=1

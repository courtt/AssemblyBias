######################## sliceplot.py #################################
# This code will loop through all of the simulations, calculate the   #
# density of each bin and then calculate the local median for a number#
# of mass bins. It then sorts the values to be above or below the med #
# and calculates P(k) above and below the med and puts that data in a #
# .txt file. It also computes the combined P(k), unsorted.            #
#######################################################################

#from mpi4py import MPI
import numpy as np
import readsnap, readfof
import MAS_library as MASL
import Pk_library as PKL
import sys,os
import matplotlib.pyplot as plt

def clustering(root_sims,folder,count,dims1,dims,bins,ptype):#ptype

    print 'REALIZATION NUMBER', count
    snapshot_fname = root_sims+'%i/snapdir_004/snap_004'%count
    snapdir = root_sims+'%i/'%count
    snapnum = 4

    # read snapshot head and obtain BoxSize, Omega_m and Omega_L    
    print '\nREADING SNAPSHOTS PROPERTIES'
    head     = readsnap.snapshot_header(snapshot_fname)
    Masses   = head.massarr*1e10 #Msun/h 
    BoxSize  = head.boxsize/1e3  #Mpc/h                  
    redshift = head.redshift
    z        = '%.3f'%redshift
    
    # read positions of cdm
    pos = readsnap.read_block(snapshot_fname,"POS ",parttype=ptype)/1e3 #Mpc/h
    rhoFunc = np.zeros((dims1,dims1,dims1), dtype=np.float32)
    MASL.CIC(pos, rhoFunc, BoxSize)
    
    # read positions and velocities of halos
    FoF   = readfof.FoF_catalog(snapdir,snapnum,long_ids=False,
                                swap=False,SFR=False)
    pos_h = FoF.GroupPos/1e3            #Mpc/h
    mass  = FoF.GroupMass*1e10          #Msun/h
    
    # compute the density in the positions of the halos
    rho_h = np.zeros(len(pos_h), dtype=np.float32)
    MASL.CIC_interp(rhoFunc, BoxSize, pos_h, rho_h)

    ############### SORT HALOS BY DENSITY #######################
    # create bins for local median     
    Mmin, Mmax = np.min(mass), np.max(mass)
    intervals  = np.logspace(np.log10(Mmin), np.log10(Mmax), bins) 
    delta_mass = np.log10(intervals[1])-np.log10(intervals[0]) 
    index      = np.int_((np.log10(mass)-np.log10(intervals[0]))/delta_mass) 
    
    
    # create empty arrays since they vary in length
    indexAbove, indexBelow = [],[]
    massAbove, massBelow, rhoAbove, rhoBelow = [],[],[],[]
    
    # for each bin, find halos above and below the median
    for i in xrange(bins):
        index2 = np.where(index==i)[0]
        rho_median = np.median(rho_h[index2]) # index2 stores the indeces of what's in each bin
        rho_diff = rho_h[index2]- rho_median
        
        for j in xrange(len(rho_diff)):
            if rho_diff[j] >= 0:
                #massAbove.append(mass[index2[j]])
                #rhoAbove.append(rho[index2[j]])
                indexAbove.append(index2[j])
            else:
                #massBelow.append(mass[index2[j]])
                #rhoBelow.append(rho[index2[j]])    
                indexBelow.append(index2[j])
                
    #massAbove  = np.array(massAbove);  massBelow  = np.array(massBelow)
    #rhoAbove   = np.array(rhoAbove);   rhoBelow   = np.array(rhoBelow)
    indexAbove = np.array(indexAbove);  indexBelow = np.array(indexBelow)
                    
    print 'Halos above = %d\nHalos below = %d'%(len(indexAbove),len(indexBelow))

    sliceIndex = np.where(pos_h[:,2]<5.0)[0]
    print len(sliceIndex)

    sliceAbove = np.where(pos_h[indexAbove,2]<5.0)[0]
    #sliceAbove = np.where(pos_h[indexAbove] == pos_h[sliceIndex])[0]
    print len(sliceAbove)

    sliceBelow = np.where(pos_h[indexBelow,2]<5.0)[0]
    #sliceBelow = np.where(pos_h[indexBelow] == pos_h[sliceIndex])[0]
    print len(sliceBelow)

    print len(sliceBelow)+len(sliceAbove)
    
    nmass = root_sims[-6:-1]
    
    if ptype == 1: splittype = 'CDM'
    if ptype == 2: splittype = 'Neutrino'


    plt.clf()
    plt.scatter(pos_h[sliceIndex,0],pos_h[sliceIndex,1],s=5*(mass[sliceIndex]/np.min(mass[sliceIndex])))
    plt.axis([0,1000,0,1000])
    plt.xlabel('x [Mpc/h]')
    plt.ylabel('y [Mpc/h]')
    plt.title('Slice of %s Resolution, No Split'%nmass)
    #plt.show()
    plt.savefig('halos_total_%s_%s.png'%(nmass,splittype))
    
    plt.clf()
    plt.axis([0,1000,0,1000])
    plt.xlabel('x [Mpc/h]')
    plt.ylabel('y [Mpc/h]')
    plt.scatter(pos_h[indexAbove[sliceAbove],0],pos_h[indexAbove[sliceAbove],1],c='r',s=5*(mass[indexAbove[sliceAbove]])/np.min(mass[sliceIndex]),alpha=0.5)
    plt.scatter(pos_h[indexBelow[sliceBelow],0],pos_h[indexBelow[sliceBelow],1],c='b',s=5*(mass[indexBelow[sliceBelow]])/np.min(mass[sliceIndex]),alpha = 0.5)
    plt.title('Slice of %s Resolution, %s Split'%(nmass,splittype))
    plt.savefig('halos_split_%s_%s.png'%(nmass,splittype))
    #plt.show()



'''
################ COMPUTE DENSITY AND P(k) #############################

    # compute density field of halos above
    delta_haa = np.zeros((dims,dims,dims), dtype=np.float32)
    MASL.CIC(pos_h[indexAbove],delta_haa,BoxSize)
    delta_haa /= np.mean(delta_haa,dtype=np.float64);  delta_haa -= 1.0
    print '%.3e < delta_haa < %.3e'%(np.min(delta_haa),np.max(delta_haa))

    # compute density field of halos below
    delta_hbb = np.zeros((dims,dims,dims), dtype=np.float32)
    MASL.CIC(pos_h[indexBelow],delta_hbb,BoxSize)
    delta_hbb /= np.mean(delta_hbb,dtype=np.float64);  delta_hbb -= 1.0
    print '%.3e < delta_hbb < %.3e'%(np.min(delta_hbb),np.max(delta_hbb))
    
    # compute density field of all halos 
    index_total = np.hstack([indexAbove,indexBelow])
    delta_h = np.zeros((dims,dims,dims), dtype=np.float32)
    MASL.CIC(pos_h, delta_h, BoxSize)
    delta_h /= np.mean(delta_h,dtype=np.float64);  delta_h -= 1.0
    print '%.3e < delta_h < %.3e'%(np.min(delta_h),np.max(delta_h))
'''


root_total = '/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2/'
root_output = '/mnt/xfs1/home/tcourt/test/'
############################## INPUT ##################################
dims1 = 512 #grid size to estimate overdensities at halos locations
dims  = 512  #grid size to compute halo Pk
bins  = 40
i = 1
folder = 'sliceplot'
files = 100
#ptype = 2  #1-CDM, 2-NU
#######################################################################


for root_folder in ['0.0eV/','0.15eV/','0.6eV/']:
    root_sims = root_total + root_folder
    print root_sims[-6:-1]
    print root_folder

    for ptype in [1,2]:
        print ptype
        if root_folder=='0.0eV/' and ptype==2:
            continue

        clustering(root_sims,folder,i,dims1,dims,bins,ptype)

        # create folder if it does not exist
    #if myrank==0:
    #    if not(os.path.exists(folder+'Pk_haa_100.txt')):  
    #        os.system('mkdir '+folder)
    #    else:
    #        print 'FOLDER ALREADY COMPLETED'
    #        continue

        # find the numbers each cpu will work
    #numbers = np.where(np.arange(files)%nprocs==myrank)[0]+1

    # do a loop over all realizations of each cpu
    #for i in numbers:

    #    print 'cpu %d working with realization %2d'%(myrank,i)



######################## 100split.py ##################################
# This code will loop through all of the simulations, calculate the   #
# density of each bin and then calculate the local median for a number#
# of mass bins. It then sorts the values to be above or below the med #
# and calculates P(k) above and below the med and puts that data in a #
# .txt file. It also computes the combined P(k), unsorted.            #
#######################################################################

from mpi4py import MPI
import numpy as np
import readsnap, readfof
import MAS_library as MASL
import Pk_library as PKL
import sys,os

###### MPI DEFINITIONS ###### 
comm=MPI.COMM_WORLD
nprocs=comm.Get_size()
myrank=comm.Get_rank()


def clustering(root_sims,folder,count,ptype,dims1,dims,bins):

    print 'REALIZATION NUMBER', count
    snapshot_fname = root_sims+'%i/snapdir_004/snap_004'%count
    snapdir = root_sims+'%i/'%count
    snapnum = 4

    # read snapshot head and obtain BoxSize, Omega_m and Omega_L    
    print '\nREADING SNAPSHOTS PROPERTIES'
    head     = readsnap.snapshot_header(snapshot_fname)
    BoxSize  = head.boxsize/1e3  #Mpc/h                  
    redshift = head.redshift
    z        = '%.3f'%redshift
    
    # read positions 
    pos = readsnap.read_block(snapshot_fname,"POS ",parttype=ptype)/1e3 #Mpc/h

    # compute density field
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
    
    # compute auto- and cross-power spectra
    Pk = PKL.XPk([delta_haa,delta_hbb], BoxSize, axis=0, 
                 MAS=['CIC','CIC'], threads=1)

    # subtract shot-noise from halos
    Pk.Pk[:,0,0] = Pk.Pk[:,0,0] - BoxSize**3*1.0/len(pos_h[indexAbove])
    Pk.Pk[:,0,1] = Pk.Pk[:,0,1] - BoxSize**3*1.0/len(pos_h[indexBelow])

    # find coefficients for total Pk
    c1 = len(indexAbove)*1.0/len(pos_h)
    c2 = len(indexBelow)*1.0/len(pos_h)

    # save results to file
    fout1 = folder+'Pk_haa_%i.txt'%count
    fout2 = folder+'Pk_hbb_%i.txt'%count
    fout3 = folder+'Pk_hab_%i.txt'%count
    fout4 = folder+'Pk_haa+hbb+2hab_%i.txt'%count
    np.savetxt(fout1, np.transpose([Pk.k3D, Pk.Pk[:,0,0]]))
    np.savetxt(fout2, np.transpose([Pk.k3D, Pk.Pk[:,0,1]]))
    np.savetxt(fout3, np.transpose([Pk.k3D, Pk.XPk[:,0,0]]))
    np.savetxt(fout4, 
               np.transpose([Pk.k3D, c1**2*Pk.Pk[:,0,0]+c2**2*Pk.Pk[:,0,1]+2*c1*c2*Pk.XPk[:,0,0]]))

    # compute density field of all halos 
    index_total = np.hstack([indexAbove,indexBelow])
    delta_h = np.zeros((dims,dims,dims), dtype=np.float32)
    MASL.CIC(pos_h, delta_h, BoxSize)
    delta_h /= np.mean(delta_h,dtype=np.float64);  delta_h -= 1.0
    print '%.3e < delta_h < %.3e'%(np.min(delta_h),np.max(delta_h))
    
    # compute power spectrum of all halos 
    Pk_h = PKL.Pk(delta_h, BoxSize, axis=0, MAS='CIC', threads=1)
    Pk_h.Pk[:,0] = Pk_h.Pk[:,0] - BoxSize**3*1.0/len(pos_h)
    fout = folder+'Pk_htt_%i.txt'%count
    np.savetxt(fout, np.transpose([Pk_h.k3D,Pk_h.Pk[:,0]]))


root_total = '/mnt/ceph/users/fvillaescusa/Neutrino_simulations/Sims_Dec16_2/'
root_output = '/mnt/xfs1/home/tcourt/test/'
############################## INPUT ##################################
#dims1 = 1024 #grid size to estimate overdensities at halos locations
dims  = 512  #grid size to compute halo Pk
bins  = 40

files = 100
#ptype = 2  #1-CDM, 2-NU
#######################################################################


for root_folder in ['0.0eV/','0.15eV/','0.6eV/']:
    root_sims = root_total + root_folder
    print root_folder
    for dims1 in [64,128,256,512,1024]:
        print dims1
        for ptype in [1,2]:
            print ptype
            if root_folder=='0.0eV/' and ptype==2:
                continue


            # name of output folder
            if ptype==1: folder = root_output+'results_cross_c_%s_%i/'%(root_folder[:-1],dims1)
            if ptype==2: folder = root_output+'results_cross_nu_%s_%i/'%(root_folder[:-1],dims1)
            

            # create folder if it does not exists
            if myrank==0:
                if not(os.path.exists(folder+'Pk_haa_100.txt')):  
                    os.system('mkdir '+folder)
                else:
                    print 'FOLDER ALREADY COMPLETED'
                    continue

            # find the numbers each cpu will work
            numbers = np.where(np.arange(files)%nprocs==myrank)[0]+1

            # do a loop over all realizations of each cpu
            for i in numbers:

                print 'cpu %d working with realization %2d'%(myrank,i)
                clustering(root_sims,folder,i,ptype,dims1,dims,bins)


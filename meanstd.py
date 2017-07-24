##################################################################
# This program takes the 100 files made by powerspectrum.py and  #
# it finds the mean and standard deviation of each of them then  #
# puts the results into a .txt file for plotting in GNUplot      #
#                                                                #
##################################################################

import numpy as np
#import matplotlib.pyplot as plt
import sys, os
############### setting up the array ##############################

root = '/mnt/xfs1/home/tcourt/test/'

for root_folder in ['0.0eV/','0.6eV/']:
    print root_folder

    for dims1 in [64,128,256,512,1024]:
        print dims1
        folder = 'random_%s_%i/'%(root_folder[:-1],dims1)
    
        k, Pk = np.loadtxt(root+folder+'Pk_htt_1.txt',unpack=True)

        bins = len(Pk)
        Pk_haa_matrix   = np.zeros((100,bins),dtype=np.float32)
        Pk_hab_matrix   = np.zeros((100,bins),dtype=np.float32)
        Pk_haabb_matrix = np.zeros((100,bins),dtype=np.float32)
        Pk_hbb_matrix   = np.zeros((100,bins),dtype=np.float32)
        Pk_htt_matrix   = np.zeros((100,bins),dtype=np.float32)
        Pk_m_matrix     = np.zeros((100,bins),dtype=np.float32)
        b_haa_matrix    = np.zeros((100,bins),dtype=np.float32)
        b_hab_matrix    = np.zeros((100,bins),dtype=np.float32)
        b_haabb_matrix  = np.zeros((100,bins),dtype=np.float32)
        b_hbb_matrix    = np.zeros((100,bins),dtype=np.float32)
        b_htt_matrix    = np.zeros((100,bins),dtype=np.float32)
   
    ######################### Halos Above #################################

        for i in xrange(100):

            ##### find file names #####
            file_aa   = root+folder+'Pk_haa_%i.txt'%(i+1)
            file_aabb = root+folder+'Pk_haa+hbb+2hab_%i.txt'%(i+1)
            file_ab   = root+folder+'Pk_hab_%i.txt'%(i+1)
            file_bb   = root+folder+'Pk_hbb_%i.txt'%(i+1)
            file_tt   = root+folder+'Pk_htt_%i.txt'%(i+1)
            file_m    = root+'meanstd/Pk%i.txt'%(i+1)

            #### read files ####
            k, Pk_h_aa   = np.loadtxt(file_aa,unpack=True)
            k, Pk_h_bb   = np.loadtxt(file_bb,unpack=True)
            k, Pk_h_tt   = np.loadtxt(file_tt,unpack=True)
            k, Pk_h_aabb = np.loadtxt(file_aabb,unpack=True)
            k, Pk_h_ab   = np.loadtxt(file_ab,unpack=True)
            k, Pk_m      = np.loadtxt(file_m,unpack=True)
            
            #### fill big matrix ####
            Pk_haa_matrix[i]   = Pk_h_aa
            Pk_hab_matrix[i]   = Pk_h_ab
            Pk_haabb_matrix[i] = Pk_h_aabb
            Pk_hbb_matrix[i]   = Pk_h_bb
            Pk_htt_matrix[i]   = Pk_h_tt
            Pk_m_matrix[i]     = Pk_m
            b_haa_matrix[i]    = np.sqrt(Pk_h_aa/Pk_m)
            b_hab_matrix[i]    = np.sqrt(Pk_h_ab/Pk_m)
            b_haabb_matrix[i]  = np.sqrt(Pk_h_aabb/Pk_m)
            b_hbb_matrix[i]    = np.sqrt(Pk_h_bb/Pk_m)
            b_htt_matrix[i]    = np.sqrt(Pk_h_tt/Pk_m)
           
            # write results to file
            f1 = open(root+folder+"results_haa.txt","w")
            f2 = open(root+folder+"results_hab.txt","w")
            f3 = open(root+folder+"results_haabb.txt","w")
            f4 = open(root+folder+"results_hbb.txt","w")
            f5 = open(root+folder+"results_htt.txt","w")
            f6 = open(root+folder+"results_m.txt","w")
            f7 = open(root+folder+"results_bhaa.txt","w")
            f8 = open(root+folder+"results_bhab.txt","w")
            f9 = open(root+folder+"results_bhaabb.txt","w")
            f10= open(root+folder+"results_bhbb.txt","w")
            f11= open(root+folder+"results_bhtt.txt","w")

            for i in xrange(bins):

                f1.write(str(k[i])+" "+str(np.mean(Pk_haa_matrix[:,i]))+" "+\
             str(np.std(Pk_haa_matrix[:,i]))+"\n")

                f2.write(str(k[i])+" "+str(np.mean(Pk_hab_matrix[:,i]))+" "+\
            str(np.std(Pk_hab_matrix[:,i]))+"\n")

                f3.write(str(k[i])+" "+str(np.mean(Pk_haabb_matrix[:,i]))+" "+\
             str(np.std(Pk_haabb_matrix[:,i]))+"\n")

                f4.write(str(k[i])+" "+str(np.mean(Pk_hbb_matrix[:,i]))+" "+\
             str(np.std(Pk_hbb_matrix[:,i]))+"\n")

                f5.write(str(k[i])+" "+str(np.mean(Pk_htt_matrix[:,i]))+" "+\
            str(np.std(Pk_htt_matrix[:,i]))+"\n")

                f6.write(str(k[i])+" "+str(np.mean(Pk_m_matrix[:,i]))+" "+\
             str(np.std(Pk_m_matrix[:,i]))+"\n")

                f7.write(str(k[i])+" "+str(np.mean(b_haa_matrix[:,i]))+" "+\
             str(np.std(b_haa_matrix[:,i]))+"\n")

                f8.write(str(k[i])+" "+str(np.mean(b_hab_matrix[:,i]))+" "+\
            str(np.std(b_hab_matrix[:,i]))+"\n")

                f9.write(str(k[i])+" "+str(np.mean(b_haabb_matrix[:,i]))+" "+\
             str(np.std(b_haabb_matrix[:,i]))+"\n")

                f10.write(str(k[i])+" "+str(np.mean(b_hbb_matrix[:,i]))+" "+\
            str(np.std(b_hbb_matrix[:,i]))+"\n")

                f11.write(str(k[i])+" "+str(np.mean(b_htt_matrix[:,i]))+" "+\
             str(np.std(b_htt_matrix[:,i]))+"\n")

            f1.close();  f2.close();  f3.close();  f4.close()
            f5.close();  f6.close();  f7.close();  f8.close()
            f9.close();  f10.close();  f11.close()














"""
    k, Pk_h_a = np.loadtxt(file_a,unpack=True)
    print'file %i read'%(i+1)
    Pk_h_above_matrix[i]=Pk_h_a

    f_a = open("results_halos_above.txt","w")
    for i in xrange(bins):
        f_a.write(str(k[i])+" "+str(np.mean(Pk_h_above_matrix[:,i]))+" "+str(np.std(Pk_h_above_matrix[:,i]))+"\n")
    f_a.close()
    #####################################

    ############ Pk below ###############
    file_b ='Pk_h_below_local%i.txt'%(i+1)
    k, Pk_h_b = np.loadtxt(file_b,unpack=True)
    print 'file %i read'%(i+1)
    Pk_h_below_matrix[i]=Pk_h_b

    f_b = open("results_halos_below.txt","w")
    for i in xrange(bins):
        f_b.write(str(k[i])+" "+str(np.mean(Pk_h_below_matrix[:,i]))+" "+str(np.std(Pk_h_below_matrix[:,i]))+"\n")
    f_b.close()
    #####################################

    
    file_t ='Pk_h_total_local%i.txt'%(i+1)
    k, Pk_h_t = np.loadtxt(file_a,unpack=True)
    print'file %i read'%(i+1)
    Pk_h_total_matrix[i]=Pk_h_t

    f_t = open("results_halos_total.txt","w")
    for i in xrange(bins):
        f_t.write(str(k[i])+" "+str(np.mean(Pk_h_total_matrix[:,i]))+" "+str(np.std(Pk_h_total_matrix[:,i]))+"\n")
    f_t.close()
    #####################################



    #file_bias ='b_h%i.txt'%(i+1)
    #k, b_h = np.loadtxt(file,unpack=True)
    print 'file %i read'%(i+1)
    file_pk ='/mnt/xfs1/home/tcourt/test/meanstd/Pk%i.txt'%(i+1)
    k, Pk = np.loadtxt(file_pk,unpack=True)
    #print 'file %i read'%(i+1)
    Pk_matrix[i]=Pk

    file_b ='Pk_h_below_local%i.txt'%(i+1)
    k, Pk_h_b = np.loadtxt(file_b,unpack=True)
    #print 'file %i read'%(i+1)
    Pk_h_below_matrix[i]=Pk_h_b

    file_a ='Pk_h_above_local%i.txt'%(i+1)
    k, Pk_h_a = np.loadtxt(file_a,unpack=True)
    #print 'file %i read'%(i+1)
    Pk_h_above_matrix[i]=Pk_h_a

    b_h_matrixa[i]=np.sqrt(Pk / Pk_h_a)
    b_h_matrixb[i]=np.sqrt(Pk / Pk_h_b)
        
    f_bias = open("results_bias_above.txt","w")
    for i in xrange(bins):
        f_bias.write(str(k[i])+" "+str(np.mean(b_h_matrixa[:,i]))+" "+str(np.std(b_h_matrixa[:,i]))+"\n")
    f_bias.close()

    f_bias2 = open("results_bias_below.txt","w")
    for i in xrange(bins):
        f_bias2.write(str(k[i])+" "+str(np.mean(b_h_matrixb[:,i]))+" "+str(np.std(b_h_matrixb[:,i]))+"\n")
    f_bias2.close()
"""


                

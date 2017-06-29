##################################################################
# This program takes the 100 files made by powerspectrum.py and  #
# it finds the mean and standard deviation of each of them then  #
# puts the results into a .txt file for plotting in GNUplot      #
#                                                                #
##################################################################

import numpy as np
#import matplotlib.pyplot as plt
import sys, os
resolution = 1024
############### setting up the array ##############################

root = '/mnt/xfs1/home/tcourt/test/'
folder = 'results_%i/'%resolution

k, Pk = np.loadtxt(root+folder+'Pk_h_total_local1.txt',unpack=True)

bins = len(Pk)
Pk_matrix         = np.zeros((100,bins),dtype=np.float32)
Pk_h_above_matrix = np.zeros((100,bins),dtype=np.float32)
Pk_h_below_matrix = np.zeros((100,bins),dtype=np.float32)
Pk_h_total_matrix = np.zeros((100,bins),dtype=np.float32)
Pk_m_total_matrix = np.zeros((100,bins),dtype=np.float32)
b_h_matrixa       = np.zeros((100,bins),dtype=np.float32)
b_h_matrixb       = np.zeros((100,bins),dtype=np.float32)
b_h_matrixt       = np.zeros((100,bins),dtype=np.float32)


    
######################### Halos Above #################################

for i in xrange(100):

    ##### find file names #####
    file_a = root+folder+'Pk_h_above_local%i.txt'%(i+1)
    file_b = root+folder+'Pk_h_below_local%i.txt'%(i+1)
    file_t = root+folder+'Pk_h_total_local%i.txt'%(i+1)
    file_m = root+'meanstd/Pk%i.txt'%(i+1)

    #### read files ####
    k, Pk_h_a = np.loadtxt(file_a,unpack=True)
    k, Pk_h_b = np.loadtxt(file_b,unpack=True)
    k, Pk_h_t = np.loadtxt(file_t,unpack=True)
    k, Pk_m   = np.loadtxt(file_m,unpack=True)

    #### fill big matrix ####
    Pk_h_above_matrix[i] = Pk_h_a
    Pk_h_below_matrix[i] = Pk_h_b
    Pk_h_total_matrix[i] = Pk_h_t
    Pk_m_total_matrix[i] = Pk_m
    b_h_matrixa[i]       = np.sqrt(Pk_h_a/Pk_m)
    b_h_matrixb[i]       = np.sqrt(Pk_h_b/Pk_m)
    b_h_matrixt[i]       = np.sqrt(Pk_h_t/Pk_m)


# write results to file
f1 = open(root+folder+"results_halos_above.txt","w")
f2 = open(root+folder+"results_halos_below.txt","w")
f3 = open(root+folder+"results_halos_total.txt","w")
f4 = open(root+folder+"results_matter.txt","w")
f5 = open(root+folder+"results_bias_halos_above.txt","w")
f6 = open(root+folder+"results_bias_halos_below.txt","w")
f7 = open(root+folder+"results_bias_halos_total.txt","w")
for i in xrange(bins):

    f1.write(str(k[i])+" "+str(np.mean(Pk_h_above_matrix[:,i]))+" "+\
             str(np.std(Pk_h_above_matrix[:,i]))+"\n")

    f2.write(str(k[i])+" "+str(np.mean(Pk_h_below_matrix[:,i]))+" "+\
            str(np.std(Pk_h_below_matrix[:,i]))+"\n")

    f3.write(str(k[i])+" "+str(np.mean(Pk_h_total_matrix[:,i]))+" "+\
             str(np.std(Pk_h_total_matrix[:,i]))+"\n")

    f4.write(str(k[i])+" "+str(np.mean(Pk_m_total_matrix[:,i]))+" "+\
             str(np.std(Pk_m_total_matrix[:,i]))+"\n")

    f5.write(str(k[i])+" "+str(np.mean(b_h_matrixa[:,i]))+" "+\
             str(np.std(b_h_matrixa[:,i]))+"\n")

    f6.write(str(k[i])+" "+str(np.mean(b_h_matrixb[:,i]))+" "+\
             str(np.std(b_h_matrixb[:,i]))+"\n")

    f7.write(str(k[i])+" "+str(np.mean(b_h_matrixt[:,i]))+" "+\
             str(np.std(b_h_matrixt[:,i]))+"\n")


f1.close();  f2.close();  f3.close();  f4.close()
f5.close();  f6.close();  f7.close()














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


                

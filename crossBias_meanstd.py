##################################################################
# This program takes the 100 files made by 100split.py and 0     #
# it finds the mean and standard deviation of each of them then  #
# puts the results into a .txt file for plotting in GNUplot      #
#                                                                #
##################################################################
import numpy as np
#import matplotlib.pyplot as plt
import sys, os
############### setting up the array #############################

root = '/mnt/xfs1/home/tcourt/test/'

for root_folder in ['0.0eV/','0.15eV/','0.6eV/']:
    print root_folder

    for dims1 in [32,64,128,256]:
        print dims1
        folder = 'CDM_bias_%s_%i/'%(root_folder[:-1],dims1)
    
        k, Pk = np.loadtxt(root+folder+'Pk_htt_1.txt',unpack=True)

        bins = len(Pk)
        #for names in ['cc','haa','haa-c','haa+hbb+2hab','hab','hbb','hbb-c','htt','htt-c']:
            #'Pk_'+ names + '_matrix' = np.zeros((100,bins),dtype=np.float32)
        Pk_cc_matrix    = np.zeros((100,bins),dtype=np.float32)
        Pk_haa_matrix   = np.zeros((100,bins),dtype=np.float32)
        Pk_haa_c_matrix = np.zeros((100,bins),dtype=np.float32)
        Pk_haabb_matrix = np.zeros((100,bins),dtype=np.float32)
        Pk_hab_matrix   = np.zeros((100,bins),dtype=np.float32)
        Pk_hbb_matrix   = np.zeros((100,bins),dtype=np.float32)
        Pk_hbb_c_matrix = np.zeros((100,bins),dtype=np.float32)
        Pk_htt_matrix   = np.zeros((100,bins),dtype=np.float32)
        Pk_htt_c_matrix = np.zeros((100,bins),dtype=np.float32)
        bc_cc_matrix    = np.zeros((100,bins),dtype=np.float32)
        bc_haa_matrix   = np.zeros((100,bins),dtype=np.float32)
        bc_haa_c_matrix = np.zeros((100,bins),dtype=np.float32)
        bc_haabb_matrix = np.zeros((100,bins),dtype=np.float32)
        bc_hab_matrix   = np.zeros((100,bins),dtype=np.float32)
        bc_hbb_matrix   = np.zeros((100,bins),dtype=np.float32)
        bc_hbb_c_matrix = np.zeros((100,bins),dtype=np.float32)
        bc_htt_matrix   = np.zeros((100,bins),dtype=np.float32)
        bc_htt_c_matrix = np.zeros((100,bins),dtype=np.float32)
   
    ######################### Halos Above #################################
        for i in xrange(100):

            ##### find file names #####
            file_cc    = root+folder+'Pk_cc_%i.txt'%(i+1)
            file_haa   = root+folder+'Pk_haa_%i.txt'%(i+1)
            file_haa_c = root+folder+'Pk_haa-c_%i.txt'%(i+1)
            file_haabb = root+folder+'Pk_haa+hbb+2hab_%i.txt'%(i+1)
            file_hab   = root+folder+'Pk_hab_%i.txt'%(i+1)
            file_hbb   = root+folder+'Pk_hbb_%i.txt'%(i+1)
            file_hbb_c = root+folder+'Pk_hbb-c_%i.txt'%(i+1)
            file_htt   = root+folder+'Pk_htt_%i.txt'%(i+1)
            file_htt_c = root+folder+'Pk_htt-c_%i.txt'%(i+1)

            #### read files ####
            k, Pk_cc    = np.loadtxt(file_cc,unpack=True)
            k, Pk_haa   = np.loadtxt(file_haa,unpack=True)   
            k, Pk_haa_c = np.loadtxt(file_haa_c,unpack=True) 
            k, Pk_haabb = np.loadtxt(file_haabb,unpack=True)
            k, Pk_hab   = np.loadtxt(file_hab,unpack=True)   
            k, Pk_hbb   = np.loadtxt(file_hbb,unpack=True)   
            k, Pk_hbb_c = np.loadtxt(file_hbb_c,unpack=True) 
            k, Pk_htt   = np.loadtxt(file_htt,unpack=True)   
            k, Pk_htt_c = np.loadtxt(file_htt_c,unpack=True) 
            
            #### fill big matrix ####
            Pk_cc_matrix[i]    = Pk_cc
            Pk_haa_matrix[i]   = Pk_haa
            Pk_haa_c_matrix[i] = Pk_haa_c
            Pk_haabb_matrix[i] = Pk_haabb
            Pk_hab_matrix[i]   = Pk_hab
            Pk_hbb_matrix[i]   = Pk_hbb
            Pk_hbb_c_matrix[i] = Pk_hbb_c
            Pk_htt_matrix[i]   = Pk_htt
            Pk_htt_c_matrix[i] = Pk_htt_c
            bc_haa_matrix[i]   = np.sqrt(Pk_haa/Pk_cc)
            bc_haa_c_matrix[i] = Pk_haa_c/Pk_cc
            bc_haabb_matrix[i] = np.sqrt(Pk_haabb/Pk_cc)
            bc_hab_matrix[i]   = np.sqrt(Pk_hab/Pk_cc)
            bc_hbb_matrix[i]   = np.sqrt(Pk_hbb/Pk_cc)
            bc_hbb_c_matrix[i] = Pk_hbb_c/Pk_cc
            bc_htt_matrix[i]   = np.sqrt(Pk_htt/Pk_cc)
            bc_htt_c_matrix[i] = Pk_htt_c/Pk_cc

           
            # write results to file
            f1 = open(root+folder+'results_haa_c.txt','w')
            f2 = open(root+folder+'results_hbb_c.txt','w')
            f3 = open(root+folder+'results_htt_c.txt','w')
            f4 = open(root+folder+'results_bhaa_c.txt','w')
            f5 = open(root+folder+'results_bhbb_c.txt','w')
            f6 = open(root+folder+'results_bhtt_c.txt','w')
            f7 = open(root+folder+'results_hcc.txt','w')


            for i in xrange(bins):

                f1.write(str(k[i])+" "+str(np.mean(Pk_haa_c_matrix[:,i]))+" "+\
             str(np.std(Pk_haa_c_matrix[:,i]))+"\n")

                f2.write(str(k[i])+" "+str(np.mean(Pk_hbb_c_matrix[:,i]))+" "+\
             str(np.std(Pk_hbb_c_matrix[:,i]))+"\n")

                f3.write(str(k[i])+" "+str(np.mean(Pk_htt_c_matrix[:,i]))+" "+\
            str(np.std(Pk_htt_c_matrix[:,i]))+"\n")

                f4.write(str(k[i])+" "+str(np.mean(bc_haa_c_matrix[:,i]))+" "+\
             str(np.std(bc_haa_c_matrix[:,i]))+"\n")

                f5.write(str(k[i])+" "+str(np.mean(bc_hbb_c_matrix[:,i]))+" "+\
            str(np.std(bc_hbb_c_matrix[:,i]))+"\n")

                f6.write(str(k[i])+" "+str(np.mean(bc_htt_matrix[:,i]))+" "+\
             str(np.std(bc_htt_matrix[:,i]))+"\n")

                f7.write(str(k[i])+" "+str(np.mean(Pk_cc_matrix[:,i]))+" "+\
             str(np.std(Pk_cc_matrix[:,i]))+"\n")

            f1.close();  f2.close();  f3.close();  f4.close()
            f5.close();  f6.close();  f7.close()





''';  f8.close()
            f9.close();  f10.close();  f11.close()

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
'''

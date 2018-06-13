#!/usr/bin/python
# -*- coding: utf-8 -*-

# This program use cross-correlation coefficient to select high quality RF
# for the flowing stacking
# For each station, we only keep those RFs, with cross-correlation coefficient
# great than cce and numbers of greater than cce not less rat% of total RFs

from obspy import read
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

import errno
import os
import glob
import shutil


#----------------------------
# Begin main program
#----------------------------



if __name__=='__main__':

    curdir = os.getcwd()
    
    # cce define the cross-correlation coefficient between two RF
    # the greater of this value means the more familiar of the two RF
    cce=0.8
    
    # the percent of total number of all RF, this ratio used to select
    # high quality RF for stacking
    rat=0.20
    
    print(curdir)
    # Remove previous low quality RF
    if os.path.exists("LowQRF"):
        shutil.rmtree("LowQRF")
   
    if not os.path.exists("LowQRF"):
        try:
           os.makedirs("LowQRF")
        except OSError as exception:
           if exception.errno != errno.EEXIST:
              raise
    # here we plot the correlation matrix of all the RF
    b=[]
    for f0 in glob.glob("*.1.5"):
        st0=read(f0)
        b.append(st0[0].data)

    c=np.array(b)
    try:
        cm=np.corrcoef(c)
        plt.imshow(cm,interpolation='nearest')
        plt.colorbar()
        plt.savefig("cc.png")
        plt.close()
    except ValueError:
        print("Something wrong when plot RF cross-correlation matrix!")
        pass
        
    # Get the total number of RF    
    fnum=len(glob.glob("*.1.5"))
    
    # empty array used to keep the low quality RF file name
    bdrf=[]
    
    for f1 in glob.glob("*.1.5"):
        i=0
        for f2 in glob.glob("*.1.5"):
            if f1!=f2:
                #print(f1,f2)
                st1=read(f1)
                st2=read(f2)
                cm=np.corrcoef(st1[0].data,st2[0].data)
                if cm[0,1] > cce:
                    i=i+1
        
        # For low quality RF, keep it's file name into a array
        if i <= int(fnum*rat):
            bdrf.append(f1)

    # finally, we move the low quality RF into LowQRF folder
    for f in bdrf:
        cmd="mv "+f+" ./LowQRF"
        print(cmd)
        os.system(cmd)        
    
    print("--------Finishi data process!-----")


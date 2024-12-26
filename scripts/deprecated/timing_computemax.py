#!/usr/bin/python3

from sys import argv #command line arguments
import numpy as np

inputfile=str(argv[1])

realtime,workime,memusage,output=np.loadtxt(inputfile, skiprows=1, unpack=True)#, usecols=([0,1,2,3]))

maxr=np.max(realtime)
maxw=np.max(workime)
maxm=np.max(memusage)
maxo=np.max(output)


print ('{0:16}{1:16.2f}{2:16.2f}{3:16}'.format(maxr,maxw,maxm/1024.0,maxo))

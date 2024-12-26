#!/usr/bin/python2

import numpy as np
import warnings
import matplotlib.pyplot as plt

import evaluation_module as em
from evaluation_module import g
from evaluation_module import ph






#===================
def main():
#===================

    global g

    em.get_output_info()
    em.read_info_files()
    x, y, z, clumpid = em.read_particle_data()


    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    ax.scatter(x, y, c='k', s=0.1)
    plt.show()




#==============================
if __name__ == "__main__":
#==============================
    main()



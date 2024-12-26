#!/usr/bin/python2


import numpy as np
import matplotlib.pyplot as plt
import evaluation_module as em
from evaluation_module import g
from evaluation_module import ph




#===================
def main():
#===================
    import matplotlib.colors as colors

    global g

    nc = 2000

    em.get_output_info()
    em.read_info_files()
    x, y, z, clumpid = em.read_particle_data()

    xind = np.digitize(x, bins=np.linspace(0,1,nc))
    yind = np.digitize(y, bins=np.linspace(0,1,nc))
    zind = np.digitize(z, bins=np.linspace(0,1,nc))

    # Note that I histogram (y, x) not (x, y) for the plot:
    # Imshow apparently needs y-axis to be first index
    # other option: histogram(x, y), then use np.swapaxes()
    particle_density, xedges, yedges = np.histogram2d(y, x, bins=nc, range=((0,1),(0,1)))
    particle_density[particle_density==0] = 0.1


    fig = plt.figure(figsize=(10,10), dpi=300)
    ax = fig.add_subplot(111)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.imshow(particle_density,
        origin='lower',
        extent=(0,1,0,1),
        norm=colors.LogNorm(),
        interpolation='bicubic',
        cmap='inferno')



    em.read_data()
    ax.scatter(g.galaxy_pos[:,0], g.galaxy_pos[:,1], c='white',s=1, lw=0)






    plt.savefig('galaxy_plot.png', dpi=300)




#==============================
if __name__ == "__main__":
#==============================
    main()


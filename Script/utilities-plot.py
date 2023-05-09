# -*- coding: utf-8 -*-
"""
Created on May 8 2023
by Kazuko Fuchi

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import rc
rc('text', usetex=True)
font = {'family' : 'normal',
        'size'   : 16}
plt.rc('font', **font)

def plot_fields(X,Y,vals,figname,dirname='',vmin=None,vmax=None,fig_resolution=200):
    fig = plt.figure()
    if vmin:
        plt.contourf(X,Y,vals,cmap='jet',vmin=vmin,vmax=vmax)
    else:
        plt.contourf(X,Y,vals,cmap='jet')
    plt.xlabel('x')
    plt.ylabel('y')
    fig.savefig(dirname+figname,dpi=fig_resoluion)

    return


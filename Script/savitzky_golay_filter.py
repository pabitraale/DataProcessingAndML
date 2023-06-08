import vtk
import numpy as np
import matplotlib.pyplot as plt
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from matplotlib import cm
from plotly import graph_objects as go
import numpy as np
import scipy
from scipy import signal
import savitzkygolay
from mpl_toolkits.mplot3d import Axes3D


def savitzky_filter_smooth(x, y,dataWithNoise,orginalData):
    data_fit = savitzkygolay.SGolayFilter2(window_size = 13, poly_order=3)(dataWithNoise)
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    ax.plot_wireframe(x, y, data_fit, linewidths=0.8, color='black', cstride = 10, rstride = 10)
    #ax.scatter(x, y, dataWithNoise, s=10, c='r')

    ax.plot_surface(x, y, data_fit, linewidth=0, cmap='cool')
    ax.plot_surface(x, y, orginalData, color ='green', linewidth=0, alpha =.40)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.show()




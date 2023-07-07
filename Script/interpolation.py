

import vtk
import numpy as np
import matplotlib.pyplot as plt
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtkmodules.vtkCommonColor import vtkNamedColors
#from mayavi import mlab
from scipy import interpolate
from matplotlib import cm
from scipy import signal
import savitzky_golay_filter
from mpl_toolkits.mplot3d import Axes3D


#open file and read unstructured grid
reader = vtk.vtkDataSetReader()
#reader.SetFileName("DataProcessingAndML\Data\cavity-coarse_100.vtk")
reader.SetFileName("DataProcessingAndML\Data\cavity_100[9331].vtk")
#reader.SetFileName("../Data/cavity_100[9331].vtk")
reader.ReadAllVectorsOn()
reader.ReadAllColorScalarsOn()
reader.Update()
grid = reader.GetOutput()

row, col = 20, 20

#get centroid of cell
cellCentersFilter = vtk.vtkCellCenters()
cellCentersFilter.SetInputData(grid)
cellCentersFilter.VertexCellsOn()
cellCentersFilter.Update()

print("\ncell center: ", cellCentersFilter.GetOutput().GetNumberOfPoints())
numberOfCenterInCells = cellCentersFilter.GetOutput().GetNumberOfPoints()

centroid = []
for i in range(0, numberOfCenterInCells ):
    centroid.append(cellCentersFilter.GetOutput().GetPoint(i)) #getting the centroid from each cell

centroidArr = np.array(centroid)
print("----cell centroids----")
print( centroidArr)

#get centroid x, y and z
centroidX = centroidArr[:, 0]
centroidY = centroidArr[:, 1]
centroidZ = centroidArr[:, 2]

cell_Data = grid.GetCellData()
cellData_p = cell_Data.GetArray('p') #get cell data value p
Cell_data_p = vtk_to_numpy(cellData_p)

cellData = cell_Data.GetArray('U')#get cell data value U
cell_data_u = vtk_to_numpy(cellData)

#get Cell U vector x, y and z
cell_X = cell_data_u[:, 0] 
cell_Y = cell_data_u[:, 1]
cell_Z = cell_data_u[:, 2]

#function to generate noise
def randomNoise(data):
    minValue = min(data)* 0.1
    maxValue = max(data)*0.1
    store_noise = []
    for i in range(0, data.size):
        rand_noise = np.random.uniform(minValue, maxValue) + data[i]
        store_noise.append(rand_noise)
    noise = np.array(store_noise)
    return noise

#get random noise for cell data p and U(x, y, z)
noise_p = randomNoise(Cell_data_p)
noise_X = randomNoise(cell_X)
noise_Y = randomNoise(cell_Y)
noise_Z = randomNoise(cell_Z)

#reshape centroid and p (cell) value to 20*20
centX = centroidX.reshape(row,col)
centY = centroidY.reshape(row,col)
centZ = centroidZ.reshape(row,col)

#reshaping for p, U(x, y, and z)
p = Cell_data_p.reshape(row,col)
Xcell = cell_X.reshape(row,col)
Ycell = cell_Y.reshape(row,col)
Zcell = cell_Z.reshape(row,col)

#reshape noise data
noiseP = noise_p.reshape(row, col)
noiseX = noise_X.reshape(row, col)
noiseY = noise_Y.reshape(row, col)
noiseZ = noise_Z.reshape(row, col)


#plot figures
def plot_fig(x, y, data, plotdir='../Data/',figname=''):
    z = data.reshape(row,col)
    fig, ax= plt.subplots(1,1)
    cp = ax.contourf(x,y, z)
    fig.colorbar(cp)
    ax.set_title('contours plot')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()
    #fig.savefig(plotdir+figname+'.png',dpi=200)

#plot figure for cell p value
plot_fig(centX, centY, Cell_data_p)
plot_fig(centX, centY, cell_X)
plot_fig(centX, centY, cell_Y)
plot_fig(centX, centY, cell_Z)
#plot figure for noise data p
plot_fig(centX, centY, noise_p, figname='p_w_noise')
plot_fig(centX, centY, noise_X, figname='ux_w_noise')
plot_fig(centX, centY, noise_Y, figname='uy_w_noise')
plot_fig(centX, centY, noise_Z, figname='uz_w_noise')

def smooth_interpolate(cenX, cenY, values, numberOfPoints, plotdir='../Data/',figname='sample'):
    centX_data = cenX[0, :]
    centY_data = cenY[:, 0]

    step = (max(centX_data)-(min(centX_data)))/numberOfPoints

    centDataX = np.arange(min(centX_data),  max(centX_data), step)
    centDataY = np.arange(min(centY_data),  max(centY_data), step)

    newX, newY = np.meshgrid(centDataX, centDataY)
    dataP = interpolate.bisplrep(cenY, cenX, values, s = 1)
    newP =  interpolate.bisplev(centDataX, centDataY, dataP)
    return newX, newY, newP

    # EDIT: Return the filtered data or interpolator from the function and do the plotting outside of the function. Do the plotting using plot_fig function.

def plt_fig_interpolate(newX, newY, newP, cenX, cenY, values, plotdir='../Data/',figname=''):
    fig = plt.figure()
    #lims = dict(cmap='RdBu_r', vmin = min(centX_data), vmax=max(centX_data))
    plt.pcolormesh(newX, newY, newP)
    plt.colorbar()
    plt.title("Interpolated function")
    plt.show()
    #fig.savefig(plotdir+figname+'_filtered'+'.png',dpi=200)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('without filter')
    surf = ax.plot_surface(cenX, cenY, values, cmap=cm.cool)
    plt.show()
    #fig.savefig(plotdir+figname+'_surf_data.png',dpi=200)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('filtered')
    surf = ax.plot_surface(newX, newY, newP, cmap=cm.cool)
    plt.show()
    #fig.savefig(plotdir+figname+'_surf_filtered'+'.png',dpi=200)

#smoothing data using interpolation 
totalDataPoints = 100

interplote_data_p = smooth_interpolate(centX, centY, p, totalDataPoints, figname='p')
plt_fig_interpolate(interplote_data_p[0], interplote_data_p[1], interplote_data_p[2], centX, centY, p)

interplote_data_x = smooth_interpolate(centX, centY, Xcell, totalDataPoints, figname='ux')
plt_fig_interpolate(interplote_data_x[0], interplote_data_x[1], interplote_data_x[2], centX, centY, Xcell)

interplote_data_y = smooth_interpolate(centX, centY, Ycell, totalDataPoints, figname='uy')
plt_fig_interpolate(interplote_data_y[0], interplote_data_y[1], interplote_data_y[2], centX, centY, Ycell)

interplote_data_z = smooth_interpolate(centX, centY, Zcell, totalDataPoints, figname='uz')
plt_fig_interpolate(interplote_data_z[0], interplote_data_z[1], interplote_data_z[2], centX, centY, Zcell)

#savitzky noise smoothing filter

filtered_data_p = savitzky_golay_filter.savitzky_filter_smooth(centX, centY,noiseP, p)
savitzky_golay_filter.plt_fig_savitzky(centX, centY, filtered_data_p, p)

filtered_data_X = savitzky_golay_filter.savitzky_filter_smooth(centX, centY,noiseX, Xcell)
savitzky_golay_filter.plt_fig_savitzky(centX, centY,filtered_data_X, Xcell)

filtered_data_Y = savitzky_golay_filter.savitzky_filter_smooth(centX, centY,noiseY, Ycell)
savitzky_golay_filter.plt_fig_savitzky(centX, centY,filtered_data_Y, Ycell)

filtered_data_Z = savitzky_golay_filter.savitzky_filter_smooth(centX, centY,noiseZ, Zcell)
savitzky_golay_filter.plt_fig_savitzky(centX, centY,filtered_data_Z, Zcell)

# EDIT: Make plots of the savitzky-golay filtred data. Save figures using different file names.


import vtk
import numpy as np
import matplotlib.pyplot as plt
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtkmodules.vtkCommonColor import vtkNamedColors

from utilities_plot import *

def random_noise(nx,ny,vmax,vmin,frac=0.1):
    noise = frac*(vmax-vmin)*np.random.rand(nx,ny)
    return noise 

def sin_noise(X,Y,dx,nx,ny,vmax,vmin,frac=0.1):
    r = frac*(vmax-vmin)*np.random.rand(nx,ny)
    s = np.sin(X/dx) 
    noise = s*r
    return noise 

#open file and read unstructured grid
reader = vtk.vtkUnstructuredGridReader()
#reader.SetFileName("DataProcessingAndML\Data\cavity-coarse_100.vtk")
#reader.SetFileName("../Data/cavity-coarse_100.vtk") # for Linux
#nx, ny = 2, 2
reader.SetFileName("../Data/cavity_100[9331].vtk") # for Linux
nx, ny = 20, 20
reader.Update()
grid = reader.GetOutput()

#get number of cells
numberOfCellsInGrid = grid.GetNumberOfCells()
print("number of cell: ", numberOfCellsInGrid)

##cell data 
cell_data = grid.GetCellData()
cell_array = cell_data.GetArray('cellID')

cell_data_p = cell_data.GetArray('p')
cell_data_u = cell_data.GetArray('U')

#get centroid of cell
cellCentersFilter = vtk.vtkCellCenters()
cellCentersFilter.SetInputData(grid)
cellCentersFilter.VertexCellsOn()
cellCentersFilter.Update()

numberOfCenterInCells = cellCentersFilter.GetOutput().GetNumberOfPoints()

centroid = []
for i in range(0, numberOfCenterInCells ):
    centroid.append(cellCentersFilter.GetOutput().GetPoint(i))
centroidArr = np.array(centroid)

#Plot fields
xc = centroidArr[:,0]
yc = centroidArr[:,1]
Xc = xc.reshape(nx,ny) # cell center x-coords
Yc = yc.reshape(ny,ny) # cell center y-coords

p = vtk_to_numpy(cell_data_p)
u1 = vtk_to_numpy(cell_data_u)[:,0]
u2 = vtk_to_numpy(cell_data_u)[:,1]
u3 = vtk_to_numpy(cell_data_u)[:,2]

vals_list = {'p':p,'u1':u1,'u2':u2,'u3':u3}
fieldnames = [val for val in vals_list.keys()]
dirname = '../Data/'
for i in range(len(fieldnames)):
    fieldname = fieldnames[i]
    vals = vals_list[fieldname]
    V = vals.reshape(nx,ny)
    vmin, vmax = min(vals), max(vals) 

    # Original field
    figname = 'field_'+fieldname
    plot_fields(Xc,Yc,V,figname,dirname=dirname,vmin=vmin,vmax=vmax)

    # Field + random noise 
    frac = 0.1 # noise scailng fraction
    noise_rand = random_noise(nx,ny,vmax,vmin,frac=frac)
    V_noise1 = V + noise_rand 
    figname = 'field_random_noise_'+fieldname
    plot_fields(Xc,Yc,V_noise1,figname,dirname=dirname,vmin=vmin,vmax=vmax)

    # Field + sine-distributed random noise 
    figname = 'field_sin_noise_'+fieldname
    dx = max(xc)-min(xc)
    noise_sin = sin_noise(Xc,Yc,dx,nx,ny,vmax,vmin,frac=frac)
    V_noise2 = V + noise_sin
    plot_fields(Xc,Yc,V_noise2,figname,dirname=dirname,vmin=vmin,vmax=vmax)


import vtk
import numpy as np
import matplotlib.pyplot as plt
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtkmodules.vtkCommonColor import vtkNamedColors

from utilities_plot import *

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
    figname = 'field_'+fieldname
    vals = vals_list[fieldname]
    V = vals.reshape(nx,ny)
    vmin, vmax = min(vals), max(vals) 
    plot_fields(Xc,Yc,V,figname,dirname=dirname,vmin=vmin,vmax=vmax)

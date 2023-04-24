import vtk
import numpy as np
import matplotlib.pyplot as plt
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk


'''print("hello world")

a = np.arange(15).reshape(3,5)
print(a)

plt.plot([1,2,3,4])
plt.ylabel('some number')
plt.show()'''

#open file and read unstructured grid
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("Data\cavity-coarse_100.vtk")
reader.Update()
grid = reader.GetOutput()
print(grid)
#get number of cells
print("number of cell: ", grid.GetNumberOfCells())

#read points
vtk_points = grid.GetPoints()
print(vtk_points.GetData())
xyz3d = vtk_to_numpy(vtk_points.GetData())
print(xyz3d)
#print(xyz3d.size)

print("get cell arrys")
print("=============================")

print("==========================================")


#cell data 
print("printing cell data")
cell_data = grid.GetCellData()
cell_arrya = cell_data.GetArray('cellID')
print(vtk_to_numpy(cell_arrya))
print('\n')

cell_data_p = cell_data.GetArray('p')
print(vtk_to_numpy(cell_data_p))
print('\n')

cell_data_u = cell_data.GetArray('U')
print(vtk_to_numpy(cell_data_u))
print('\n')

#point data
print("printing point data")
point_data = grid.GetPointData()
points_data_p = point_data.GetArray('p')
print(vtk_to_numpy(points_data_p))

print('\n')
points_data_u = point_data.GetArray('U')
print(vtk_to_numpy(points_data_u))


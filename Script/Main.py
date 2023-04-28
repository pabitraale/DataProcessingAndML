import vtk
import numpy as np
import matplotlib.pyplot as plt
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk


#open file and read unstructured grid
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("DataProcessingAndML\Data\cavity-coarse_100.vtk")
reader.Update()
grid = reader.GetOutput()

#get number of cells
numberOfCellsInGrid = grid.GetNumberOfCells()
print("number of cell: ", numberOfCellsInGrid)

#read points
vtk_points = grid.GetPoints()
xyz3d = vtk_to_numpy(vtk_points.GetData())

print("number of points: ", grid.GetNumberOfPoints())
print('POINTS  %s  float'% (grid.GetNumberOfPoints()))
print(xyz3d)
print(xyz3d.shape)
print(type(xyz3d))
x = xyz3d[:, 0]
print(x)
print(x.shape)
'''i = 0
while(i < grid.GetNumberOfPoints()):
    print("Next point: %d, %s " % (i+1, grid.GetPoint(i)))
    i = i+1'''

#find cell id
verticesInSingleCell = []
# print index 0 to find out the total number of vertices in each cell
for i in range(0, numberOfCellsInGrid):
    currentCell = grid.GetCell(i)
    verticesInSingleCell.append(currentCell.GetNumberOfPoints())
print("\nTotal vetices in each cell")
print(verticesInSingleCell)

#get vertices of each cell and store in cellConnectivityPoints
#cellConnectivityPoints = []
cellVertices = []

for i in range ( 0, numberOfCellsInGrid):
    currentCell = grid.GetCell(i)
    #print(currentCell.GetNumberOfPoints())
    arr = [currentCell.GetNumberOfPoints()]
    print("array: ", arr) 
    for j in range (0, currentCell.GetNumberOfPoints()):
        tempPointStore = currentCell.GetPointId(j)
        #print(" temp: ", tempPointStore)
        #cellConnectivityPoints.append(tempPointStore )
        arr.append(tempPointStore)# = np.array(cellConnectivityPoints)
    
    cellVertices.append(arr)
    
#print(cellVertices[0])       
for i in range (0, numberOfCellsInGrid):
    print(cellVertices[i])
   

print('\nCELL %s ' % numberOfCellsInGrid)
print("get cell arrys")
print("=============================")
# get cell 0 and print it's types 
cell0 = grid.GetCell(0)
print("cell types: ", cell0.GetCellType()) #prints cell types
print("first point: ", cell0.GetNumberOfPoints()) #print first point
print("=============================")
#get total points of cells (CELL 4 36)

print("total points: ", grid.GetNumberOfPoints())



#create a list to store cell types
storeCellTypes = []
for i in range (0, numberOfCellsInGrid):
    storeCellTypes.append(grid.GetCellType(i))
print("cell types: ", storeCellTypes)

  
#cell data 
print("printing cell data")
print("cellID 1 4 int")
cell_data = grid.GetCellData()
cell_array = cell_data.GetArray('cellID')
print(vtk_to_numpy(cell_array))
print('\n')

print("p 1 float")
cell_data_p = cell_data.GetArray('p')
print(vtk_to_numpy(cell_data_p))
print('\n')

print("U 3 4 float")
cell_data_u = cell_data.GetArray('U')
print(vtk_to_numpy(cell_data_u))
print('\n')

#point data
print("printing point data")
point_data = grid.GetPointData()
points_data_p = point_data.GetArray('p')
print("p 1 18 float")
print(vtk_to_numpy(points_data_p))

print('\n')
points_data_u = point_data.GetArray('U')
print("U 3 18 float")
print(vtk_to_numpy(points_data_u))

#create the mapper that corresponds the objects of the vtk file
mapper = vtk.vtkDataSetMapper()
mapper.SetInputData(grid)

#Create the Actor
actor = vtk.vtkActor()
actor.GetProperty().SetRepresentationToWireframe()
actor.SetMapper(mapper)




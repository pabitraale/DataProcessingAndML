import vtk
import numpy as np
import matplotlib.pyplot as plt
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtkmodules.vtkCommonColor import vtkNamedColors

#open file and read unstructured grid
reader = vtk.vtkUnstructuredGridReader()
#reader.SetFileName("DataProcessingAndML\Data\cavity_100[9331].vtk")
reader.SetFileName("../Data/cavity_100[9331].vtk")
reader.Update()
grid = reader.GetOutput()
#print(reader.GetOutput()) #this print vtk unStructured grid6
#print(reader.GetOutputPort()) # this print vtk Algorithm output

#get number of cells
numberOfCellsInGrid = grid.GetNumberOfCells()
print("number of cell: ", numberOfCellsInGrid)

#read points
vtk_points = grid.GetPoints()
xyz3d = vtk_to_numpy(vtk_points.GetData())

print("number of points: ", grid.GetNumberOfPoints())
print('POINTS  %s  float'% (grid.GetNumberOfPoints()))
print(xyz3d)
print("x value: ", xyz3d[:, 0])


#get color form vtk
color = vtk.vtkNamedColors()
output_port = reader.GetOutputPort() #creating output port and assigning reader outputport
scalar_range = grid.GetScalarRange()

#find cell id
verticesInSingleCell = []
#Get the total number of vertices in each cell
for i in range(0, numberOfCellsInGrid):
    currentCell = grid.GetCell(i)
    verticesInSingleCell.append(currentCell.GetNumberOfPoints())
#print("\nVetices in each cell")
#print(verticesInSingleCell)

#get vertices of each cell and store in cellConnectivityPoints
cellVertices = []

for i in range ( 0, numberOfCellsInGrid):
    currentCell = grid.GetCell(i)
    arr = [currentCell.GetNumberOfPoints()]
    for j in range (0, currentCell.GetNumberOfPoints()):
        tempPointStore = currentCell.GetPointId(j)
        arr.append(tempPointStore)
    
    cellVertices.append(arr)
    
#print(cellVertices[0])  
# index 0 (value 8) represent 8 vertices in each cell
# after index 0, they represent number in each points
print("\nCells point values")    
for i in range (0, numberOfCellsInGrid):
    print(cellVertices[i])
   

print('\nCELL %s ' % numberOfCellsInGrid)
# get cell 0 and print it's types 
cell0 = grid.GetCell(0)
print("cell types: ", cell0.GetCellType()) #prints cell types
print("first point: ", cell0.GetNumberOfPoints()) #print first point

#get total points of cells (CELL 4 36)
print("total points: ", grid.GetNumberOfPoints())

#create a list to store cell types
storeCellTypes = []
for i in range (0, numberOfCellsInGrid):
    storeCellTypes.append(grid.GetCellType(i))
print("cell types: ", storeCellTypes)

  
#cell data 
print("\nprinting cell data")
print("cellID 1 4 int")
cell_data = grid.GetCellData()
cell_array = cell_data.GetArray('cellID')
print(vtk_to_numpy(cell_array))

print("\np 1 float")
cell_data_p = cell_data.GetArray('p')
print(vtk_to_numpy(cell_data_p))

print("\nU 3 4 float")
cell_data_u = cell_data.GetArray('U')
print(vtk_to_numpy(cell_data_u))


#point data
print("\nprinting point data")
point_data = grid.GetPointData()
points_data_p = point_data.GetArray('p')
print("p 1 18 float")
print(vtk_to_numpy(points_data_p))

points_data_u = point_data.GetArray('U')
print("\nU 3 18 float")
print(vtk_to_numpy(points_data_u))

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
print("centroid x: ", centroidArr[:, 0])



#Display cell centers
centerMapper = vtk.vtkDataSetMapper()
centerMapper.SetInputConnection(cellCentersFilter.GetOutputPort())
centerActor = vtk.vtkActor()
centerActor.SetMapper(centerMapper)
centerActor.GetProperty().SetPointSize(10)
centerActor.GetProperty().SetColor(color.GetColor3d('White'))



#create the mapper that corresponds the objects of the vtk file
# This creates grid
'''mapper = vtk.vtkDataSetMapper()
mapper.SetInputConnection(output_port)
mapper.SetScalarRange(scalar_range)

#Create the Actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetRepresentationToWireframe() 
#actor.GetProperty().SetColor(color.GetColor3d('tomato'))
actor.GetProperty().EdgeVisibilityOn()'''


#contour filter
contour = vtk.vtkContourGrid()
contour.SetInputConnection(output_port)

#store different data from cell data or point data
#getDataArray = grid.GetPointData().GetArray(0) # getting points data p
getDataArray = grid.GetPointData().GetArray(1) # getting points data vector U

getXMinMax = grid.GetPointData().GetArray(1).GetRange(0) #get min and max value of x
print("x min and max : ", getXMinMax)




'''getDataPointX = vtk_to_numpy(getDataArray)[:, 0]
print("x: ", getDataPointX)
print("length of x: ", getDataPointX.shape)

getDataPointY = vtk_to_numpy(getDataArray)[:, 1]
print("y: ", getDataPointY)
print("length of y: ", getDataPointY.shape)

getDataPointZ = vtk_to_numpy(getDataArray)[:, 2]
print("Z: ", getDataPointZ)
print("length of Z: ", getDataPointZ.shape)'''


#getDataArray = grid.GetCellData().GetArray(1) # getting cell data p
#getDataArray = grid.GetCellData().GetArray(2) # getting cell data vector U
'''print("getData Arrya cell: ", vtk_to_numpy(getDataArray))
getCellDataX = vtk_to_numpy(getDataArray)[:, 0] # x coordinate of U data vector
print("cell data x: ", getCellDataX)
print("cell data shape: ", getCellDataX.shape)
getCellDataY = vtk_to_numpy(getDataArray)[:, 1] # x coordinate of U data vector
print("cell data x: ", getCellDataY)
print("cell data shape: ", getCellDataY.shape)
getCellDataZ = vtk_to_numpy(getDataArray)[:, 2] # x coordinate of U data vector
print("cell data x: ", getCellDataZ)
print("cell data shape: ", getCellDataZ.shape)'''

#get min value  and max value
scalarRange = getDataArray.GetValueRange()

#generate 10 value between min value and max value
contour.GenerateValues(10, scalarRange[0], scalarRange[1])
contour.Update()

#creat look up table and set number in grid 
lut = vtk.vtkLookupTable()
lut.SetNumberOfColors(100)
lut.SetHueRange(0.655, 0.0)
lut.Build()

contourMapper = vtk.vtkDataSetMapper()
contourMapper.SetInputConnection(output_port)
#contourMapper.GetInput().GetPointData().SetActiveScalars('p') # for data point p
#contourMapper.GetInput().GetPointData().SetActiveScalars('U') # for data point vector U
#contourMapper.GetInput().GetCellData().SetActiveScalars('p')  # for cell data p
#contourMapper.GetInput().GetCellData().SetActiveScalars('U') # for cell data vector U


contourMapper.ScalarVisibilityOn()
contourMapper.SetScalarRange(scalarRange)
contourMapper.SetLookupTable(lut)
contourMapper.SetInterpolateScalarsBeforeMapping(1)

#set scalar bar
scalar_bar = vtk.vtkScalarBarActor()
scalar_bar.SetOrientationToHorizontal()
scalar_bar.SetLookupTable(lut)

#create grid actor
contourActor = vtk.vtkActor()
contourActor.SetMapper(contourMapper)
contourActor.GetProperty().EdgeVisibilityOff()


#Create the Renderer
contextView = vtk.vtkContextView()
renderer =contextView.GetRenderer()
renderer.AddActor(contourActor)
#renderer.AddActor(actor)
#renderer.AddActor(centerActor) #add center actor
renderer.SetBackground(color.GetColor3d('indigo'))

#Create the RendererWindow
renderer_window = vtk.vtkRenderWindow()
renderer_window.AddRenderer(renderer)
renderer_window.SetSize(700, 600)
renderer_window.SetWindowName("contour")

#Create the renderWindow interactor and siaplay the vtk file
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renderer_window)

#create scalar bar
scalarBarWidge = vtk.vtkScalarBarWidget()
scalarBarWidge.SetInteractor(interactor)
scalarBarWidge.SetScalarBarActor(scalar_bar)
scalarBarWidge.On()

#initializeing and starting window
interactor.Initialize()
interactor.Render()
interactor.Start()




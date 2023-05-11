import vtk
import numpy as np
import matplotlib.pyplot as plt
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtkmodules.vtkCommonColor import vtkNamedColors


#open file and read unstructured grid
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("DataProcessingAndML\Data\cavity-coarse_100.vtk")
reader.Update()
grid = reader.GetOutput()
#print("vtk unStructured gird :")
print(reader.GetOutput())
#print("vtk Algorithm output:")
print(reader.GetOutputPort())



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

#get color form vtk
color = vtk.vtkNamedColors()
output_port = reader.GetOutputPort()
scalar_range = grid.GetScalarRange()
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




#get centroid of cell
'''cellCentersFilter = vtk.vtkCellCenters()
cellCentersFilter.SetInputData(grid)
cellCentersFilter.VertexCellsOn()
cellCentersFilter.Update()

print("cell center: ", cellCentersFilter.GetOutput().GetNumberOfPoints())
numberOfCenterInCells = cellCentersFilter.GetOutput().GetNumberOfPoints()

centroid = []
for i in range(0, numberOfCenterInCells ):
    #print("points: ", cellCentersFilter.GetOutput().GetPoint(i))
    centroid.append(cellCentersFilter.GetOutput().GetPoint(i))

centroidArr = np.array(centroid)
print("----cell centroids----")
print( centroidArr)
print(centroidArr.shape)

#Display cell centers
centerMapper = vtk.vtkDataSetMapper()
centerMapper.SetInputConnection(cellCentersFilter.GetOutputPort())
centerActor = vtk.vtkActor()
centerActor.SetMapper(centerMapper)
centerActor.GetProperty().SetPointSize(10)
centerActor.GetProperty().SetColor(color.GetColor3d('White'))'''



#create the mapper that corresponds the objects of the vtk file
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
getDataArray = grid.GetPointData().GetArray(0) # getting points data p

#get min value  and max value
scalarRange = getDataArray.GetValueRange()

#generate 10 value between min value and max value
contour.GenerateValues(10, scalarRange[0], scalarRange[1])
contour.Update()


#create outline
'''outline = vtk.vtkOutlineFilter()
outline.SetInputConnection(reader.GetOutputPort())

outlineMapper = vtk.vtkDataSetMapper()
outlineMapper.SetInputConnection(outline.GetOutputPort())

outlineActor = vtk.vtkActor()
outlineActor.SetMapper(outlineMapper)
outlineActor.GetProperty().SetColor(color.GetColor3d("Red"))
outlineActor.GetProperty().SetLineWidth(3.0)'''

#creat look up table and set number in grid 
lut = vtk.vtkLookupTable()
lut.SetNumberOfColors(5)
lut.SetHueRange(0.655, 0.0)
lut.Build()

contourMapper = vtk.vtkDataSetMapper()
contourMapper.SetInputConnection(output_port)
#contourMapper.GetInput().GetPointData().SetActiveScalars('p') # for data point p
#contourMapper.GetInput().GetPointData().SetActiveScalars('U') # for data point vector U
contourMapper.GetInput().GetCellData().SetActiveScalars('p')  # for cell data p
#contourMapper.GetInput().GetCellData().SetActiveScalars('U') # for cell data vector U
contourMapper.ScalarVisibilityOn()
contourMapper.SetScalarRange(scalar_range)
contourMapper.SetLookupTable(lut)
contourMapper.SetInterpolateScalarsBeforeMapping(1)

#set scalar bar
scalar_bar = vtk.vtkScalarBarActor()
scalar_bar.SetOrientationToHorizontal()
scalar_bar.SetLookupTable(lut)



#Create the Renderer
renderer = vtk.vtkRenderer()
#renderer.AddActor(actor)
#renderer.AddActor(centerActor) #add center actor
#renderer.AddActor(contourActor) # add contour actor
#renderer.AddActor(outlineActor)
renderer.SetBackground(color.GetColor3d('white'))
renderer.ResetCamera()
renderer.GetActiveCamera().Azimuth(30)
renderer.GetActiveCamera().Elevation(30)
#renderer.SetLayer(0)


#Create the RendererWindow
renderer_window = vtk.vtkRenderWindow()
renderer_window.AddRenderer(renderer)
renderer_window.SetSize(500, 500)
renderer_window.SetWindowName("contour")

#Create the renderWindow interactor and siaplay the vtk file
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renderer_window)
interactor.Initialize()
interactor.Render()
interactor.Start()




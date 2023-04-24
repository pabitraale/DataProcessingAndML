import vtk
import numpy as np
import matplotlib.pyplot as plt


print("hello world")

a = np.arange(15).reshape(3,5)
print(a)

plt.plot([1,2,3,4])
plt.ylabel('some number')
plt.show()

reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName()
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

with open('scatter.txt') as f:
    l = f.readlines() 
    data = [v.split(' ') for v in l]
    data = np.asarray(data)

rank = data[:,0].astype(int)
cmap = plt.get_cmap('jet', np.max(rank))
rank = rank / np.max(rank)

print(rank)
points=data[:,1::].astype(float)

x,y,z = points[:,0],points[:,1],points[:,2]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c=rank, cmap=cmap)
plt.show()


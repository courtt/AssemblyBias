import numpy as np
import matplotlib.pyplot as plt

pos = np.random.random((50,3))
pos = pos*2 - 1
pos_a = np.where(pos[:,2]>0)[0]
pos_b = np.where(pos[:,2]<0)[0]
slicer = np.where(np.logical_and(pos[:,2]<0.3, pos[:,2]>-0.3))[0]
slice_a = np.where(pos[pos_a,2]<0.3)[0]
slice_b = np.where(pos[pos_b,2]>-0.3)[0]
mass = np.random.random((50,1))
mass = mass*100
print mass
print len(pos[slicer])
print len(pos[pos_a[slice_a]])
print len(pos[pos_b[slice_b]])

plt.scatter(pos[slicer,0],pos[slicer,1],s=mass[slicer])
#plt.show()
plt.savefig('01.png')

plt.clf()
plt.scatter(pos[pos_a[slice_a],0],pos[pos_a[slice_a],1],c='r',s=1*(mass[pos_a[slice_a]]/np.min(mass[pos_a[slice_a]]))+10,a=0.4)
plt.scatter(pos[pos_b[slice_b],0],pos[pos_b[slice_b],1],c='b',s=1*(mass[pos_b[slice_b]]/np.min(mass[pos_b[slice_b]]))+10,a=0.4)
#plt.show()
plt.savefig('02.png')

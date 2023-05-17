import numpy as np
from numpy import random
import matplotlib.pyplot as plt

x = []
y = []
v = []
r = []

e = 0.2
G = 4*(np.pi**2)
M = 1010

for i in range(12):
    j = 12+3*i
    angles = np.linspace(0,1,j)*(np.pi*2)
    x.extend(2*i*np.cos(angles))
    y.extend(2*i*np.sin(angles))
    r.extend(np.full(shape=j, fill_value=2.0*i))

r = np.array(r, float)
v = np.sqrt(G*M*r/(r**2+e**2))

plt.scatter(x,y)
plt.savefig("initial.jpg")

print(len(r))
sourceFile = open('initial.txt', 'w')
for i in range(len(r)):
    print('{:f}\t{:f}\t{:f}\t{:f}'.format(x[i], y[i], v[i], r[i]), file = sourceFile)
sourceFile.close()

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

n = 1
x = [[[] for i in range(3)] for j in range(n)]
t = []

f = open('leapint.out')
content = f.readlines()
f.close()
for line in content:
    s = line.split()
    i = int(s[0])
    if i == 0:
        t.append(float(s[1]))
    for j in range(3):
        x[i][j].append(float(s[2 + j]))

fig = plt.figure(figsize=(12, 10))
ax = plt.axes(projection='3d')
for i in range(n):
    ax.plot3D(x[i][0], x[i][1], x[i][2])
ax.plot3D(0, 0, 0, 'ko')
ax.set_xlim((-22, 2))
ax.set_ylim((-2, 30))
ax.set_zlim((-12, 2))
ax.set_xlabel('x (AU)')
ax.set_ylabel('y (AU)')
ax.set_zlabel('z (AU)')
plt.title('Halley Comet (1P/Halley) in 27000 days (~ a Halley year)')
plt.show()
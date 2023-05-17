import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

n = 9
x = [[[] for i in range(3)] for j in range(n)]
t = []
planet = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']

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
    ax.plot3D(x[i][0], x[i][1], x[i][2], label=planet[i])
ax.set_xlim((-50, 50))
ax.set_ylim((-50, 50))
ax.set_zlim((-50, 50))
ax.set_xlabel('x (AU)')
ax.set_ylabel('y (AU)')
ax.set_zlabel('z (AU)')
plt.title('Solar system in 100000 days (~ a pluto year)')
plt.legend()
plt.savefig('plot_p2.png')
# plt.show()

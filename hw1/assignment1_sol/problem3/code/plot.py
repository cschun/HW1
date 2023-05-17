import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

n = 9
x = [[[] for i in range(3)] for j in range(n)]
t = []
kepler_x = [[[] for i in range(3)] for j in range(n)]

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
print('ndata =', len(x[0][0]))

f = open('kepler.out')
content = f.readlines()
f.close()
for line in content:
    s = line.split()
    i = int(s[0])
    for j in range(3):
        kepler_x[i][j].append(float(s[1 + j]))

# plot
fig = plt.figure(figsize=(12, 10))
ax = plt.axes(projection='3d')
plt.title('Leapfrog Simulation vs Kepler Analytic Orbit: Mars (100 revolution)')
# plot leapfrog and kepler
i = 3
ax.plot3D(x[i][0], x[i][1], x[i][2], label="Leapfrog")
ax.plot3D(kepler_x[i][0], kepler_x[i][1], kepler_x[i][2], label="Kepler")
# plot the sun
ax.plot3D(0, 0, 0, 'ko')
#ax.set_xlim((-50, 50))
#ax.set_ylim((-50, 50))
#ax.set_zlim((-50, 50))
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
plt.legend()
plt.savefig('plot_p4b4.png')
#plt.show()
 
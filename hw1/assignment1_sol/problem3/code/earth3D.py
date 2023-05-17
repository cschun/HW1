import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
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
ndata = len(x[0][0])
print('ndata =', ndata)

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
ax = p3.Axes3D(fig)
plt.title('Earth')

# initialize scatters
scatter = ax.scatter([x[2][0][0]], [x[2][1][0]], [x[2][2][0]], c='m')

# set figure properties
ax.scatter(0, 0, 0, c='k')
ax.plot(kepler_x[2][0], kepler_x[2][1], kepler_x[2][2])
ax.set_xlim3d([-1.2, 1.2])
ax.set_xlabel('X (AU)')
ax.set_ylim3d([-1.2, 1.2])
ax.set_ylabel('Y (AU)')
ax.set_zlim3d([-0.0001, 0.0001])
ax.set_zlabel('Z (AU)')
#ax.view_init(30, 30)

def update(idata, x, scatter):
    scatter._offsets3d = ([x[2][0][idata]], [x[2][1][idata]], [x[2][2][idata]])
    return scatter

ani = animation.FuncAnimation(fig, update, frames=ndata, fargs=(x, scatter), interval=16, blit=False)
ani.save('plot_p3c.mp4')
#plt.show()
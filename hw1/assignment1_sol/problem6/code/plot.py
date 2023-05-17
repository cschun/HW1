import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
import numpy as np

n = 2
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
ndata = len(x[0][0])
print('ndata =', ndata)

# create figure and axes
fig = plt.figure(figsize=(12, 10))
ax = p3.Axes3D(fig)

# initialize scatters
scatters = [ax.scatter([x[i][0][0]], [x[i][1][0]], [x[i][2][0]]) for i in range(n)]

# set figure properties
ax.scatter(0, 0, 0, 'ko')
ax.set_xlim3d([-22, 2])
ax.set_ylim3d([-2, 30])
ax.set_zlim3d([-12, 2])
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
ax.view_init(30, 30)

def update(idata, x, t, scatters):
    #ax.set_title('Halley Comet (1P/Halley) in 27000 days (~ a Halley year)\nnow at {:.0f} years {:.0f} days'.format(t[idata]/365, t[idata]%365))
    for i in range(len(scatters)):
        scatters[i]._offsets3d = ([x[i][0][idata]], [x[i][1][idata]], [x[i][2][idata]])
    return scatters

ani = animation.FuncAnimation(fig, update, frames=ndata, fargs=(x, t, scatters), interval=10, blit=False)
ani.save('plot_p6.mp4')
ani.save('plot_p6.gif')
#plt.show()
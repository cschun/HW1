import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3

f = open('leapint_p1d.out')
content = f.readlines()
f.close()

nline = len(content)
print('read {:d} lines of data'.format(nline))

x = [[], [], []]
v = [[], [], []]
t = []
for line in content:
    [tt, ii, xx, vv] = line.split()
    tt = float(tt)
    ii = int(ii)
    xx = float(xx)
    vv = float(vv)
    if tt not in t:
        t.append(tt)
    x[ii].append(xx)
    v[ii].append(vv)

fig = plt.figure(figsize=(11,5))
fig.suptitle('Problem 1 (d)')
ax1 = plt.subplot(121)
ax1.set_xlabel('x')
ax1.set_ylabel('v')
ax1.plot(x[0], v[0], 'b.')
ax1.plot(x[1], v[1], 'r.')
ax1.plot(x[2], v[2], 'g.')
ax1.grid(which='both')
ax2 = plt.subplot(122, projection='3d')
ax2.set_xlabel('x')
ax2.set_ylabel('v')
ax2.set_zlabel('t')
ax2.scatter(x[0], v[0], t, 'b.')
ax2.scatter(x[1], v[1], t, 'r.')
ax2.scatter(x[2], v[2], t, 'g.')
fig.savefig('plot_p1d.png')
plt.close()

fig = plt.figure()
fig.suptitle('Problem 1 (d)')
ax1 = plt.subplot(111)
ax1.set_xlabel('x')
ax1.set_ylabel('v')
xdata, ydata = [], []
line, = ax1.plot([], [], 'b.')

def init():
    ax1.set_xlim(-4, 21)
    ax1.set_ylim(-2, 3.2)
    return line,

def animate(i):
    xdata = [x[_][i] for _ in range(3)]
    ydata = [v[_][i] for _ in range(3)]
    ax1.set_title('t = {:.4f}'.format(t[i]))
    line.set_data(xdata, ydata)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=65, interval=130, blit=True)
anim.save('plot_p1d.mp4')

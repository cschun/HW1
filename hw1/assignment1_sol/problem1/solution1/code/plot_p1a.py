import matplotlib.pyplot as plt

f = open('leapint_p1a.out')
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

fig = plt.figure(figsize=(11, 5))
fig.suptitle('Problem 1 (a)')
ax1 = fig.add_subplot(121)
ax1.set_xlabel('x')
ax1.set_ylabel('v')
ax1.plot(x[0], v[0], 'b.')
ax1.plot(x[1], v[1], 'r.')
ax1.plot(x[2], v[2], 'g.')
ax1.grid(which='both')
ax2 = fig.add_subplot(122, projection='3d')
ax2.set_xlabel('x')
ax2.set_ylabel('v')
ax2.set_zlabel('t')
ax2.scatter(x[0], v[0], t, 'b.')
ax2.scatter(x[1], v[1], t, 'r.')
ax2.scatter(x[2], v[2], t, 'g.')
plt.savefig('plot_p1a.png')
plt.close()


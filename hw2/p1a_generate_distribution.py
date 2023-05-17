import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 

# initial conditions
M = 2.0 * 10**11        # in unit of Msolar
R = 1500.0              # in unit of parsec
G = 4.49879 * 10**-3    # in unit of pc^3 / (Msolar * Myr^2)
N = 10**4               # number of points
R_MAX = 10 * R          # maximum radius considered
V_MAX = 1000            # maximum speed considered

# derived parameters
m = M / N               # mass of each point
def U(r):
    return - G * M * R**-1 * (1 + (r / R)**2)**-0.5
def E(r, v):
    return U(r) + 0.5 * v**2
def f(r, v):
    E_ = E(r, v)
    return 24 * 2**0.5 / (7 * np.pi**3) * G**-5 * M**-4 * R**2 * (- E_ * np.heaviside(-E_, 0))**3.5
def rho(r):
    return 3 / (4 * np.pi) * M * R**-3 * (1 + (r / R)**2)**-2.5

'''
# the phase space distribution plot
r = np.linspace(0, 2*R)
v = np.linspace(300, 600)
#r = np.linspace(0, 10 * R)
#v = np.linspace(0, 1000)
a, b = np.meshgrid(r, v)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(r, v, f(a, b) * a**2 * b**2, 100)
ax.set_xlabel('r')
ax.set_ylabel('v')
ax.set_zlabel('f * r^2 * v^2')
plt.show()
'''

# find a distribution
def random_direction_vector(r):
    phi = 2 * np.pi * np.random.rand()
    cos_theta = -1 + 2 * np.random.rand()
    sin_theta = (1 - cos_theta**2)**0.5
    x = r * sin_theta * np.cos(phi)
    y = r * sin_theta * np.sin(phi)
    z = r * cos_theta
    return x, y, z
fout = open('plummer_initial_condition.txt', 'w')
count_total = 0
count = 0
fr2v2_max = 1200
while count < N:
    r = R_MAX * np.random.rand()
    v = V_MAX * np.random.rand()
    fr2v2 = f(r, v) * r**2 * v**2
    threshold = fr2v2 / fr2v2_max
    assert(threshold < 1)
    if np.random.rand() < threshold:
        # found a point!
        count += 1
        x, y, z = random_direction_vector(r)
        vx, vy, vz = random_direction_vector(v)
        fout.write('{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n'.format(x, y, z, vx, vy, vz))
    # show progress
    count_total += 1
    if count % 100 == 0:
        print('found {:d} points, total attempts: {:d}'.format(count, count_total))
fout.close()

    
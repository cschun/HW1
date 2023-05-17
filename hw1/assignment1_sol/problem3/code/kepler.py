import numpy as np
import math
from scipy.spatial.transform import Rotation as R

n = 9 # number of planets
x = np.zeros((3, n))
v = np.zeros((3, n))
# in units of AU, day, solar mass
G = 2.959159 * 10**(-4) 
npoint = 100 # number of points each orbit

f = open('initial_condition.txt')
content = f.readlines()
f.close()
for i in range(n):
    line = content[i].split()
    for j in range(3):
        x[j, i] = float(line[j])
        v[j, i] = float(line[j+3])

# calculate orbit constants
x = x.T
v = v.T
for i in range(n):
    Lhat = np.cross(x[i], v[i])
    # angular momentum
    L = np.linalg.norm(Lhat)
    # normal vector of the orbital plane
    Lhat /= L
    # total energy
    E = 0.5 * np.linalg.norm(v[i])**2 - G / np.linalg.norm(x[i])
    # major axis
    a = - 0.5 * G / E
    # eccentricity
    ecc = (1 + 2 * E * L**2 / G**2)**0.5
    # initial distance
    r0 = np.linalg.norm(x[i])
    # initial direction
    rhat0 = x[i] / r0
    # initial angle
    theta0 = np.arccos((1 - a * (1 - ecc**2) / r0) / ecc)
    if np.dot(x[i], v[i]) > 0:
        theta0 = - theta0

    # calculate the orbit
    theta = np.linspace(0, 2 * math.pi, num=npoint)
    r = a * (1 - ecc**2) / (1 - ecc * np.cos(theta))
    # output in cartesian coordinates
    for j in range(npoint):
        rotation_vector = Lhat * (theta[j] - theta0)
        rotation = R.from_rotvec(rotation_vector)
        rotated_vector = rotation.apply(rhat0 * r[j])
        print('{:d}\t{:e}\t{:e}\t{:e}'.format(i, rotated_vector[0], rotated_vector[1], rotated_vector[2]))
    
        

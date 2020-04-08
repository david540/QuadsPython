#!/usr/bin/env python3

from math import sqrt, acos, cos, sin, asin
import numpy as np
import scipy.sparse.linalg
import cmath
from process import *

def getVectorFieldsBad(dz):
    u, v = [], []
    if len(dz) == 0:
        return u, v
    mean_value = sum(dz) / len(dz)
    mean_angle = cmath.phase(mean_value) / 4
    for z in dz:
        angle = cmath.phase(z) / 4
        while abs(mean_angle - angle) > np.pi/4:
            if angle < mean_angle:
                angle += np.pi/2
            else:
                angle -= np.pi/2
        u.append(cmath.rect(1, angle))
        v.append(cmath.rect(1, (angle + np.pi/2)))
    return u, v

def rotateFieldAboveXvalue(mean_x_value, a, b, faces, vertices):
    for t in range(len(faces)):
        centerX = sum([vertices[n - 1][0] for n in faces[t]]) / 3
        if centerX > mean_x_value:
            a[t] *= cmath.rect(1, np.pi / 2)
            b[t] *= cmath.rect(1, np.pi / 2)




def getVectorFields(dz, neighbours):
    u, v = [0 for i in range(len(dz))], [0 for i in range(len(dz))]
    if len(dz) == 0:
        return u, v
    mean_value = sum(dz) / len(dz)
    actual_angle = cmath.phase(mean_value) / 4
    alreadyDone = [False for i in range(len(dz))]
    idAct = 0
    compteur = 0
    idToTest = [[0, 0]]
    maxi = 0
    while compteur < len(neighbours):
        idAct, actual_angle = idToTest[0]
        if alreadyDone[idAct] == True:
            del idToTest[0]
            continue
        alreadyDone[idAct] = True
        angle = cmath.phase(dz[idAct]) / 4
        if abs(angle) > abs(maxi):
            maxi = angle
        while abs(actual_angle - angle) > np.pi/4:
            if angle < actual_angle:
                angle += np.pi/2
            else:
                angle -= np.pi/2
        u[idAct] = cmath.rect(1, angle)
        v[idAct] = cmath.rect(1, (angle + np.pi/2))
        for n in neighbours[idAct]:
            idToTest.append([n, angle])
        #actual_angle = angle
        del idToTest[0]
        compteur+= 1

    return u, v

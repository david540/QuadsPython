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




def getVectorFields(dz, neighbours, faces):
    u, v = [0 for i in range(len(dz))], [0 for i in range(len(dz))]
    if len(dz) == 0:
        return u, v
    mean_value = sum(dz) / len(dz)
    actual_angle = cmath.phase(mean_value) / 4
    alreadyDone = [False for i in range(len(dz))]
    arbreCouvrant = []
    idAct = 0
    compteur = 0
    idToTest = [[123, cmath.phase(dz[123]) / 4]]
    maxi = 0
    while compteur <= len(neighbours):
        idAct, actual_angle = idToTest[0]
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
            if alreadyDone[n]:
                continue

            alreadyDone[n] = True
            idToTest.append([n, angle])
            idCommonEdge = []
            for i in faces[n]:
                if i in faces[idAct]:
                    idCommonEdge.append(i - 1)
            arbreCouvrant.append(tuple(sorted(idCommonEdge)))
        #actual_angle = angle
        del idToTest[0]
        compteur+= 1

    return u, v, arbreCouvrant

def zipArbreCouvrant(arbreCouvrant, neighbours, singularities):
    pileVertice = []
    count = [0 for i in range(len(neighbours))]
    for i in range(len(neighbours)):
        for n in neighbours[i]:
            key = min(i, n), max(i, n)
            if key not in arbreCouvrant:
                count[i] += 1
        if count[i] == 1:
            pileVertice.append(i)
    while len(pileVertice) > 0:
        curVert = pileVertice[0]
        if curVert in singularities or count[curVert] == 0:
            del(pileVertice[0])
            continue
        for n in neighbours[curVert]:
            key = min(curVert, n), max(curVert, n)
            if key not in arbreCouvrant:
                arbreCouvrant.append(key)
                count[curVert] -= 1
                if count[curVert] != 0:
                    print("ERROR ZIP 940", count[curVert])
                del(pileVertice[0])
                count[n] -= 1
                if count[n] == 1:
                    pileVertice.append(n)

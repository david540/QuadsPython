#!/usr/bin/env python3

from math import sqrt, acos, cos, sin, asin
import numpy as np
import scipy.sparse.linalg
import cmath
from process import *

def defineBoundary(dz):
    idsInside = []
    idsBoundary = []
    for i in range(len(dz)):
        #print(normals[i])
        if dz[i] == 0:
            idsInside.append(i)
        else:
            idsBoundary.append(i)
    return idsInside, idsBoundary

def getMiddleOfVertice(p1, p2):
    return (p2[0] + p1[0])/2, (p2[1] + p1[1])/2

def dist(p, v):
    return (v[0] - p[0]) ** 2 + (v[1] - p[1]) ** 2

def fieldOnAnyPoint(p, vertices, u, gap = 0.001, alpha=0.5):
    uz = 0
    z = complex(p[0], p[1])

    for i in range(len(vertices)):
        distance = pow(dist(vertices[i], p), alpha)
        if distance < gap:
            uz = u[i]
            break
        uz += u[i] / distance
    return uz / abs(uz)

def verifJustesse(A, b, X, neighbours, idsBoundary):

    A_X = A * X
    A_X = A_X.tolist()
    b = b.tolist()
    som = 0
    for i in range(len(neighbours)):
        som += (A_X[i][0] - b[i][0])
    print(som)

def normalize(u, recadre = True):
    length = sqrt(u[0] ** 2 + u[1] ** 2)
    if length == 0:
        return [0, 0]


    vect = [u[0]/length, u[1]/length]

    if recadre:
        while vect[0] <= 0 or vect[1] < 0:
            vect[0], vect[1] = -vect[1], vect[0]
    return vect

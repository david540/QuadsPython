#!/usr/bin/env python3

from math import sqrt, acos, cos, sin, asin
import numpy as np
import scipy.sparse.linalg
import cmath
from process import *


def solveLinearEquationTriangleCenter(dz, idsBoundary, facesIdentification, neifaces):
    dz = [dz[i] for i in range(len(dz))]
    frames = []
    for i in range(len(dz)):
        if dz[i] == 0:
            frames.append(complex(1, 0))
        else:
            frames.append(dz[i])
    #dz = [complex(1,0) for j in range(len(neighbours))]

    A = []
    b = []

    for i in range(len(neifaces)):
        A.append([0 for _ in range(len(neifaces))])
        b.append([0])
        if i in idsBoundary:
            A[i][i] = 1
            b[i][0] = frames[i]
        else:
            if(len(neifaces[i]) == 0):
                print("wow")
            for j in neifaces[i]:
                if j in idsBoundary:
                    b[i][0] += frames[j]
                else:
                    A[i][j] = -1
            A[i][i] = len(neifaces[i])

    #print(b)
    #print(A[-4])
    A = np.matrix(A)
    b = np.matrix(b)

    #print(A)

    X = np.linalg.inv(A) * b

    dz = X.transpose().tolist()[0]
        #print(X)
    #X = ((np.linalg.inv(A)*b).transpose().tolist()[0])
    #verifJustesse(A, b, X, neighbours, idsBoundary)
    for i in range(len(neifaces)):
        #normals[i] = normalize([cos(cmath.phase(X[i]) / 4), sin(cmath.phase(X[i])/ 4)])
        #normals[i] = normalize([cos(cmath.phase(X[i]) / 4), sin(cmath.phase(X[i])/ 4)])
        if dz[i] == 0:
            dz[i] = X[i, 0]

    return dz

def solveLinearEquationComplex(dz, neighbours, idsInside, idsBoundary):
    dz = [dz[i] for i in range(len(dz))]
    frames = []
    for i in range(len(dz)):
        #print(normals[i])
        if dz[i] == 0:
            frames.append(complex(1, 0))
        else:
            frames.append(dz[i])
            #print(normals[i], cmath.rect(1, normals[i][0]))

    #dz = [complex(1,0) for j in range(len(neighbours))]

    A = [[0 for j in range(len(neighbours))] for i in range(len(neighbours))]
    b = [[0] for j in range(len(neighbours))]

    for i in range(len(neighbours)):
        if i in idsBoundary:
            A[i][i] = 1
            b[i][0] = frames[i]
        else:
            for j in neighbours[i]:
                if j in idsBoundary:
                    b[i][0] += frames[j]
                else:
                    A[i][j] = -1
            A[i][i] = len(neighbours[i])

    #print(b)
    A = np.matrix(A)
    b = np.matrix(b)

    X = np.linalg.inv(A) * b


        #print(X)
    #X = ((np.linalg.inv(A)*b).transpose().tolist()[0])
    #verifJustesse(A, b, X, neighbours, idsBoundary)
    for i in range(len(neighbours)):
        #normals[i] = normalize([cos(cmath.phase(X[i]) / 4), sin(cmath.phase(X[i])/ 4)])
        #normals[i] = normalize([cos(cmath.phase(X[i]) / 4), sin(cmath.phase(X[i])/ 4)])
        dz[i] = X[i, 0]

    return dz


def solveLinearEquationComplexMoindresCarres(dz, neighbours, idsInside, idsBoundary):
    C = 100
    dz = [dz[i] for i in range(len(dz))]
    frames = []
    for i in range(len(dz)):
        #print(normals[i])
        if dz[i] == 0:
            frames.append(complex(1, 0))
        else:
            frames.append(dz[i])
            #print(normals[i], cmath.rect(1, normals[i][0]))

    #dz = [complex(1,0) for j in range(len(neighbours))]

    A = []
    b = []

    for i in range(len(neighbours)):
        A.append([0 for _ in range(len(neighbours))])
        b.append([0])
        for j in neighbours[i]:
            if j in idsBoundary:
                b[(len(b) - 1)][0] += frames[j]
            else:
                A[(len(A) - 1)][j] = -1
        A[(len(A) - 1)][i] = len(neighbours[i])

    for i in range(len(neighbours)):
        if i in idsBoundary:
            A.append([0 for _ in range(len(neighbours))])
            b.append([0])
            A[(len(A) - 1)][i] = C
            b[(len(b) - 1)][0] = C * frames[i]

    print(len(A), len(A[0]))

    X = np.linalg.lstsq(A, b)[0]

    for i in range(len(neighbours)):
        #normals[i] = normalize([cos(cmath.phase(X[i]) / 4), sin(cmath.phase(X[i])/ 4)])
        #normals[i] = normalize([cos(cmath.phase(X[i]) / 4), sin(cmath.phase(X[i])/ 4)])
        dz[i] = X[i, 0]

    return dz

def solveLinearEquationComplexAtA(dz, neighbours, idsInside, idsBoundary):
    C = 100
    dz = [dz[i] for i in range(len(dz))]
    frames = []
    for i in range(len(dz)):
        #print(normals[i])
        if dz[i] == 0:
            frames.append(complex(1, 0))
        else:
            frames.append(dz[i])
            #print(normals[i], cmath.rect(1, normals[i][0]))

    #dz = [complex(1,0) for j in range(len(neighbours))]

    A = []
    b = []

    for i in range(len(neighbours)):
        A.append([0 for _ in range(len(neighbours))])
        b.append([0])
        for j in neighbours[i]:
            if j in idsBoundary:
                b[(len(b) - 1)][0] += frames[j]
            else:
                A[(len(A) - 1)][j] = -1
        A[(len(A) - 1)][i] = len(neighbours[i])

    for i in range(len(neighbours)):
        if i in idsBoundary:
            A.append([0 for _ in range(len(neighbours))])
            b.append([0])
            A[(len(A) - 1)][i] = C
            b[(len(b) - 1)][0] = C * frames[i]

    print(len(A), len(A[0]))
    A = np.matrix(A)
    b = np.matrix(b)

    AtA = A.transpose() * A
    Atb = A.transpose() * b

    AtA = np.matrix(AtA)
    Atb = np.matrix(Atb)

    X = np.linalg.inv(AtA) * Atb


    for i in range(len(neighbours)):
        #normals[i] = normalize([cos(cmath.phase(X[i]) / 4), sin(cmath.phase(X[i])/ 4)])
        #normals[i] = normalize([cos(cmath.phase(X[i]) / 4), sin(cmath.phase(X[i])/ 4)])
        dz[i] = X[i, 0]

    return dz

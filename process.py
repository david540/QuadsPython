#!/usr/bin/env python3

from itertools import permutations
from math import sqrt, acos, cos, sin, asin
import numpy as np
import scipy.sparse.linalg
import cmath

PI = 3.14159

def computeNormalsOfVertices(vertices, neighbours, contours, d=600):
    def computeCosAngle(u, v):
        return (u[0] * v[0] + u[1] * v[1]) / ( sqrt(u[0] * u[0] + u[1] * u[1]) * sqrt(v[0] * v[0] + v[1] * v[1]) )

    normals = [[0, 0] for i in range(len(vertices))]

    for i in range(len(vertices)):
        for j in neighbours[i]:
            key = min(i, j), max(i, j)
            if key in contours:
                vmid, vx, vy = computeNormalOfContour(vertices[i], vertices[j], vertices[contours[key][0]], d=d)
                if (normals[i][0] != 0 or  normals[i][1] != 0):
                    while computeCosAngle((normals[i][0], normals[i][1]), (vx, vy)) < 0.7:
                        vx, vy = -vy, vx #A COMPRENDRE POURQUOI uy et pas -uy
                normals[i][0] += vx
                normals[i][1] += vy
        normals[i] = normalize(normals[i])

    return normals

def computeNormalOfContour(vertice1, vertice2, vertice3, d=600):
    vmid = [(vertice2[i] + vertice1[i]) / 2 for i in range(2)]
    ux, uy = [(vertice2[i] - vertice1[i]) for i in range(2)]
    vx, vy = sqrt(uy**2 / (ux**2 + uy**2) * d), sqrt(ux**2 / (ux**2 + uy**2) * d)
    if(abs(vx * ux + vy * uy) > 0.001):
        vx = -vx
    if(abs(vx * ux + vy * uy) > 0.001):
        print("erreur de conception mathématique 1")
    wx, wy = [(vertice3[i] - vmid[i]) for i in range(2)]
    if ux == 0:
        if vy != 0:
            print("erreur de conception mathématique 2 :", vy)
        if vx * wx >= 0:
            vx = -vx
    elif vy * (wy - uy / ux * wx) >= 0 :
        vx, vy = -vx, -vy

    [vx, vy] = normalize([vx, vy])
    return vmid, vx, vy

def propagateConstraints(normals, neighbours, nbiters = 1):
    idsInside = []
    for i in range(len(normals)):
        if normals[i][0] == 0 and normals[i][1] == 0:
            idsInside.append(i)

    newnormals = list(normals)

    for _ in range(nbiters):
        for i in idsInside:
            somme = 0
            count = 0
            for nei in neighbours[i]:
                angle = acos(normals[nei][0])
                if count > 0 :
                    angle_mean = somme/ count
                    if angle - angle_mean > PI/4:
                        angle -= PI/2
                    elif angle - angle_mean < -PI/4:
                        angle += PI/2
                somme += angle
                count += 1

            somme /= count
            somme %= PI/2
            newnormals[i] = cos(somme), sqrt(1 - cos(somme)**2)

        for i in idsInside:
            normals[i] = normalize(newnormals[i])

def solveLinearEquation(normals, neighbours):
    C = 100
    idsInside = []
    idsBoundary = []
    frames = []
    for i in range(len(normals)):
        #print(normals[i])
        if normals[i][0] == 0 and normals[i][1] == 0:
            idsInside.append(i)
            frames.append([1, 0])
        else:
            idsBoundary.append(i)
            frames.append([cos(4 * acos(normals[i][0])), sin(4 * acos(normals[i][0]))])


    A = []
    b = []
    numLine = 0

    for i in range(len(neighbours)):
        if i in idsBoundary:
            A.append([0 for i in range(2 * len(neighbours))])
            A[numLine][2 * i] = C
            b.append([C * frames[i][0]])
            numLine += 1

            A.append([0 for i in range(2 * len(neighbours))])
            A[numLine][2 * i + 1] = C
            b.append([C * frames[i][1]])
            numLine += 1


        for j in neighbours[i]:
            A.append([0 for i in range(2 * len(neighbours))])
            A[numLine][2 * i] = sqrt(PI)
            A[numLine][2 * j] = - sqrt(PI)
            b.append([0])
            numLine += 1

            A.append([0 for i in range(2 * len(neighbours))])
            A[numLine][2 * i + 1] = sqrt(PI)
            A[numLine][2 * j + 1] = - sqrt(PI)
            b.append([0])
            numLine += 1

    #print(b)
    A = np.matrix(A)
    b = np.matrix(b)

    #print(b)


    At = A.getT()


    At_A = At * A
    At_b = At * b

    #print(At_b)

    X = np.linalg.inv(At * A) * At * b

    #X = ((np.linalg.inv(At_A)*At_b).transpose().tolist()[0])

    #X = (np.linalg.inv(At_A)*At_b)

    #X = matrixProduct3(np.linalg.inv(A), b)
    #for i in range(len(X)):
    #    X[i] %= PI/2
    #print(X)
    #test()

    #print(At_A, At_b)
    #print(At_A.shape())
    #print(At_B.shape())

    #X = np.matrix(scipy.sparse.linalg.cg(At_A, At_b)[0]).transpose()
    #print(X)

    verifJustesse(A, b, X, neighbours)

    A_X = A * X
    #X = X.transpose().tolist()[0]
    #X = b
    #for i in idsBoundary:
#        X[i] = b[i]
    #print(X,"\n\n\n")
    #for i in range(len(neighbours)):
        #print(X[2 * i], X[2 * i + 1])
        #normals[i] = normalize([X[2 * i], X[2 * i + 1]])
    #print(normals)


def solveLinearEquationComplex(normals, neighbours):
    idsInside = []
    idsBoundary = []
    frames = []
    for i in range(len(normals)):
        #print(normals[i])
        if normals[i][0] == 0 and normals[i][1] == 0:
            idsInside.append(i)
            frames.append(complex(1, 0))
        else:
            idsBoundary.append(i)
            frames.append(cmath.rect(1, 4 * acos(normals[i][0])))
            #print(normals[i], cmath.rect(1, normals[i][0]))

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



    #X = ((np.linalg.inv(A)*b).transpose().tolist()[0])
    #verifJustesse(A, b, X, neighbours, idsBoundary)
    for i in range(len(neighbours)):
        #normals[i] = normalize([cos(cmath.phase(X[i]) / 4), sin(cmath.phase(X[i])/ 4)])
        normals[i] = normalize([cos(cmath.phase(X[i]) / 4), sin(cmath.phase(X[i])/ 4)])


def verifJustesse(A, b, X, neighbours, idsBoundary):

    A_X = A * X
    A_X = A_X.tolist()
    b = b.tolist()
    som = 0
    for i in range(len(neighbours)):
        #print(A_X[i][0])
        som += (A_X[i][0] - b[i][0])
        #print(A_X[i] - b[i])

    for i in idsBoundary:
        print(X[i], b[i])


    print(som)
    #for i in range(len(neighbours)):
    #    som = A[i, i] * X[i]
    #    for j in neighbours[i]:
    #        som += A[i, j] * X[j]
    #    print(som - b[i, 0])

def findNeighboursAndContours(nbVertices, faces):
    neighbours = [[] for i in range(nbVertices)]
    normalsVect = [[0, 0] for i in range(nbVertices)]
    contours = {}
    for face in faces:
        for v1, v2, v3 in permutations([face[0][0] - 1, face[1][0] - 1, face[2][0] - 1], 3):
            key = (min(v1, v2), max(v1, v2))
            if v2 not in neighbours[v1]:
                neighbours[v1].append(v2)
                contours[key] = [v3]
            else:
                if key in contours:
                    del(contours[key])

    return neighbours, contours

def normalize(u):
    length = sqrt(u[0] ** 2 + u[1] ** 2)
    if length == 0:
        return [0, 0]


    vect = [u[0]/length, u[1]/length]

    while vect[0] <= 0 or vect[1] < 0:
        vect[0], vect[1] = -vect[1], vect[0]
    return vect

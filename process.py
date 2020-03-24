#!/usr/bin/env python3

from itertools import permutations
from math import sqrt, acos, cos, sin, asin
import numpy as np
import scipy.sparse.linalg
import cmath


def computeNormalsOfVertices(vertices, neighbours, contours, d=600):
    def computeCosAngle(u, v):
        return (u[0] * v[0] + u[1] * v[1]) / ( sqrt(u[0] * u[0] + u[1] * u[1]) * sqrt(v[0] * v[0] + v[1] * v[1]) )

    dz = [complex(0, 0) for i in range(len(vertices))]
    quadrants = [0 for i in range(len(vertices))]

    corners = []

    for i in range(len(vertices)):
        vmid, vx, vy = 0, 0, 0
        initial_z = complex(1, 0)
        for j in neighbours[i]:
            key = min(i, j), max(i, j)
            if key in contours:
                vmid, vx, vy = computeNormalOfContour(vertices[i], vertices[j], vertices[contours[key][0]], d=d)
                initial_z = complex(vx, vy)
                z = cmath.rect(1, 4 * cmath.phase(initial_z))
                dz[i] += z

        if dz[i] != 0:
            dz[i] /= abs(dz[i])
            while abs((cmath.phase(dz[i]) / 4 - cmath.phase(initial_z) + np.pi/2 * quadrants[i])) > np.pi /4:
                if cmath.phase(dz[i]) / 4 + np.pi/2 * quadrants[i] < cmath.phase(initial_z):
                    quadrants[i] += 1
                else:
                    quadrants[i] -= 1

    return dz, quadrants

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

    [vx, vy] = normalize([vx, vy], False)
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
                    if angle - angle_mean > np.pi/4:
                        angle -= np.pi/2
                    elif angle - angle_mean < -np.pi/4:
                        angle += np.pi/2
                somme += angle
                count += 1

            somme /= count
            somme %= np.pi/2
            newnormals[i] = cos(somme), sqrt(1 - cos(somme)**2)

        for i in idsInside:
            normals[i] = normalize(newnormals[i])

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



def computeIsoValueLines(neighbours, dz, vertices, idsBoundary, quadrants, eps = 0.1):
    def dist(p, v):
        return sqrt((v[0] - p[0]) ** 2 + (v[1] - p[1]) ** 2)

    startingCorner = 0
    for i in range(len(neighbours)):
        if len(neighbours[i]) == 2:
            startingCorner = i
            break

    lines = []

    cornerCount = 0
    boundId = neighbours[startingCorner][0]
    oldId = startingCorner
    oldoldId = neighbours[startingCorner][1]
    for loop in range(40):
        closestId = boundId
        p = vertices[boundId]
        lines.append([p])
        oldAngle = cmath.phase(dz[boundId]) / 4 + quadrants[boundId] * np.pi/2
        #while dist([cos(oldAngle), sin(oldAngle)], normals[boundId]) > 0.01:
        #    oldAngle += np.pi/2
        #for i in range(0, 4):

        oldAngle += np.pi
        for _ in range(300):
            distance = dist(p, vertices[closestId])
            if distance < 0.01:
                curDz = dz[closestId]
            else:
                curDz = dz[closestId] / distance
                for n in neighbours[closestId]:
                    distance = dist(p, vertices[n])
                    if distance < 0.01:
                        curDz = dz[n]
                        break
                    curDz += dz[n] / distance

            angle = cmath.phase(curDz) / 4
            while abs(angle - oldAngle) > np.pi/4:
                if angle < oldAngle:
                    angle += np.pi/2
                else:
                    angle -= np.pi/2
            angle %= 2*np.pi
            if _ == 0:
                #print(angle, oldAngle)
                angle = oldAngle
            p= [p[0] + eps * cos(angle), p[1] + eps * sin(angle)]

            #for i in range(neighbours[])

            lines[-1].append(p)

            argmin = closestId
            distmin = dist(p, vertices[closestId])
            for i in neighbours[closestId]:
                distance = dist(p, vertices[i])
                if distance < distmin :
                    distmin = distance
                    argmin = i
                    #print(argmin)

            closestId = argmin
            oldAngle = angle

        for n in neighbours[boundId]:
            if n in idsBoundary and n != oldId and n != oldoldId:
                boundId, oldId, oldoldId = n, boundId, oldId




    return lines


def verifJustesse(A, b, X, neighbours, idsBoundary):

    A_X = A * X
    A_X = A_X.tolist()
    b = b.tolist()
    som = 0
    for i in range(len(neighbours)):
        #print(A_X[i][0])
        som += (A_X[i][0] - b[i][0])
        #print(A_X[i] - b[i])

    #for i in idsBoundary:
    #    print(X[i], b[i])


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

def normalize(u, recadre = True):
    length = sqrt(u[0] ** 2 + u[1] ** 2)
    if length == 0:
        return [0, 0]


    vect = [u[0]/length, u[1]/length]

    if recadre:
        while vect[0] <= 0 or vect[1] < 0:
            vect[0], vect[1] = -vect[1], vect[0]
    return vect

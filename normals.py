#!/usr/bin/env python3

from math import sqrt, acos, cos, sin, asin
import numpy as np
import scipy.sparse.linalg
import cmath
from process import *
from itertools import combinations, permutations

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

def computeNormalOfTriangles(vertices, faces, contours, d=600):
    def computeCosAngle(u, v):
        return (u[0] * v[0] + u[1] * v[1]) / ( sqrt(u[0] * u[0] + u[1] * u[1]) * sqrt(v[0] * v[0] + v[1] * v[1]) )

    dz = [complex(0, 0) for i in range(len(faces))]
    quadrants = [0 for i in range(len(faces))]

    corners = []

    for n, ids in enumerate(faces):
        vmid, vx, vy = 0, 0, 0
        initial_z = complex(1, 0)
        for i, j, k in permutations(ids, 3):
            #print([i,j], key)
            i, j, k = i-1, j-1, k-1
            key = i, j
            if key in contours:
                initial_z = computeNormalOfContourComplex(vertices[i], vertices[j], vertices[k])
                #z = cmath.rect(1, 4 * cmath.phase(vz))
                dz[n] = cmath.rect(1, 4 * cmath.phase(initial_z))

                #dz[n] /= abs(dz[n])
                while abs((cmath.phase(dz[n]) / 4 - cmath.phase(initial_z) + np.pi/2 * quadrants[n])) > np.pi /4:
                    if cmath.phase(dz[n]) / 4 + np.pi/2 * quadrants[n] < cmath.phase(initial_z):
                        quadrants[n] += 1
                    else:
                        quadrants[n] -= 1
                break

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

def computeNormalOfContourComplex(vertice1, vertice2, vertice3):
    zu = complex(vertice2[0] - vertice1[0], vertice2[1] - vertice1[1])
    zv = cmath.rect(1, cmath.phase(zu) + np.pi/2)
    zw = complex(vertice3[0] - vertice1[0], vertice3[1] - vertice1[1])
    ztest = zv/zw
    if ztest.real > 0:
        #angle entre v et w est inférieure à pi/2, donc la normale est orienté dans le mauvais sens

        #rotation de pi
        zv *= cmath.rect(1, np.pi)

    return zv

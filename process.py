#!/usr/bin/env python3

from itertools import permutations
from math import sqrt, acos, cos


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

def propagateConstraints(normals, neighbours, nbiters = 1000):
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

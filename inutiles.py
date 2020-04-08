#!/usr/bin/env python3

from math import sqrt, acos, cos, sin, asin
import numpy as np
import scipy.sparse.linalg
import cmath
from process import *

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



def computeIsoValueLinesRemake(neighbours, dz, vertices, idsBoundary, quadrants, eps = 0.1):
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
            curDz = fieldOnAnyPoint(p, vertices, dz)
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

            oldAngle = angle

        for n in neighbours[boundId]:
            if n in idsBoundary and n != oldId and n != oldoldId:
                boundId, oldId, oldoldId = n, boundId, oldId




    return lines

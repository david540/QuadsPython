#!/usr/bin/env python3

from math import sqrt, acos, cos, sin, asin
import numpy as np
import scipy.sparse.linalg
import cmath
from process import *
from itertools import permutations, combinations


def calculCarteIndirection(vertices, faces):
    adjFaces = [[] for i in range(len(vertices))]
    for i, face in enumerate(faces):
        for idv in face:
            adjFaces[idv - 1].append(i)
    return adjFaces

def computeSingularities(faces, dz, vertices, neifaces, adjFaces):

    singularities = []
    for idv in range(len(adjFaces)):
        cur_idf = adjFaces[idv][0]
        first_idf = cur_idf
        curNeighVertice = 1
        if faces[cur_idf][curNeighVertice] - 1 == idv:
            curNeighVertice = 0
        #print(curNeighVertice)
        sumRotations = 0
        previdf = cur_idf
        #print("\n\n", idv," \n\n")

        for _ in range(len(adjFaces[idv])):
            found = False

            for i, j, k in permutations(range(3), 3):
                if faces[cur_idf][i] - 1 == idv and k != curNeighVertice:
                    nextEdge = [faces[cur_idf][k] - 1, idv]
                    for idf in adjFaces[idv]:
                        if idf != cur_idf and idf != previdf:
                            #print(idf, cur_idf, previdf)
                            if sum([faces[idf][ite] in faces[cur_idf] for ite in range(3)]) == 2:
                                for ite in range(3):
                                    if faces[idf][ite] != faces[cur_idf][curNeighVertice] and faces[idf][ite] - 1 != vertices[idv]:
                                        #print("oui", cmath.phase(dz[idf]), cmath.phase(dz[cur_idf]))
                                        rotation = (cmath.phase(dz[idf]) / 4 - cmath.phase(dz[cur_idf]) /4)
                                        #if abs(rotation) > np.pi/8:
                                        #    print(rotation, idv)
                                        while rotation < -np.pi/4:
                                            rotation += np.pi/2
                                        while rotation > np.pi/4:
                                            rotation -= np.pi/2

                                        sumRotations += rotation
                                        previdf = cur_idf
                                        cur_idf = idf
                                        curNeighVertice = ite
                                        found = True
                                        break
                        if found:
                            break
                if found:
                    break

            #print(sumRotations)
            #print(cur_idf)


            if first_idf == cur_idf:
                break

        if abs(sumRotations) > np.pi/4:
            #print(sumRotations)
            singularities.append(idv)

    #print(singularities)
    return singularities

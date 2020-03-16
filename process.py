#!/usr/bin/env python3

from itertools import permutations

def findNeighboursAndContours(nbVertices, faces):
    neighbours = [[] for i in range(nbVertices)]
    contours = [[] for i in range(nbVertices)]
    for face in faces:
        for v1, v2 in permutations([face[0][0] - 1, face[1][0] - 1, face[2][0] - 1], 2):
            if v2 not in neighbours[v1]:
                neighbours[v1].append(v2)
                contours[v1].append(v2)
            else:
                if v2 in contours[v1]:
                    contours[v1].remove(v2)

    return neighbours, contours

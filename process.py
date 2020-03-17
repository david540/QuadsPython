#!/usr/bin/env python3

from itertools import permutations

def findNeighboursAndContours(nbVertices, faces):
    neighbours = [[] for i in range(nbVertices)]
    contours = {}
    for face in faces:
        for v1, v2, v3 in permutations([face[0][0] - 1, face[1][0] - 1, face[2][0] - 1], 3):
            key = (min(v1, v2), max(v1, v2))
            if v2 not in neighbours[v1]:
                neighbours[v1].append(v2)
                contours[key] = v3
            else:
                if key in contours:
                    del(contours[key])

    return neighbours, contours

#!/usr/bin/env python3
from sys import argv

from drawer import drawObj

import numpy as np
import cmath
from random import randint, random
from parser import *
from math import cos, sin, sqrt, acos
from process import *
from scalarfield import *
from normals import *
from singularities import *
from vectorfield import *
from frame import *
from disjointset import *
from inutiles import computeIsoValueLinesRemake

def main(filename):
    print("Initialisation")

    vertices, textureCoord, normals, faces = parseObj(filename)

    neighbours, contours, facesIdentification, neifaces = findNeighboursAndContours(vertices, faces)

    print("Calcul des normales")

    dz, quadrants = computeNormalOfTriangles(vertices, faces, contours, d=100)

    idsInside, idsBoundary = defineBoundary(dz)

    print("Calcul des croix sur les sommets")

    dz= solveLinearEquationTriangleCenter(dz, idsBoundary, facesIdentification, neifaces)

    print("Calcul des champs vectoriels")

    a, b = 0, 0
    a, b = getVectorFields(dz, neifaces)
    #print(a, b)
    #mean_x_value = (max(vertices)[0] + min(vertices)[0]) / 2
    #for _ in range(2):
    #    rotateFieldAboveXvalue(mean_x_value, a, b, faces, vertices)
    #print(mean_x_value)
    #print(a)

    for i in range(len(a)):
        k = randint(1, 3)
        #print(a[i])
        #a[i] *= cmath.rect(1, k * np.pi/2)
        #b[i] *= cmath.rect(1, k * np.pi/2)

    print("Calcul des singularites")
    singularities = 0
    #singularities = computeSingularities(faces, dz)


    print("Calcul du disjoint set")
    ds, singularities = setDisjointSet(vertices, faces, neifaces, a)

    print("Calcul des champs scalaires")
    u, v = 0, 0
    u, v = scalarFieldLinearSystemDisjointSet(faces, a, b, neifaces, ds, vertices, singularities)

    #print(u, v)

    #u = scalarFieldMiddleTriangle(faces, a, vertices)

    #v = scalarFieldLinearSystem(vertices, b, neighbours)

    print("Affichage")
    drawObj(faces, vertices, neighbours, contours, dz, quadrants,\
        idsBoundary, singularities, v, u, a, b)


if __name__ == "__main__":
    filename = argv[1] if len(argv) > 1 else "obj/duck.obj"

    print(filename)
    main(filename)

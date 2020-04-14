#!/usr/bin/env python3

from tkinter import *

from math import cos, sin, sqrt, acos

import numpy as np
from normals import computeNormalOfContourComplex
from process import getMiddleOfVertice
import cmath
import colorsys
from itertools import permutations, combinations


centreX, centreY = 0, 0
canvas_width = 2000
canvas_height = 2000
scale = 1

def initGlobalVariables(vertices):
    global centreX, centreY, canvas_width, canvas_height, scale
    listXandY = [[v[i] for v in vertices] for i in range(2)]
    minX, minY, maxX, maxY = min(listXandY[0]), min(listXandY[1]), max(listXandY[0]), max(listXandY[1])
    canvas_width = 1000
    canvas_height = 1000

    dimension = 0.8
    deltaX = (maxX - minX)
    deltaY = (maxY - minY)

    debutX = (max(deltaX / canvas_width, deltaY / canvas_height) * canvas_width - dimension * deltaX) / (2 * dimension * deltaX)
    debutY = (max(deltaX / canvas_width, deltaY / canvas_height) * canvas_width - dimension * deltaY) / (2 * dimension * deltaY)

    centreX = minX - (maxX - minX) * debutX
    centreY = minY- (maxY - minY) * debutY
    scale = max((maxX - minX) / (dimension * canvas_width), (maxY - minY) / (dimension * canvas_height))

def drawCross(canvas, vertice, lengthRay, angle, color, width):
    canvas.create_line(setCoordX(vertice[0] - lengthRay * cos(angle)), setCoordY(vertice[1] + lengthRay * sin(angle)),
                    setCoordX(vertice[0] + lengthRay * cos(angle)), setCoordY(vertice[1] - lengthRay * sin(angle)),
                    fill=color, width=width)

    canvas.create_line(setCoordX(vertice[0] + lengthRay * sin(angle)), setCoordY(vertice[1] + lengthRay * cos(angle)),
                    setCoordX(vertice[0] - lengthRay * sin(angle)), setCoordY(vertice[1] - lengthRay * cos(angle)),
                    fill=color, width=width)

def drawLine(canvas, vertice1, vertice2, color, width):
    canvas.create_line(setCoordX(vertice1[0]), setCoordY(vertice1[1]),
                    setCoordX(vertice2[0]), setCoordY(vertice2[1]),
                    fill=color, width=width)

def drawMultipleLines(canvas, line, color="purple", width=3):
    for i in range(len(line) - 1):
        drawLine(canvas, line[i], line[i + 1], color, width)

def drawNormalVector(canvas, origin, dz, quadrant, color="red", width=1, d=100):
    if (dz == 0):
        return
    #u = vecteur P2 - P1, v = vecteur normal à u et orienté du côté opposé au point P3 (pour sortir du contour)
    #print("cc")
    v = cmath.rect(1, cmath.phase(dz)/ 4 + quadrant * np.pi/2)
    canvas.create_line(setCoordX(origin[0]), setCoordY(origin[1]),
                        setCoordX(origin[0]) + d * v.real, setCoordY(origin[1]) - d * v.imag,
                        fill=color, width=width)


def drawField(canvas, origin, v, color="red", width=2, d=20):
    if (v == 0):
        return
    #u = vecteur P2 - P1, v = vecteur normal à u et orienté du côté opposé au point P3 (pour sortir du contour)
    #print("cc")
    #v = cmath.rect(1, cmath.phase(dz))
    canvas.create_line(setCoordX(origin[0]), setCoordY(origin[1]),
                        setCoordX(origin[0]) + d * v.real, setCoordY(origin[1]) - d * v.imag,
                        fill=color, width=width)

def drawSquare(canvas, origin, c, color="yellow"):
    canvas.create_rectangle(setCoordX(origin[0]) - c, setCoordY(origin[1]) - c, setCoordX(origin[0]) + c, setCoordY(origin[1]) + c, fill=color)

def drawFrameField(canvas, origin, dz, color="orange", width=3, d=20):

    #PI = 3.14159

    dx = cos(cmath.phase(dz) / 4)
    dy = sin(cmath.phase(dz) / 4)

    #u = vecteur P2 - P1, v = vecteur normal à u et orienté du côté opposé au point P3 (pour sortir du contour)
    #red = int(max(0, 255* 2/PI * angle))
    #green = 0 #int(max(0, 3/PI * angle * 255 if angle < 3/PI else 255 * (1 - 3/PI * (angle - PI / 3))))
    #blue = 0 #if int(max(0, 3/PI * (angle - PI/ 3) * 255 if angle < 3/PI else 255 * (1 - 3/PI * (angle - 2*PI / 3))))
    #color = "#" + '{:02x}'.format(red) + '{:02x}'.format(green)  + '{:02x}'.format(blue)

    #print(dz, setCoordX(origin[0]) -dx * d/2)
    canvas.create_line(setCoordX(origin[0]) - dx * d/2, setCoordY(origin[1]) + dy * d/2,
                        setCoordX(origin[0]) + dx * d/2, setCoordY(origin[1]) - dy * d/2,
                        fill=color, width=width)

    canvas.create_line(setCoordX(origin[0]) + dy * d/2, setCoordY(origin[1]) + dx * d/2,
                        setCoordX(origin[0]) - dy * d/2, setCoordY(origin[1]) - dx * d/2,
                        fill=color, width=width)


def setCoordX(curX):
    #return centreX + curX * scale
    return (curX - centreX) / scale

def setCoordY(curY):
    #return centreY - curY * scale
    return canvas_height - (curY - centreY) / scale

def drawScalarField(canvas, u, v, faces, vertices):
    m = max([max(u), max(v)])*(1.001)
    variance = 1
    for t in range(len(faces)):
        center = [sum([vertices[n - 1][j] for n in faces[t]]) / 3 for j in range(2)]

        i, j, k = faces[t][0]-1, faces[t][1]-1, faces[t][2]-1
        centerij = [(vertices[i][coord] + vertices[j][coord]) / 2 for coord in range(2)]
        centerjk = [(vertices[j][coord] + vertices[k][coord]) / 2 for coord in range(2)]
        centerki = [(vertices[k][coord] + vertices[i][coord]) / 2 for coord in range(2)]
        rgb = colorsys.hsv_to_rgb((u[3*t]%(m/variance)) / (m/variance), 1, 1)
        color = "#" + '{:02x}'.format(int(rgb[0] * 255)) + '{:02x}'.format(int(rgb[1] * 255))  + '{:02x}'.format(int(rgb[2] * 255))
        canvas.create_polygon(setCoordX(vertices[i][0]), setCoordY(vertices[i][1]),
                            setCoordX(centerij[0]), setCoordY(centerij[1]),
                            setCoordX(centerki[0]), setCoordY(centerki[1]), fill=color)
        rgb = colorsys.hsv_to_rgb((u[3*t + 1]%(m/variance)) / (m/variance), 1, 1)
        color = "#" + '{:02x}'.format(int(rgb[0] * 255)) + '{:02x}'.format(int(rgb[1] * 255))  + '{:02x}'.format(int(rgb[2] * 255))
        canvas.create_polygon(setCoordX(vertices[j][0]), setCoordY(vertices[j][1]),
                            setCoordX(centerij[0]), setCoordY(centerij[1]),
                            setCoordX(centerjk[0]), setCoordY(centerjk[1]), fill=color)
        rgb = colorsys.hsv_to_rgb((u[3*t + 2]%(m/variance)) / (m/variance), 1, 1)
        color = "#" + '{:02x}'.format(int(rgb[0] * 255)) + '{:02x}'.format(int(rgb[1] * 255))  + '{:02x}'.format(int(rgb[2] * 255))
        canvas.create_polygon(setCoordX(vertices[k][0]), setCoordY(vertices[k][1]),
                            setCoordX(centerjk[0]), setCoordY(centerjk[1]),
                            setCoordX(centerki[0]), setCoordY(centerki[1]), fill=color)
        rgb = colorsys.hsv_to_rgb(sum([(u[3*t + indice]%(m/variance)) / (m/variance) for indice in range(3)]) / 3, 1, 1)
        #rgb = colorsys.hsv_to_rgb((u[3*t]%(m/variance)) / (m/variance), 1, 1)
        color = "#" + '{:02x}'.format(int(rgb[0] * 255)) + '{:02x}'.format(int(rgb[1] * 255))  + '{:02x}'.format(int(rgb[2] * 255))
        canvas.create_polygon(setCoordX(centerjk[0]), setCoordY(centerjk[1]),
                            setCoordX(centerij[0]), setCoordY(centerij[1]),
                            setCoordX(centerki[0]), setCoordY(centerki[1]), fill=color)


            #drawSquare(canvas, vertices[faces[t][k] - 1], 8, color= color)

def drawScalarFieldMiddleTriangle(canvas, u, faces):
    m = max(u)
    for i in range(len(faces)):
        center = [sum([vertices[n - 1][j] for n in faces[i]]) / 3 for j in range(2)]
        rgb = colorsys.hsv_to_rgb(u[i]/m, 1, 1)
        color = "#" + '{:02x}'.format(int(rgb[0] * 255)) + '{:02x}'.format(int(rgb[1] * 255))  + '{:02x}'.format(int(rgb[2] * 255))
        drawSquare(canvas, center, 8, color= color)

def drawNormalVectorSegm(canvas, origin, z, d=100, color="red", width=1):
    if (z == 0):
        return
    #u = vecteur P2 - P1, v = vecteur normal à u et orienté du côté opposé au point P3 (pour sortir du contour)
    #print("cc")
    #v = cmath.rect(1, cmath.phase(dz)/ 4 + quadrant * np.pi/2)
    canvas.create_line(setCoordX(origin[0]), setCoordY(origin[1]),
                        setCoordX(origin[0]) + d * z.real, setCoordY(origin[1]) - d * z.imag,
                        fill=color, width=width)


def drawObj(faces, vertices, neighbours, contours, dz, quadrants,\
    idsBoundary, singularities, u, v, a, b, arbreCouvrant):


    master = Tk()

    master.title("Quads")

    show_triangles = True
    show_normal_vectors_segments_contour = False
    show_normal_vectors_vertices_contour = False
    show_normal_vectors_faces = False
    show_scalar_field = False
    show_contour = False
    show_links_between_vertices = True
    show_random_1ring = False
    show_frame_fields = False
    show_isovaluelines = False
    show_vector_fields = False
    show_arbre_couvrant = True

    initGlobalVariables(vertices)
    xmin, ymin, zmin = [int(t) - 1 for t in min(vertices)]
    xmax, ymax, zmax = [int(t) + 1 for t in max(vertices)]

    canvas = Canvas(master,
           width=canvas_width,
           height=canvas_height)
    canvas.pack()

    canvas.create_rectangle(0, 0, canvas_width, canvas_height, fill= "white")

    if show_triangles:
        for face in faces:
            v1, v2, v3 = face[0] - 1, face[1] - 1, face[2] - 1
            canvas.create_polygon(setCoordX(vertices[v1][0]), setCoordY(vertices[v1][1]),
                                setCoordX(vertices[v2][0]), setCoordY(vertices[v2][1]),
                                setCoordX(vertices[v3][0]), setCoordY(vertices[v3][1]))

    for n, face in enumerate(faces):
        for i, j, k in permutations(face, 3):
            i, j, k = i - 1, j - 1, k - 1
            key = i, j
            if key in contours:
                if show_normal_vectors_segments_contour:
                    zv = computeNormalOfContourComplex(vertices[i], vertices[j], vertices[k])
                    #drawNormalVectorSegm(canvas, getMiddleOfVertice(vertices[i], vertices[j]), zv,)

                    drawNormalVector(canvas, getMiddleOfVertice(vertices[i], vertices[j]), dz[n], quadrants[n])
                if show_contour:
                    drawLine(canvas, vertices[i], vertices[j], "green", 5)
            #else:
            #    if show_links_between_vertices:
            #        drawLine(canvas, vertices[i], vertices[j], "grey", 1)



    if show_random_1ring :
        random_id = randint(0, len(vertices))
        for j in neighbours[random_id]:
            drawLine(canvas, vertices[random_id], vertices[j], "yellow", 2)


    #propagateConstraints(normals, neighbours)

    #if show_constraints:
    #    for i in range(len(vertices)):
    #        drawConstraint(canvas, vertices[i], z[i])

    if show_normal_vectors_vertices_contour:
        for i in range(len(vertices)):
            if i in idsBoundary:
                drawNormalVector(canvas, vertices[i], dz[i], quadrants[i])

    if show_normal_vectors_faces:
        for i in range(len(faces)):
            if i in idsBoundary:
                drawNormalVector(canvas, faces[i], dz[i], quadrants[i])


    if show_isovaluelines:
        lines = computeIsoValueLinesRemake(neighbours, dz, vertices, idsBoundary, quadrants)
        for line in lines:
            drawMultipleLines(canvas, line)



    #if show_triangles and False:


    if show_frame_fields:
        for i in range(len(faces)):
            #print([n for n in faces[i]])
            center = [sum([vertices[n - 1][j] for n in faces[i]]) / 3 for j in range(2)]
            drawFrameField(canvas, center, dz[i], d=10, width=2)

    #print(v)
    if show_scalar_field:
        drawScalarField(canvas, u, v, faces, vertices)
    #drawScalarFieldMiddleTriangle(canvas, u, faces)

    for n, face in enumerate(faces):
        for i, j, k in permutations(face, 3):
            i, j, k = i- 1, j-1, k-1
            if show_links_between_vertices:
                drawLine(canvas, vertices[i], vertices[j], "red", 2)
    if show_vector_fields:
        for i in range(len(faces)):
            center = [sum([vertices[n - 1][j] for n in faces[i]]) / 3 for j in range(2)]
            drawField(canvas, center, a[i], d=10)
            drawField(canvas, center, b[i], d=10, color="orange")

    if show_arbre_couvrant:
        for i, j in arbreCouvrant:
            drawLine(canvas, vertices[i], vertices[j], "green", 2)
    if False:
        for s1, s2 in combinations(singularities, 2):
            v1 = faces[s1[0]][s1[2]] - 1
            v2 = faces[s2[0]][s2[2]] - 1
            print(v2, len(neighbours))
            if v1 in neighbours[v2]:
                drawLine(canvas, vertices[v1], vertices[v2], "green", 5)

    mainloop()

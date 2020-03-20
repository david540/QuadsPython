#!/usr/bin/env python3

from tkinter import *

import numpy as np

from random import randint, random
from parser import parseObj
from math import cos, sin, sqrt, acos
from process import findNeighboursAndContours, computeNormalOfContour,\
computeNormalsOfVertices, propagateConstraints, solveLinearEquationComplex


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

def drawNormalVector(canvas, origin, vx, vy, color="blue", width=1, d=20):
    if (vx == 0 and vy == 0):
        return
    #u = vecteur P2 - P1, v = vecteur normal à u et orienté du côté opposé au point P3 (pour sortir du contour)

    canvas.create_line(setCoordX(origin[0]), setCoordY(origin[1]),
                        setCoordX(origin[0]) + vx, setCoordY(origin[1]) - vy,
                        fill=color, width=width)

def drawConstraint(canvas, origin, vx, vy, color="orange", width=3, d=20):

    if (vx == 0 and vy == 0):
        vx = 1

    angle = acos(vx)

    vx *= d
    vy *= d

    #PI = 3.14159

    #u = vecteur P2 - P1, v = vecteur normal à u et orienté du côté opposé au point P3 (pour sortir du contour)
    #red = int(max(0, 255* 2/PI * angle))
    #green = 0 #int(max(0, 3/PI * angle * 255 if angle < 3/PI else 255 * (1 - 3/PI * (angle - PI / 3))))
    #blue = 0 #if int(max(0, 3/PI * (angle - PI/ 3) * 255 if angle < 3/PI else 255 * (1 - 3/PI * (angle - 2*PI / 3))))
    #color = "#" + '{:02x}'.format(red) + '{:02x}'.format(green)  + '{:02x}'.format(blue)


    canvas.create_line(setCoordX(origin[0]) - vx/2, setCoordY(origin[1]) + vy/2,
                        setCoordX(origin[0]) + vx/2, setCoordY(origin[1]) - vy/2,
                        fill=color, width=width)

    canvas.create_line(setCoordX(origin[0]) + vy/2, setCoordY(origin[1]) + vx/2,
                        setCoordX(origin[0]) - vy/2, setCoordY(origin[1]) - vx/2,
                        fill=color, width=width)


def setCoordX(curX):
    #return centreX + curX * scale
    return (curX - centreX) / scale

def setCoordY(curY):
    #return centreY - curY * scale
    return canvas_height - (curY - centreY) / scale

def drawObj(filename):

    master = Tk()

    master.title("Quads")

    vertices, textureCoord, normals, faces = parseObj(filename)

    initGlobalVariables(vertices)

    neighbours, contours = findNeighboursAndContours(len(vertices), faces)

    normals = computeNormalsOfVertices(vertices, neighbours, contours, d=100)


    canvas = Canvas(master,
           width=canvas_width,
           height=canvas_height)
    canvas.pack()

    canvas.create_rectangle(0, 0, canvas_width, canvas_height, fill= "white")


    show_triangles = True
    show_normal_vectors_segments_contour = False
    show_normal_vectors_vertices_contour = False
    show_contour = True
    show_links_between_vertices = True
    show_random_1ring = True
    show_constraints = False

    if show_triangles:
        for face in faces:
            v1, v2, v3 = face[0][0] - 1, face[1][0] - 1, face[2][0] - 1
            canvas.create_polygon(setCoordX(vertices[v1][0]), setCoordY(vertices[v1][1]),
                                setCoordX(vertices[v2][0]), setCoordY(vertices[v2][1]),
                                setCoordX(vertices[v3][0]), setCoordY(vertices[v3][1]))

    for i in range(len(vertices)):
        for j in neighbours[i]:
            if j < i:
                continue
            key = i, j
            if key in contours:
                if show_normal_vectors_segments_contour:
                    vmid, vx, vy = computeNormalOfContour(vertices[i], vertices[j], vertices[contours[key][0]])
                    drawNormalVector(canvas, vmid, vx, vy)
                if show_contour:
                    drawLine(canvas, vertices[i], vertices[j], "green", 5)
            else:
                if show_links_between_vertices:
                    drawLine(canvas, vertices[i], vertices[j], "grey", 1)

    if show_random_1ring :
        random_id = randint(0, len(vertices))
        for j in neighbours[random_id]:
            drawLine(canvas, vertices[random_id], vertices[j], "yellow", 2)

    if show_normal_vectors_vertices_contour:
        for i in range(len(vertices)):
            drawNormalVector(canvas, vertices[i], normals[i][0], normals[i][1])

    #propagateConstraints(normals, neighbours)

    if show_constraints:
        for i in range(len(vertices)):
            drawConstraint(canvas, vertices[i], normals[i][0], normals[i][1])

    solveLinearEquationComplex(normals, neighbours)

    for i in range(len(vertices)):
        drawConstraint(canvas, vertices[i], normals[i][0], normals[i][1])

    #for vertice in vertices:
    #    drawCross(canvas, vertice, 10, random() * 3.14/2, "orange", 4)



    mainloop()

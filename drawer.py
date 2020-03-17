#!/usr/bin/env python3

from tkinter import *

from random import randint, random
from parser import parseObj
from math import cos, sin, sqrt
from process import findNeighboursAndContours

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
    print(deltaY, deltaX, "ok")
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

def drawNormalVector(canvas, vertice1, vertice2, vertice3, color="blue", width=1, d=600):
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

    canvas.create_line(setCoordX(vmid[0]), setCoordY(vmid[1]),
                        setCoordX(vmid[0]) + vx, setCoordY(vmid[1]) - vy,
                        fill=color, width=width)

    lenU = sqrt(ux ** 2 + uy ** 2)
    lenV = sqrt(vx ** 2 + vy ** 2)
    b1x = ux * 5/lenU - vx * 5/lenV
    b1y = uy * 5/lenU - vy * 5/lenV

    canvas.create_line(setCoordX(vmid[0]) + vx, setCoordY(vmid[1]) - vy,
                        setCoordX(vmid[0]) + vx + b1x, setCoordY(vmid[1]) -vy - b1y,
                        fill=color, width=width)

    b2x = - ux * 5/lenU - vx * 5/lenV
    b2y = - uy * 5/lenU - vy * 5/lenV
    canvas.create_line(setCoordX(vmid[0]) + vx, setCoordY(vmid[1]) - vy,
                        setCoordX(vmid[0]) + vx + b2x, setCoordY(vmid[1]) - vy - b2y,
                        fill=color, width=width)

def setCoordX(curX):
    #return centreX + curX * scale
    return (curX - centreX) / scale

def setCoordY(curY):
    #return centreY - curY * scale
    return canvas_height - (curY - centreY) / scale

def drawObj(filename):

    master = Tk()

    vertices, textureCoord, normals, faces = parseObj(filename)

    initGlobalVariables(vertices)

    print(centreX, centreY, scale)

    neighbours, contours = findNeighboursAndContours(len(vertices), faces)


    canvas = Canvas(master,
           width=canvas_width,
           height=canvas_height)
    canvas.pack()

    canvas.create_rectangle(0, 0, canvas_width, canvas_height, fill= "white")

    print(setCoordX(40))

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
                if random() > 0:
                    drawNormalVector(canvas, vertices[i], vertices[j], vertices[contours[key]])
                drawLine(canvas, vertices[i], vertices[j], "green", 5)
            else:
                drawLine(canvas, vertices[i], vertices[j], "grey", 1)

    random_id = randint(0, len(vertices))
    for j in neighbours[random_id]:
        drawLine(canvas, vertices[random_id], vertices[j], "yellow", 2)

    #for vertice in vertices:
    #    drawCross(canvas, vertice, 10, random() * 3.14/2, "orange", 4)



    mainloop()

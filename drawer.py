#!/usr/bin/env python3

from tkinter import *

from random import randint
from parser import parseObj
from process import findNeighboursAndContours

def drawLine(canvas, vertice1, vertice2, centreX, centreY, color, width):
    canvas.create_line(centreX + vertice1[0], centreY - vertice1[1],
                    centreX + vertice2[0], centreY - vertice2[1],
                    fill=color, width=width)

def drawObj(filename):

    vertices, textureCoord, normals, faces = parseObj(filename)
    neighbours, contours = findNeighboursAndContours(len(vertices), faces)


    master = Tk()

    canvas_width = 400
    canvas_height = 400

    centreX = canvas_width // 2
    centreY = canvas_height // 2

    canvas = Canvas(master,
           width=canvas_width,
           height=canvas_height)
    canvas.pack()

    canvas.create_rectangle(0, 0, 400, 400, fill= "white")

    for face in faces:
        v1, v2, v3 = face[0][0] - 1, face[1][0] - 1, face[2][0] - 1
        canvas.create_polygon(centreX + vertices[v1][0], centreY - vertices[v1][1],
                            centreX + vertices[v2][0], centreY - vertices[v2][1],
                            centreX + vertices[v3][0], centreY - vertices[v3][1])

    for i in range(len(vertices)):
        for j in neighbours[i]:
            if j in contours[i]:
                drawLine(canvas, vertices[i], vertices[j], centreX, centreY, "red", 2)
            else:
                drawLine(canvas, vertices[i], vertices[j], centreX, centreY, "grey", 1)

    random_id = randint(0, len(vertices))
    for j in neighbours[random_id]:
        drawLine(canvas, vertices[random_id], vertices[j], centreX, centreY, "yellow", 2)



    mainloop()
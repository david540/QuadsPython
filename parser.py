#!/usr/bin/env python3

from itertools import permutations

def parseObj(filename):

    def parseFace(info):
        return [int(i) for i in info.split("/")]

    with open(filename, 'r') as objFile:
        content = objFile.read()
        lines = content.split("\n")

        vertices = []
        textureCoord = []
        normals = []
        faces = []

        for line in lines:
            infos = line.split(" ")
            entete = infos[0]
            if entete == 'o':
                print("Parsing of object :", infos[1])
            elif entete == 'v':
                vertices.append([float(infos[1]), float(infos[2]), float(infos[3])])
            elif entete == 'vt':
                textureCoord.append([float(infos[1]), float(infos[2])])
            elif entete == 'vn':
                normals.append([float(infos[1]), float(infos[2]), float(infos[3])])
            elif entete == 'f':
                faces.append([parseFace(infos[1])[0], parseFace(infos[2])[0], parseFace(infos[3])[0]])

        return vertices, textureCoord, normals, faces


def findNeighboursAndContours(vertices, faces):
    nbVertices = len(vertices)
    neighbours = [[] for i in range(nbVertices)]
    neifaces = [[] for i in range(len(faces))]
    normalsVect = [[0, 0] for i in range(nbVertices)]
    contours = {}
    facesIdentification = {}
    for i, face in enumerate(faces):
        for v1, v2, v3 in permutations([face[0] - 1, face[1] - 1, face[2] - 1], 3):
            key = (min(v1, v2), max(v1, v2))
            if key in facesIdentification and i not in facesIdentification[key]:
                facesIdentification[key].append(i)
                neifaces[facesIdentification[key][0]].append(facesIdentification[key][1])
                neifaces[facesIdentification[key][1]].append(facesIdentification[key][0])
            else:
                facesIdentification[key] = [i]
            if v2 not in neighbours[v1]:
                neighbours[v1].append(v2)
                contours[key] = [v3]
            else:
                if key in contours:
                    del(contours[key])



    for i, face in enumerate(faces):
        count = 0
        listBoundVertices = []
        commonVertice = None
        for v1, v2, v3 in permutations([face[0] - 1, face[1] - 1, face[2] - 1], 3):
            key = v1, v2
            if key in contours:
                count += 1
                if v1 not in listBoundVertices:
                    listBoundVertices.append(v1)
                else:
                    commonVertice = v1
                if v2 not in listBoundVertices:
                    listBoundVertices.append(v2)
                else:
                    commonVertice = v2
        if count > 3:
            print("ERROR 239959")
        elif count == 2 or count == 3:
            if len(listBoundVertices) != 3 or commonVertice == None:
                print("ERROR 33449")
            center = [sum([vertices[n - 1][j] for n in face]) / 3 for j in range(3)]
            vertices.append(center)
            neighbours.append([])
            vertices[-1], vertices[commonVertice] = vertices[commonVertice], vertices[-1]
            for v in listBoundVertices:
                if v == commonVertice:
                    continue
                del(contours[min(v, commonVertice), max(v, commonVertice)])
                contours[v, len(vertices) - 1] = commonVertice
                faces.append([v + 1, commonVertice + 1, len(vertices)])

            neifaces.append([])
            neifaces.append([])
            for f1, f2, f3 in permutations([len(faces) - 1, len(faces) - 2, i]):
                neifaces[f1].append(f2)
                neifaces[f1].append(f3)

            for n in listBoundVertices:
                neighbours[n].append(len(vertices) - 1)
                neighbours[-1].append(n)


    #print(len(faces))
    #print(neifaces)
    return neighbours, contours, facesIdentification, neifaces

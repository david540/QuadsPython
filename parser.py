#!/usr/bin/env python3

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
                faces.append([parseFace(infos[1]), parseFace(infos[2]), parseFace(infos[3])])

        return vertices, textureCoord, normals, faces

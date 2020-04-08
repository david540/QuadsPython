#!/usr/bin/env python3

from math import sqrt, acos, cos, sin, asin
import numpy as np
import scipy.sparse.linalg
from itertools import combinations, permutations
import cmath
from process import *


def scalarFieldLinearSystemDisjointSet(faces, du, dv, neifaces, ds, vertices, singularities):

    A = []
    b = []
    C = 1
    C_inter = 1
    #r = [faces[0][0] - 1, faces[220][0] - 1]
    r = faces[0][0] - 1
    #print(list(ds.itersets()))
    s = set()
    ids = []
    nbTriangles = len(faces)
    for t in range(nbTriangles):
        for i in range(len(faces[t])):
            s.add(ds.find(3 * t + i))
    for t in range(nbTriangles):
        for i in range(len(faces[t])):
            for idV, v in enumerate(s):
                if ds.find(3 * t + i) == v:
                    ids.append(idV)
                    break

    if len(ids) != 3 * len(faces):
        print("ERROR 30304")

    lenS = len(s)

    for t in range(nbTriangles):
        for k, m in permutations(range(len(faces[t])), 2):
            i = faces[t][k] - 1
            A.append([0 for _ in range(2 * lenS)])
            if i == r:
                A[-1][i] = C
                b.append([0])
                A.append([0 for _ in range(2 * lenS)])
                A[-1][lenS + i] = C
                b.append([0])
                continue
            j = faces[t][m] - 1
            z = complex(vertices[i][0] - vertices[j][0], vertices[i][1] - vertices[j][1])
            if z == 0:
                print("ERROR 49430")
                continue
            somZ = du[t] / (z/abs(z))
            #b.append([0])
            reference = ds.find(3 * t + m)
            A[-1][ids[3 * t + m]] = -1
            A[-1][ids[3 * t + k]] = 1
            b.append([sqrt(dist(vertices[i], vertices[j])) * somZ.real])

            A.append([0 for _ in range(2 * lenS)])
            somZ = dv[t] / (z/abs(z))
            A[-1][lenS + ids[3 * t + m]] = -1
            A[-1][lenS + ids[3 * t + k]] = 1
            b.append([sqrt(dist(vertices[i], vertices[j])) * somZ.real])

    #print("gestion des singularites")
    #for t, n, k, m, val in singularities:
    #    if val == 0:
    #        print("ERROR 430534")
    #    elif abs(val) == 1:
    #        idG, idD = ids[3 * t + k], ids[3 * n + m]
    #        if val == -1:
    #            idG, idD = idD, idG
    #        A.append([0 for _ in range(2 * len(s))])
    #        A[-1][idG] = C_inter
    #        A[-1][lenS + idD] = C_inter
    #        b.append([0])
    #        A.append([0 for _ in range(2 * len(s))])
    #        A[-1][lenS + idG] = C_inter
    #        A[-1][idD] = -C_inter
    #        b.append([0])
    #    elif abs(val) == 2:
    #        idG, idD = ids[3 * t + k], ids[3 * n + m]
    #        if val == -2:
    #            idG, idD = idD, idG
    #        A.append([0 for _ in range(2 * len(s))])
    #        A[-1][idG] = C_inter
    #        A[-1][idD] = C_inter
    #        b.append([0])
    #        A.append([0 for _ in range(2 * len(s))])
    #        A[-1][lenS + idG] = C_inter
    #        A[-1][lenS + idD] = C_inter
    #        b.append([0])

    print("resolution du système linéaire")
    A = np.matrix(A)
    b = np.matrix(b)

    AtA = A.transpose() * A
    Atb = A.transpose() * b

    X = np.linalg.inv(AtA) * Atb

    sol = X.transpose().tolist()[0]

    mini = min(sol)
    for i in range(len(sol)):
        sol[i] -= mini

    u_coins_triangles = [0 for i in range(3 * len(faces))]
    v_coins_triangles = [0 for i in range(3 * len(faces))]
    for t in range(len(faces)):
        for k in range(len(faces[t])):
            vref = ds.find(3 * t + k)
            for i, v in enumerate(s):
                if v == vref:
                    u_coins_triangles[3 * t + k] = sol[i]
                    v_coins_triangles[3 * t + k] = sol[i + lenS]
                    break
    print("champs u et v calculés")
    return u_coins_triangles, v_coins_triangles

def scalarFieldMiddleTriangle(faces, du, vertices):
    A = []
    b = []
    C = 1
    r = faces[4][0] - 1
    sol = [0 for i in range(len(vertices))]
    #u[j] = 0
    #duref = du[4]
    for t in range(len(faces)):
        for k, m in permutations(range(len(faces[t])), 2):
        #for k in range(len(faces[t])):
            i = faces[t][k] - 1
            A.append([0 for _ in range(len(vertices))])
            if i == r:
                A[-1][i] = C
                b.append([0])
                continue
            j = faces[t][m] - 1
            z = complex(vertices[i][0] - vertices[j][0], vertices[i][1] - vertices[j][1])
            if z == 0:
                print("ERROR 49430")
                continue
            somZ = du[t] / (z/abs(z))
            #A.append([0 for _ in range(len(vertices))])
            A[-1][i] = 1
            A[-1][j] = -1
            b.append([sqrt(dist(vertices[i], vertices[j])) * somZ.real])
            #sol[i] = sol[j] + sqrt(dist(vertices[i], vertices[j])) * somZ.real


    print(b[:3])
    A = np.matrix(A)
    b = np.matrix(b)

    AtA = A.transpose() * A
    Atb = A.transpose() * b

    X = np.linalg.inv(AtA) * Atb

    sol = X.transpose().tolist()[0]

    u = [0 for i in range(len(vertices))]
    mini = min(sol)
    for i in range(len(sol)):
        if i == j:
            continue
        #print(count[i])
        u[i] = sol[i] - mini

    return u

def naiveScalarFieldMiddleTriangle(faces, du, vertices):
    A = []
    b = []
    C = 1
    j = faces[4][0] - 1
    sol = [0 for i in range(len(vertices))]
    #u[j] = 0
    duref = du[4]
    for t in range(len(faces)):
        #for k, m in combinations(range(len(faces[t])), 2):
        for k in range(len(faces[t])):
            i = faces[t][k] - 1
            A.append([0 for _ in range(len(vertices))])
            if i == j:
                A[-1][j] = C
                b.append([0])
                continue
            #j = faces[t][m] - 1
            z = complex(vertices[i][0] - vertices[j][0], vertices[i][1] - vertices[j][1])
            if z == 0:
                print("ERROR 49430")
                continue
            somZ = duref / (z/abs(z))
            #A.append([0 for _ in range(len(vertices))])
            A[-1][i] = 1
            A[-1][j] = -1
            b.append([sqrt(dist(vertices[i], vertices[j])) * somZ.real])
            #sol[i] = sol[j] + sqrt(dist(vertices[i], vertices[j])) * somZ.real


    print(b[:3])
    A = np.matrix(A)
    b = np.matrix(b)

    AtA = A.transpose() * A
    Atb = A.transpose() * b

    X = np.linalg.inv(AtA) * Atb

    sol = X.transpose().tolist()[0]

    u = [0 for i in range(len(vertices))]
    mini = min(sol)
    for i in range(len(sol)):
        if i == j:
            continue
        #print(count[i])
        u[i] = sol[i] - mini

    return u

def scalarField(vertices, du, neighbours):
    length = len(vertices)
    u = [0 for i in range(length)]
    alreadyDone = [False for i in range(len(du))]
    idAct = 0
    compteur = 0
    idToTest = [0]
    while compteur < len(neighbours):
        idAct = idToTest[0]
        if alreadyDone[idAct] == True:
            del idToTest[0]
            continue
        alreadyDone[idAct] = True
        for n in neighbours[idAct]:
            z = complex(vertices[n][0] - vertices[idAct][0], vertices[n][1] - vertices[idAct][1])
            if z == 0:
                u[n] = 0
                continue
            somZ = du[idAct] / (z/abs(z))
            u[n] = u[idAct] + sqrt(dist(vertices[n], vertices[idAct])) * somZ.real
            idToTest.append(n)
        del idToTest[0]
        compteur+= 1

    mini = min(u)
    for i in range(length):
        u[i] -= mini
    return u

def scalarFieldLinearSystem(vertices, du, neighbours):
    length = len(vertices)
    A = []
    b = []
    alreadyDone = [False for i in range(len(du))]
    idAct = 0
    compteur = 0
    idToTest = [0]
    while compteur < len(neighbours):
        idAct = idToTest[0]
        if alreadyDone[idAct] == True:
            del idToTest[0]
            continue
        alreadyDone[idAct] = True
        for n in neighbours[idAct]:
            if alreadyDone[n] == True:
                continue
            z = complex(vertices[n][0] - vertices[idAct][0], vertices[n][1] - vertices[idAct][1])
            if z == 0:
                u[n] = 0
                continue
            somZ = du[idAct] / (z/abs(z))
            A.append([0 for _ in range(length)])
            A[-1][n] = 1
            A[-1][idAct] = -1

            b.append([sqrt(dist(vertices[n], vertices[idAct])) * somZ.real])
            #u[n] = u[idAct] + sqrt(dist(vertices[n], vertices[idAct])) * somZ.real
            idToTest.append(n)
        del idToTest[0]
        compteur+= 1

    A = np.matrix(A)
    b = np.matrix(b)

    AtA = A.transpose() * A
    Atb = A.transpose() * b

    AtA = np.matrix(AtA)
    Atb = np.matrix(Atb)

    X = np.linalg.inv(AtA) * Atb

    u = X.transpose().tolist()[0]

    mini = min(u)
    for i in range(length):
        u[i] -= mini
    return u

def scalarFieldAdaptive(vertices, du, neighbours):
    length = len(vertices)
    sizeSteps = 1
    u = [0 for i in range(length)]
    alreadyDone = [False for i in range(len(du))]
    idAct = 0
    compteur = 0
    idToTest = [0]
    while compteur < len(neighbours):
        #print(idToTest)
        idAct = idToTest[0]
        if alreadyDone[idAct] == True:
            del idToTest[0]
            continue
        alreadyDone[idAct] = True
        for n in neighbours[idAct]:
            if alreadyDone[n]:
                continue
            #alreadyDone[n] = True
            p = list(vertices[idAct])
            reachedNeigh = False
            uCur = 0
            countCur = 0
            signAdd = 0
            #print("\nnew")
            while not reachedNeigh:
                duAct = fieldOnAnyPoint(p, vertices, du, gap= 0.00001, alpha = 0.5)
                #print(duAct)
                z = complex(vertices[n][0] - p[0], vertices[n][1] - p[1])
                if z == 0:
                    print("ERROR 101")
                    u[n] = 0
                    continue
                somZ = duAct / (z/abs(z))
                distancetoH = np.dot([z.real, z.imag], [duAct.real, duAct.imag])
                if signAdd == 0 :
                    signAdd = np.sign(distancetoH)
                #print(distancetoH)
                if abs(distancetoH) <= sizeSteps + 0.001:
                    reachedNeigh = True
                    uCur += sqrt(dist(vertices[n], p)) * somZ.real #* (distancetoH / sizeSteps)
                    countCur += distancetoH / sizeSteps
                else:
                    uCur += sqrt(dist(vertices[n], p)) * somZ.real #* (distancetoH / sizeSteps)
                    countCur += 1
                    p[0] += duAct.real * sizeSteps * signAdd
                    p[1] += duAct.imag * sizeSteps * signAdd

            u[n] = u[idAct] + uCur #/ countCur
            idToTest.append(n)
        del idToTest[0]
        compteur+= 1

    mini = min(u)
    for i in range(length):
        u[i] -= mini
    return u

def naiveScalarField(vertices, du):
    n = len(vertices)
    u = [0 for i in range(n)]
    j = 4
    u[j] = 0
    for i in range(n):
        z = complex(vertices[i][0] - vertices[j][0], vertices[i][1] - vertices[j][1])
        if z == 0:
            u[i] = 0
            continue
        somZ = du[j] / (z/abs(z))
        u[i] = u[j] + sqrt(dist(vertices[i], vertices[j])) * somZ.real

    mini = min(u)
    for i in range(n):
        u[i] -= mini
    return u, j


def naiveScalarFieldMiddleTriangle(faces, du, vertices):
    n = len(vertices)
    u = [0 for i in range(n)]
    j = faces[4][0] - 1
    u[j] = 0
    duref = du[4]

    for t in range(len(faces)):
        for k, m in combinations(range(len(faces[t])), 2):
            i = faces[t][k] - 1
            z = complex(vertices[i][0] - vertices[j][0], vertices[i][1] - vertices[j][1])
            if z == 0:
                u[i] = 0
                continue
            somZ = duref / (z/abs(z))
            #print(somZ)
            u[i] = u[j] + sqrt(dist(vertices[i], vertices[j])) * somZ.real

    print(u)
    mini = min(u) - 1
    for i in range(n):
        u[i] -= mini
    return u

def lessnaiveScalarFieldMiddleTriangle(faces, du, vertices):
    n = len(vertices)
    u = [0 for i in range(n)]
    r = faces[4][0] - 1
    u[r] = 0
    duref = du[4]
    #print(du)
    for t in range(len(faces)):
        for k, m in combinations(range(len(faces[t])), 2):
            i = faces[t][k] - 1
            j = faces[t][m] - 1
            duact = du[t]
            if u[j] == 0:
                j = r
                duact = duref
            z = complex(vertices[i][0] - vertices[j][0], vertices[i][1] - vertices[j][1])
            if z == 0:
                u[i] = 0
                continue
            somZ = du[t] / (z/abs(z))
            #print(somZ)
            u[i] = u[j] + sqrt(dist(vertices[i], vertices[j])) * somZ.real

    print(u)
    mini = min(u) - 1
    for i in range(n):
        u[i] -= mini
    return u

def scalarFieldBad(vertices, du, x, y, xmin, ymin, exp = 3, eps = 0.01):
    #xi = vertices[:][0]
    #xmax = max(xi)
    #xmax = x
    somme = 0
    for i in range(len(vertices)):
        if x > vertices[i][0] + eps/2:
            somme += du[i] / (exp - 1) * (pow((eps/2 + abs(vertices[i][1] - y)), 1 - exp) - pow(vertices[i][0] - xmin + abs(vertices[i][1] - y), 1-exp)) \
            + eps * du[i] + du[i] / (exp - 1) * (pow((eps/2 + abs(vertices[i][1] - y)), 1 - exp) - pow((x - vertices[i][0] + abs(vertices[i][1] - y)), 1 - exp))
        elif x > vertices[i][0] - eps/2:
            somme += du[i] / (exp - 1) * (pow((eps/2 + abs(vertices[i][1] - y)), 1 - exp) - pow(vertices[i][0] - xmin + abs(vertices[i][1] - y), 1-exp)) \
            + (x - (vertices[i][0] - eps/2)) * du[i]
        else:
            somme += du[i] / (exp - 1) * (pow((vertices[i][0] - x + abs(vertices[i][1] - y)), 1 - exp) - pow(vertices[i][0] - xmin + abs(vertices[i][1] - y), 1-exp))

    return somme

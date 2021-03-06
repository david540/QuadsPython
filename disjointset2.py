import numpy as np
import cmath
from process import *
from itertools import permutations, combinations

class DisjointSetWithSign:

    # constructor knowing the number of elements
    def __init__(self, num):
        self.m_ids = [i for i in range(num)]
        self.m_signs = [True for i in range(num)]
#        print(num)

    #merge two component
    def merge(self, a, b, sameSign):

        assert (a>=0 and a<len(self.m_ids)), "Error : Merging component that is not in range(num)"
        assert (b>=0 and b<len(self.m_ids)), "Error : Merging component that is not in range(num)"

        rootA = self.findRoot(a)
        rootB = self.findRoot(b)

        if rootA == rootB:
            return

        if not self.m_signs[a]:
            sameSign = not sameSign
        if not self.m_signs[b]:
            sameSign = not sameSign



        self.m_ids[rootA] = self.m_ids[rootB]

        self.m_signs[rootA] = ((self.m_signs[rootA] == self.m_signs[rootB]) == sameSign)

    def getNumSets(self):

        for i in range(len(self.m_ids)):
            self.findRoot(i)

        #count the different roots

        num = 0

        for i in range(len(self.m_ids)):

            if i == self.m_ids[i]:
                num += 1

        return num

    def getSetsId(self):

        for i in range(len(self.m_ids)):
            self.findRoot(i)

        #prepare the correspondance root Id => setId

        idToSetId = [0 for i in range(len(self.m_ids))]
        numSets = 0

        for i in range(len(self.m_ids)):
            if i == self.m_ids[i]:
                idToSetId[i] = numSets
                numSets += 1

        for i in range(len(self.m_ids)):

            idToSetId[i] = idToSetId[self.m_ids[i]]

        return idToSetId

    #return the number of values
    def __len__(self):

        return len(self.m_ids)

    #return the parent indices
    def getParentIds(self):
        return self.m_ids

    #return the signs, note must be called/used after calling getSetsId and before doing any merge */
    def getSigns(self):
        return self.m_signs

    def findRoot(self, i):


        actSameSign = self.m_signs[i]
        #find the root
        rootId = i

        while rootId != self.m_ids[rootId]:
            actSameSign = (actSameSign==self.m_signs[self.m_ids[rootId]])
            rootId = self.m_ids[rootId]

        # connect all the path to root

        rootSameSign = actSameSign
        actSameSign = True

        while i != rootId:

            newI = self.m_ids[i]

            self.m_ids[i] = rootId

            newSameSign = actSameSign == self.m_signs[i]

            self.m_signs[i] = actSameSign==rootSameSign

            i=newI

            actSameSign = newSameSign

        return rootId

def setDisjointSet(vertices, faces, neifaces, a, b, arbreCouvrant, contours):

    nbTriangles = len(faces)
    nbCoinsTriangles = 3 * nbTriangles
    ds = DisjointSetWithSign(2 * nbCoinsTriangles)
    listeCasDecoupe = []
    positiveSide = []
    for t in range(nbTriangles):
        for n in neifaces[t]:
            ids = []
            listVId = findCommonVertices(faces[t], faces[n])
            key = tuple(sorted([listVId[0][0], listVId[1][0]]))
            if key not in arbreCouvrant:
                for c in range(2):
                    listVId[c].append(t)
                    listVId[c].append(n)
                    if t < n:
                        listeCasDecoupe.append(listVId[c])
                continue
            for i in range(len(faces[t])):
                for j in range(len(faces[n])):
                    if faces[t][i] - 1 == faces[n][j] - 1:
                        diffPhase = cmath.phase(a[t]) - cmath.phase(a[n])
                        while diffPhase > np.pi:
                            diffPhase -= 2 * np.pi
                        while diffPhase < -np.pi:
                            diffPhase+= 2 * np.pi
                    #    if abs(diffPhase) < np.pi/4:
                    #else:
                        val = 0
                        while diffPhase + val * np.pi/2 < -np.pi/4:
                            val += 1
                        while diffPhase + val * np.pi/2 > np.pi/4:
                            val -= 1
                        #singularities.append(tuple([t, n, i, j, val]))
                        ids.append([i, j, val])

            #if len(ids) != 2:
            #    print("ERROR 504304")
                #return
            #print(id1, id2, parentIds[id1], parentIds[id2], numCoord)

            for k in range(len(ids)):
                val = ids[k][2]
                idG, idD = 3 * t + ids[k][0], 3 * n + ids[k][1]
                if(faces[t][ids[k][0]] != faces[n][ids[k][1]]):
                    print("oh non pas comme ça")
                if val < 0:
                    idG, idD = idD, idG
                uG, uD, vG, vD = idG, idD, nbCoinsTriangles + idG, nbCoinsTriangles + idD
                #print(abs(val))

                if val == 0:
                    ds.merge(uG, uD, True)
                    ds.merge(vG, vD, True)
                    #print(uG, uD, ds.getParentIds()[uG], ds.getParentIds()[uD])

                elif abs(val) == 1:
                    #print("ok")
                    ds.merge(uG, vD, False)
                    ds.merge(vG, uD, True)
                elif abs(val) == 2:
                    #print('ok')
                    ds.merge(uG, uD, False)
                    ds.merge(vG, vD, False)

    tolerance = 0.3
    for t in range(nbTriangles):
        for k, m in combinations(range(len(faces[t])), 2):
            if m == k:
                continue
            v1, v2 = min(faces[t][k] - 1, faces[t][m] - 1), max(faces[t][k] - 1, faces[t][m] - 1)
            dist = sqrt((vertices[v2][0] - vertices[v1][0]) ** 2 + (vertices[v2][1] - vertices[v1][1]) ** 2)
            if (v1, v2) in contours:
                if abs(a[t].real * (vertices[v2][0] - vertices[v1][0]) / dist + a[t].imag * (vertices[v2][1] - vertices[v1][1]) / dist) < tolerance:
                    #aaa = 0
                    ds.merge(3 * t + k, 3 * t + m, True)
                elif abs(b[t].real * (vertices[v2][0] - vertices[v1][0]) / dist + b[t].imag * (vertices[v2][1] - vertices[v1][1]) / dist) < tolerance:
                    #aaa = 0
                    ds.merge(nbCoinsTriangles + 3 * t + k, nbCoinsTriangles + 3 * t + m, True)
                else:
                    print(t, a[t].real, vertices[v2][0] - vertices[v1][0])
                    print("ERROR 320492304")

    setsId = ds.getSetsId()
    print("fin ds")


    return ds, listeCasDecoupe


def setDisjointSetCoupe(vertices, faces, neifaces, arbreCouvrant, ds):
    nbTriangles = len(faces)
    nbCoinsTriangles = 3 * nbTriangles
    ds_coupe = DisjointSetWithSign(ds.getNumSets())
    alreadyDone = [False for t in range(nbTriangles)]
    setsId = ds.getSetsId()
    for numCoord in range(2):
        for t in range(nbTriangles):
            for n in neifaces[t]:
                listVId = findCommonVertices(faces[t], faces[n])
                key = tuple(sorted([listVId[0][0], listVId[1][0]]))
                if key not in arbreCouvrant:
                    for i, k, m in listVId:
                        ds_coupe.merge(setsId[numCoord * nbCoinsTriangles + 3 * t + k], setsId[numCoord * nbCoinsTriangles + 3 * n + m], False)
                        ds_coupe.getSetsId()

                        #print(ds_coupe.getSigns()[setsId[numCoord * nbCoinsTriangles + 3 * t + k]],
                        #ds_coupe.getSigns()[setsId[numCoord * nbCoinsTriangles + 3 * n + m]])
                    ds_coupe.merge(setsId[numCoord * nbCoinsTriangles + 3 * t + listVId[0][1]],
                    setsId[numCoord * nbCoinsTriangles + 3 * t + listVId[1][1]], True)
                    ds_coupe.getSetsId()
                    #print(ds_coupe.getSigns()[setsId[numCoord * nbCoinsTriangles + 3 * t + listVId[0][1]]],\
                    #ds_coupe.getSigns()[setsId[numCoord * nbCoinsTriangles + 3 * t + listVId[1][1]]])
    return ds_coupe

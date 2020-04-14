import numpy as np
import cmath

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

        self.connectToRoot(a)
        self.connectToRoot(b)

        self.m_ids[self.m_ids[a]] = self.m_ids[b]

        self.m_signs[self.m_ids[a]] = ((self.m_signs[a] == self.m_signs[b]) == sameSign)

    def getNumSets(self):

        for i in range(len(self.m_ids)):
            self.connectToRoot(i)

        #count the different roots

        num = 0

        for i in range(len(self.m_ids)):

            if i == self.m_ids[i]:
                num += 1

        return num

    def getSetsId(self):

        for i in range(len(self.m_ids)):
            self.connectToRoot(i)

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

    def connectToRoot(self, i):

        #find the root
        rootId = self.m_ids[i]

        while rootId != self.m_ids[rootId]:
            rootId = self.m_ids[rootId]

        # connect all the path to root

        while i != rootId:

            newI = self.m_ids[i]

            self.m_ids[i] = rootId

            self.m_signs[i] = self.m_signs[i] == self.m_signs[rootId]

            i=newI


def setDisjointSet(vertices, faces, neifaces, a):
    nbTriangles = len(faces)
    nbCoinsTriangles = 3 * nbTriangles
    ds = DisjointSetWithSign(2 * nbCoinsTriangles)
    for t in range(nbTriangles):
        for n in neifaces[t]:
            ids = []
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
                #val = 0
                idG, idD = 3 * t + ids[k][0], 3 * n + ids[k][1]
                if(faces[t][ids[k][0]] != faces[n][ids[k][1]]):
                    print("oh non pas comme ça")
                if val < 0:
                    idG, idD = idD, idG
                uG, uD, vG, vD = idG, idD, nbCoinsTriangles + idG, nbCoinsTriangles + idD
                if val == 0:
                    ds.merge(uG, uD, True)
                    ds.merge(vG, vD, True)
                    #print(uG, uD, ds.getParentIds()[uG], ds.getParentIds()[uD])

                if abs(val) == 1:
                    ds.merge(uG, vD, False)
                    ds.merge(vG, uD, True)

                if abs(val) == 2:
                    ds.merge(uG, uD, False)
                    ds.merge(vG, vD, False)

    print("fin ds")
    return ds

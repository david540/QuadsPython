from disjoint_set import DisjointSet
import numpy as np
import cmath

def setDisjointSet(vertices, faces, neifaces, a):
    ds = DisjointSet()
    singularities = []
    nbTriangles = len(faces)
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
            for k in range(len(ids)):
                val = ids[k][2]

                idG, idD = 3 * t + ids[k][0], 3 * n + ids[k][1]
                if val < 0:
                    idG, idD = idD, idG
                if val == 0:
                    ds.union(idG, idD)
                    ds.union(nbTriangles + idG, nbTriangles + idD)

                    ds.union(2 * nbTriangles + idG, 2 * nbTriangles + idD)
                    ds.union(3 * nbTriangles + idG, 3 * nbTriangles + idD)

                if abs(val) == 1:
                    ds.union(idG, 3 * nbTriangles + idD)
                    ds.union(nbTriangles + idG, 2 * nbTriangles + idD)

                    ds.union(2 * nbTriangles + idG, idD)
                    ds.union(3 * nbTriangles + idG, nbTriangles + idD)

                if abs(val) == 2:
                    ds.union(idG, nbTriangles + idD)
                    ds.union(nbTriangles + idG, idD)

                    ds.union(2 * nbTriangles + idG, 3 * nbTriangles + idD)
                    ds.union(3 * nbTriangles + idG, 2 * nbTriangles + idD)

    return ds, singularities

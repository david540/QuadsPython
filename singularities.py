#!/usr/bin/env python3

from math import sqrt, acos, cos, sin, asin
import numpy as np
import scipy.sparse.linalg
import cmath
from process import *

def computeSingularities(faces, dz):
    singularities = []
    for i, face in enumerate(faces):
        id1 = face[0][0] - 1
        angle_ref = cmath.phase(dz[id1]) / 4
        id2 = face[1][0] - 1
        angle_2 = cmath.phase(dz[id2]) / 4
        id3 = face[2][0] - 1
        angle_3 = cmath.phase(dz[id3]) / 4
        while (angle_2 - angle_ref) > np.pi/4:
            angle_2 -= np.pi/2
        while (angle_2 - angle_ref) < -np.pi/4:
            angle_2 += np.pi/2

        while (angle_3 - angle_ref) > np.pi/4:
            angle_3 -= np.pi/2
        while (angle_3 - angle_ref) < -np.pi/4:
            angle_3 += np.pi/2


        if abs(angle_2 - angle_3) > np.pi/4:
            singularities.append(i)
        #else:
            #print(angle_2 - angle_3)
    return singularities

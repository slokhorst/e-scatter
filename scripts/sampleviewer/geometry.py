import numpy as np


class Triangle():
    def __init__(self, v1, v2, v3, mat_i, mat_o):
        self.v1 = np.array(v1)
        self.v2 = np.array(v2)
        self.v3 = np.array(v3)
        e1 = self.v2 - self.v1
        e2 = self.v3 - self.v1
        normal = np.cross(e1, e2)
        self.n = normal/np.linalg.norm(normal)
        self.mat_i = mat_i
        self.mat_o = mat_o

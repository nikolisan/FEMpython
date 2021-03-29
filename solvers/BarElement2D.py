import numpy as np
from numpy.linalg import norm

class BarElement2D():
    def __init__(self, id, fromNode, toNode, area, stiff, nodes):
        self._id = id
        self._fromNode = fromNode
        self._toNode = toNode
        self._stiff = stiff
        self._area = area
        self._stiff = stiff
        self._nodes = nodes

    def __str__(self):
        return f'Element {self.elid}: F: {self.fromPoint} T: {self.toPoint} \tL: {self.length:5.2f} \tA: {self.area} \tE: {self.stiff:3.2E}'

    @property
    def elid(self):
        return self._id
    
    @property
    def area(self):
        return self._area
    
    @property
    def stiff(self):
        return self._stiff

    @property
    def fromPoint(self):
        return np.array(self._nodes[self._fromNode])

    @property
    def toPoint(self):
        return np.array(self._nodes[self._toNode])
    
    @property
    def length(self):
        return norm(self.toPoint+self.fromPoint)


if __name__ == '__main__':
    nodes = {0: [0,0], 1:[0,5]}
    e = BarElement2D(0, 0, 1, 1, 30e6, nodes)
    print(e)
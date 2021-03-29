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
        return f'Element {self.elid}: F: {self.fromPoint} T: {self.toPoint} L: {self.length:5.2f} A: {self.area} E: {self.stiff:3.2E} Dofs: {self.dofs}'

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

    @property
    def dofs(self):
        ''' returns:  list containing the index of the dof
            example:  fromNode = 0 -> u=0, v=1
                      toNode   = 2 -> u=4, v=5'''
        return [self._fromNode*2, self._fromNode*2+1, self._toNode*2, self._toNode*2+1]


if __name__ == '__main__':
    nodes = {0: [0,0], 1: [5,0], 2: [0,5]}
    e = BarElement2D(id=0, fromNode=1, toNode=2, area=1, stiff=30e6,nodes=nodes)
    print(e.dofs)
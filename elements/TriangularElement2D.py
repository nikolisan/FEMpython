import numpy as np
from scipy.linalg import inv
from math import sqrt


class TriangularElement2D():
    '''Wrapper class for 2D triangular elements (plane stress, plane strain)
        Parameters:
            id: (int) The id of the element
            points: np.array(3x2) Array of the node coordinates
            nodes: (list(int)) List containing the ids of the vertices
            thickness: (float) The thickness of the element
            stiffness: (float) Young's modulus of the material
            poisson: (float) Poisson's ration of the material
    '''
    def __init__(self, id, points, nodes, thickness, stiffness, poisson):
        self._id = id
        self._thickness = thickness
        self._stiffness = stiffness
        self._poisson = poisson
        self._points = points
        self._nodes = nodes
        self._displacements = None
        self._stress = None
        self._strain = None

    @property
    def id(self):
        return self._id

    @property
    def points(self):
        return self._points

    @property
    def nodes(self):
        return self._nodes

    @property
    def stiff(self):
        return self._stiffness

    @property
    def poisson(self):
        return self._poisson

    @property
    def thickness(self):
        return self._thickness

    @property
    def area(self):
        # points = np.array([ [xi, yi], [xj, yj], [xk, yk] ])
        xi = self.points[0,0]
        xj = self.points[1,0]
        xk = self.points[2,0]
        yi = self.points[0,1]
        yj = self.points[1,1]
        yk = self.points[2,1]
        side_a = sqrt((xj-xi)**2 + (yj-yi)**2)
        side_b = sqrt((xk-xj)**2 + (yk-yj)**2)
        side_c = sqrt((xk-xi)**2 + (yk-yi)**2)
        s = (side_a + side_b + side_c) / 2
        return sqrt(s*(s-side_a)*(s-side_b)*(s-side_c))

    @property
    def centroid(self):
        return self.points.mean(axis=0)

    @property
    def a_matrix(self):
        # relates strains at any points within the element to nodal displacements
        # points = np.array([ [xi, yi], [xj, yj], [xk, yk] ])
        xi = self.points[0,0]
        xj = self.points[1,0]
        xk = self.points[2,0]
        yi = self.points[0,1]
        yj = self.points[1,1]
        yk = self.points[2,1]
        A = np.array([
            [1, xi, yi, 0, 0, 0],
            [0, 0 ,0 , 1, xi, yi],
            [1, xj, yj, 0, 0, 0],
            [0, 0 ,0 , 1, xj, yj],
            [1, xk, yk, 0, 0, 0],
            [0, 0 ,0 , 1, xk, yk],
        ])
        return A
    
    @property
    def c_matrix(self):
        # relates strains at any points within the element to nodal displacements
        # points = np.array([ [xi, yi], [xj, yj], [xk, yk] ])
        C = np.array([
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 1],
            [0, 0, 1, 0, 1, 0]
        ])
        return C

    @property
    def b_matrix(self):
        return self.c_matrix.dot(inv(self.a_matrix))

    @property
    def dofs(self):
        dofs_list = []
        for node in self.nodes:
            dofs_list.append(node*2)
            dofs_list.append(node*2+1)
        return dofs_list

    @property
    def displacements(self):
        return self._displacements

    @displacements.setter
    def displacements(self, u):
        self._displacements = u

    @property
    def stress(self):
        return self._stress

    @stress.setter
    def stress(self, s):
        self._stress = s

    @property
    def strain(self):
        return self._strain

    @strain.setter
    def strain(self, s):
        self._strain = s



class TriangularPlaneStressElement(TriangularElement2D):
    '''2D Plane Stress Element
        Parameters:
            id: (int) The id of the element
            points: np.array(3x2) Array of the node coordinates
            thickness: (float) The thickness of the element
            stiffness: (float) Young's modulus of the material
            poisson: (float) Poisson's ration of the material
    '''
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def elasticity_matrix(self):
        Cd = self.stiff / (1-self.poisson**2)
        D = np.array([[1, self.poisson, 0], [self.poisson, 1, 0], [0, 0, ((1-self.poisson)/2)]])
        return Cd*D
        # return D

    @property
    def stiffness_matrix(self):
        B = self.b_matrix
        D = self.elasticity_matrix
        t = self.thickness
        A = self.area
        K = np.transpose(B).dot(D).dot(B)
        return t*A*K
        # return K



if __name__ == '__main__':
    points = np.array([[0,0],[2,0],[0,1]])
    # points = np.array([[2,0],[1.5,1.5],[0,1]])
    element = TriangularPlaneStressElement(1, points, [1, 3, 2], 0.1, 3e9, 0.333)
    print('AREA: ', element.area)
    print(element.nodes)
    print(element.centroid)
    with np.printoptions(precision=3, suppress=True):
        print(element.b_matrix.T)
        print(element.elasticity_matrix)
        print(element.b_matrix)
        print(element.stiffness_matrix)
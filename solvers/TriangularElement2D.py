import numpy as np
from math import sqrt


class TriangularElement2D():
    '''Wrapper class for 2D triangular elements (plane stress, plane strain)
        Parameters:
            id: (int) The id of the element
            points: np.array(3x2) Array of the node coordinates
            thickness: (float) The thickness of the element
            stiffness: (float) Young's modulus of the material
            poisson: (float) Poisson's ration of the material
    '''
    def __init__(self, id, points, thickness, stiffness, poisson):
        self._id = id
        self._thickness = thickness
        self._stiffness = stiffness
        self._poisson = poisson
        self._points = points

    @property
    def points(self):
        return self._points

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
        xk = self.points[1,0]
        yi = self.points[0,1]
        yj = self.points[1,1]
        yk = self.points[2,1]
        side_a = sqrt((xj-xi)**2 + (yj-yi)**2)
        side_b = sqrt((xk-xj)**2 + (yk-yj)**2)
        side_c = sqrt((xk-xi)**2 + (yk-yi)**2)
        s = (side_a + side_b + side_c) / 2
        return sqrt(s*(s-side_a)*(s-side_b)*(s-side_c))

    @property
    def b_matrix(self):
        # relates strains at any points within the element to nodal displacements
        # points = np.array([ [xi, yi], [xj, yj], [xk, yk] ])
        xi = self.points[0,0]
        xj = self.points[1,0]
        xk = self.points[1,0]
        yi = self.points[0,1]
        yj = self.points[1,1]
        yk = self.points[2,1]
        a1 = xk - xj
        b1 = yj - yk
        a2 = xi - xk
        b2 = yk - yi
        a3 = xj - xi
        b3 = yi - yj
        Cb = 1 / (2*self.area)
        B = np.array([[b1, 0, b2, 0, b3, 0], [0, a1, 0, a2, 0, a3], [a1, b1, a2, b2, a3, b3]])
        return Cb*B


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

    @property
    def stiffness_matrix(self):
        B = self.b_matrix
        D = self.elasticity_matrix
        t = self.thickness
        A = self.area
        K = np.transpose(B).dot(D).dot(B)
        return t*A*K


if __name__ == '__main__':
    points = np.array([[0,0],[1,0],[0,1]])
    element = TriangularPlaneStressElement(1, points, 0.1, 3e9, 0.333)
    print(element.stiffness_matrix)
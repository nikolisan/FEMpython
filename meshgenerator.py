import os
from loguru import logger
import numpy as np
import matplotlib.pyplot as plt


# TODO:
# [ ] Create dictionaries of all the mesh properties: elements, nodes
# [ ] Make plot non-blocking to the main thread
# [ ] Implement a way to predifine specific nodes

class UniformMeshQuad4():
    '''Create a uniform mesh of a rectangular domain

        Parameters:
            project_name: (str) The name of the current project
            lx, ly: (float) co-ordinates of the lower left corner of the domain
            ux, uy: (float) co-ordinates of the upper right corner of the domain
            nx, ny: number of elements in the x and y axis respectively
            plot: (bool) Show the mesh after generation (Default = True)
    '''
    def __init__(self, project_name, lx, ly, ux, uy, nx, ny, plot=True):
        self.project_name = project_name
        self.plot = plot
        self.lx = lx ; self.ly = ly
        self.ux = ux ; self.uy = uy
        self.nx = nx ; self.ny = ny
        self.check_files()

    def create_meshgrid(self, lx, ly, ux, uy, nx, ny):
        xx = np.linspace(lx, ux, nx+1)
        yy = np.linspace(ly, uy, ny+1)
        x, y = np.meshgrid(xx, yy)
        verts = []
        for yi in np.transpose(y)[0]:
            for xi in x[0]:
                verts.append([xi, yi])
        logger.success('Meshgrid generated')
        return (x, y, verts)

    def plot_mesh_quad4(self):
        '''Create plot a generated mesh'''
        nodes = self.nodes
        plt.plot(self.x, self.y, '--', color='slategray')
        plt.plot(np.transpose(self.x), np.transpose(self.y), '--', color='slategray')
        plt.plot(self.x, self.y, 'o', color='darkorange')
        for i, vert in enumerate(self.verts):
            plt.text(vert[0], vert[1], str(i), color='black', fontsize=10)
        for element in self.elements:
            elid = element[0]
            lln = nodes[element[1]]
            urn = nodes[element[3]]
            dx = (urn[1] - lln[1]) / 2
            dy = (urn[2] - lln[2]) / 2
            xpos = lln[1] + dx
            ypos = lln[2] + dy
            plt.text(xpos, ypos, str(elid), color='navy', fontsize=8)
        plt.show()
    
    def create_mesh_arrays(self):
        '''Generate the arrays holding the node ids, the nodes, and the elements

        nodes_id array: shape of the grid, each cell holds the id for each node
        nodes array: [ [node_id, x, y] ]
        elements arrat: [ [element_id, list(connected nodes_id)] ]
        '''
        node_ids = np.zeros_like(self.x).astype(int)
        nodes = np.zeros((self.nx+1)*(self.ny+1)*3).reshape(((self.nx+1)*(self.ny+1), 3))
        elements = np.zeros(self.nx*self.ny*5).reshape((self.nx*self.ny, 5)).astype(int)
        id = 0
        for i in range(self.ny+1):
            for j in range(self.nx+1):
                node_ids[i,j] = id
                nodes[id] = [id, self.x[i, j], self.y[i, j]]
                id += 1
        id = 0
        for i in range(self.ny):
            for j in range(self.nx):
                elements[id] = [id, node_ids[i,j], node_ids[i,j+1], node_ids[i+1,j+1], node_ids[i+1,j]]
                id += 1
        logger.success('Arrays created')
        return (node_ids, nodes, elements)

    def save_mesh_file(self):
        with open(os.path.join(os.getcwd(), self.project_name, 'mesh.dat'), 'w') as f:
            f.write(f'Project: {self.project_name}\n')
            f.write(f'----------------------------------------')
            f.write('\nLower left corner coordinates:\n')
            f.write(f'{self.lx:{4.7}}\t{self.ly:{4.7}}')
            f.write('\nUpper right corner coordinates:\n')
            f.write(f'{self.ux:{4.7}}\t{self.uy:{4.7}}')
            f.write('\nNumber of subdivisions in x, y:\n')
            f.write(f' {self.nx}\t {self.ny}')
            f.write(f'\n----------------------------------------\n')
    
    def save_arrays_to_text(self):
        try:
            self.save_mesh_file()
            with open(os.path.join(os.getcwd(), self.project_name, 'mesh.dat'), 'a+') as f:
                f.write("Node ids:\n")
                for node_id in self.node_ids:
                    for id in node_id:
                        f.write(f'{id}\t')
                    f.write('\n')
                f.write("Nodes:\n")
                for node in self.nodes:
                    f.write(f'{node[0]:{3}.{5}} \t {node[1]:{3}.{5}} \t {node[2]:{3}.{5}}')
                    f.write('\n')
                f.write("Elements:\n")
                for element in self.elements:
                    f.write(f'{element[0]}\t{element[1]}\t{element[2]}\t{element[3]}\t{element[4]}')
                    f.write('\n')
            logger.success('Mesh data saved')
        except Exception as e:
            logger.error(f'Cannot save the file\n{e}')

    def generate_mesh(self):
        '''Genarates the meshgrid and the appropriate arrays'''
        try:
            self._x, self._y, self._verts = self.create_meshgrid(self.lx, self.ly, self.ux, self.uy, self.nx, self.ny)
            self._node_ids, self._nodes, self._elements = self.create_mesh_arrays()
            self.save_arrays_to_text()
            if self.plot:
                self.plot_mesh_quad4()
        except Exception as e:
            logger.error(f'Cannot generate mesh\n{e}')

    def check_files(self):
        project_path = os.path.join(os.getcwd(), self.project_name)
        if not os.path.isdir(project_path):
            os.mkdir(project_path)
            logger.info(f'Created folder {project_path}')

    @property
    def x(self):
        return self._x
    
    @property
    def y(self):
        return self._y
    
    @property
    def verts(self):
        return self._verts

    @property
    def node_ids(self):
        return self._node_ids
    
    @property
    def nodes(self):
        return self._nodes
    
    @property
    def elements(self):
        return self._elements


if __name__ == '__main__' :
    project_name = 'test-project'
    nx = 4 ; ny = 5
    lx = 0.0 ; ly = 0.0
    ux = 4.0 ; uy = 5.0
    
    mesh = UniformMeshQuad4(project_name, lx, ly, ux, uy, nx, ny, plot=False)
    mesh.generate_mesh()
    mesh.plot_mesh_quad4()
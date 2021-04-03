import triangle as tr
import numpy as np
import matplotlib.pyplot as plt


class TriangularMesh():
    '''
    model = {nodes, boundary_nodes, segments}
    '''
    def __init__(self, model:dict, options='pq0D', show_plot=False):
        self._model = model
        self._options = options
        self.show_plot = show_plot
        self.create_mesh()

    def set_points(self):
        self._points = np.asarray([node for node in self._model['nodes'].values()])
        return self._points

    def set_boundary_nodes(self):
        self._boundary_nodes = np.asarray([b for b in self._model['boundary_nodes'].values()])
        return self._boundary_nodes

    def set_segments(self):
        self._segments = np.asarray([seg for seg in self._model['segments'].values()])
        return self._segments

    def set_holes(self):
        self._holes = np.asarray([hole for hole in self._model['holes'].values()])
        return self._holes

    def export_mesh(self):
        model = self._model
        t = self._t
        new_model = {}
        nodes = {}
        elements = {}
        # TODO: check if the nodes in the original model have the same index
        #       as with the mesh
        for i, vert in enumerate(t['vertices'].tolist()):
            nodes[i] = vert
        for i, tri in enumerate(t['triangles'].tolist()):
            elements[i] = tri

        forces = model['forces']
        rdofs = model['rdofs']
        for i in range(len(model['nodes']), len(nodes)):
            forces[i] = [0, 0]
            rdofs[i] = [0, 0]

        new_model['name'] = model['name']
        new_model['nodes'] = nodes
        new_model['elements'] = elements
        new_model['rdofs'] = rdofs
        new_model['forces'] = forces
        new_model['stiff'] = model['stiff']
        new_model['poisson'] = model['poisson']
        new_model['thick'] = model['thick']
        new_model['ndofs'] = 2*len(nodes)

        return new_model

    def create_mesh(self):
        self.set_points()
        self.set_boundary_nodes()
        self.set_segments()
        if 'holes' in self._model:
            self.set_holes()
            geometry = {'vertices': self._points, 'vertex_markers': self._boundary_nodes, 'segments': self._segments, 'holes': self._holes}
        else:
            geometry = {'vertices': self._points, 'vertex_markers': self._boundary_nodes, 'segments': self._segments}
        
        self._t = tr.triangulate(geometry, self._options)
        self._mesh = self.export_mesh()
        tr.compare(plt, geometry, self._t)
        if self.show_plot:
            plt.show()
        return self._mesh

    @property
    def model(self):
        return self._mesh


if __name__ == '__main__':
    def set_geo_data():
        name = 'Test Model'
        # Nodal coordinates in (m)
        nodes = {
            0: [0, 0],
            1: [2, 0],
            2: [0, 1],
            3: [1.5, 1.5],
            # Hole nodes
            4: [1, 0.5],
            5: [1.5, 0.5],
            6: [1.5, 1],
            7: [1, 1]
        }
        boundary_nodes = {
            0: [1],
            1: [1],
            2: [1],
            3: [1],
            # Hole nodes have boundaries
            4: [1],
            5: [1],
            6: [1],
            7: [1]
        }
        segments = {
            # Segments are the lines which define the boundaries
            # It must form closed loops
            # id (int): [from_node_id, to_node_id] (int)
            # --- Main Region ---
            0: [0, 1],
            1: [1, 3],
            2: [3, 2],
            3: [2, 0],
            # ---- Hole ----
            4: [4, 5],
            5: [5, 6],
            6: [6, 7],
            7: [7, 4]
        }
        holes = {
            # If the model contains holes define a point for each hole
            0: [1.25, 0.75]
        }
        # Dictionary containing the supports. 1: restrained dof, 0: free dof (x, y)
        restrained_dofs = {
            0: [1, 1],
            1: [0, 0],
            2: [1, 1],
            3: [0, 0],
            4: [0, 0],
            5: [0, 0],
            6: [0, 0],
            7: [0, 0]
        }
        # Dictionary containing the external forces. 2 for each node [Fx, Fy]. (N)
        forces = {
            0: [0, 0],
            1: [-30e3, -40e3],
            2: [0, 0],
            3: [0, -30e3],
            4: [0, 0],
            5: [0, 0],
            6: [0, 0],
            7: [0, 0]
        }
        # Total number of degrees of freedom. Two for each node.
        ndofs = 2*len(nodes)
        # Young's modulus of the plate material (assume same material for all the elements)(N/m2)
        E = 3e9
        # Poisson's ratio of the plate material (assume same material for all the elements)
        nu = 1/3
        # Plate thickness (m)
        t = 0.1

        assert len(nodes) == len(restrained_dofs) == len(forces) == len(boundary_nodes)

        return {"name": name, "nodes": nodes, "boundary_nodes": boundary_nodes, "segments": segments, "holes": holes, "rdofs": restrained_dofs, "forces": forces, "stiff": E, "poisson": nu, "thick": t, "ndofs":ndofs}
    
    geo = set_geo_data()
    mesh = TriangularMesh(geo, 'pqa0.05D', show_plot=True)

    model = mesh.model
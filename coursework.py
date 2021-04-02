from fempy.solvers.PlaneTriangular2D import PlaneTriangular2D as Solver
from fempy.mesh.TriangularMesh import TriangularMesh
from matplotlib import pyplot as plt

def set_geo_data():
    name = 'Test Model'
    # Nodal coordinates in (m)
    nodes = {
        0: [0, 0],
        1: [2, 0],
        2: [2, 1],
        3: [1, 1],
        4: [0, 1]
    }
    boundary_nodes = {
        0: [1],
        1: [1],
        2: [1],
        3: [1],
        4: [1]
    }
    segments = {
        # Segments are the lines which define the boundaries
        # It must form closed loops
        # id (int): [from_node_id, to_node_id] (int)
        # --- Main Region ---
        0: [0, 1],
        1: [1, 2],
        2: [2, 3],
        3: [3, 4],
        4: [4, 0]
    }
    # Dictionary containing the supports. 1: restrained dof, 0: free dof (x, y)
    restrained_dofs = {
        0: [1, 1],
        1: [0, 1],
        2: [0, 0],
        3: [0, 0],
        4: [0, 0]
    }
    # Dictionary containing the external forces. 2 for each node [Fx, Fy]. (N)
    forces = {
        0: [0, 0],
        1: [0, 0],
        2: [0, 0],
        3: [0, -2e6],
        4: [0, 0]
    }
    # Total number of degrees of freedom. Two for each node.
    ndofs = 2*len(nodes)
    # Young's modulus of the plate material (assume same material for all the elements)(N/m2)
    E = 100e9
    # Poisson's ratio of the plate material (assume same material for all the elements)
    nu = 0.3
    # Plate thickness (m)
    t = 0.1

    assert len(nodes) == len(restrained_dofs) == len(forces) == len(boundary_nodes)

    return {"name": name, "nodes": nodes, "boundary_nodes": boundary_nodes, "segments": segments, "rdofs": restrained_dofs, "forces": forces, "stiff": E, "poisson": nu, "thick": t, "ndofs":ndofs}

geo = set_geo_data()
mesh = TriangularMesh(geo, 'pqa0.001')
model = mesh.model
solver = Solver(model, 'Stress')
solver.plot_model(show_nodes=False, show_text=False)
solver.solve()
solver.plot_stresses(show_nodes=False, show_text=False)
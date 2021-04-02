from fempy.solvers.PlaneTriangular2D import PlaneTriangular2D as Solver
from fempy.mesh.TriangularMesh import TriangularMesh
from matplotlib import pyplot as plt

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
mesh = TriangularMesh(geo, 'pa.001')
model = mesh.model
solver = Solver(model, 'Stress')
solver.plot_model(show_nodes=False, show_text=False)
solver.solve()
solver.plot_stresses(show_nodes=False, show_text=False)
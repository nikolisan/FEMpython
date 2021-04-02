from fempy.solvers.PlaneTriangular2D import PlaneTriangular2D as Solver

def set_geo_data():
    name = 'Test Model'
    # Nodal coordinates in (m)
    nodes = {
        0: [0, 0],
        1: [2, 0],
        2: [0, 1],
        3: [1.5, 1.5]
    }
    # Each element list contains the index of its nodes
    elements = {
        0: [0, 1, 2],
        1: [1, 3, 2]
    }
    # Dictionary containing the supports. 1: restrained dof, 0: free dof (x, y)
    restrained_dofs = {
        0: [1, 1],
        1: [0, 0],
        2: [1, 1],
        3: [0, 0]
    }
    # Dictionary containing the external forces. 2 for each node [Fx, Fy]. (N)
    forces = {
        0: [0, 0],
        1: [-30e3, -40e3],
        2: [0, 0],
        3: [0, -30e3]
    }
    # Total number of degrees of freedom. Two for each node.
    ndofs = 2*len(nodes)
    # Young's modulus of the plate material (assume same material for all the elements)(N/m2)
    E = 3e9
    # Poisson's ratio of the plate material (assume same material for all the elements)
    nu = 1/3
    # Plate thickness (m)
    t = 0.1

    assert len(nodes) == len(restrained_dofs) == len(forces)

    return {"name": name, "nodes": nodes, "elements": elements, "rdofs": restrained_dofs, "forces": forces, "stiff": E, "poisson": nu, "thick": t, "ndofs":ndofs}

model = set_geo_data()
solver = Solver(model, 'Stress')
solver.plot_model()
solver.solve()
solver.print_results()
solver.plot_stresses()

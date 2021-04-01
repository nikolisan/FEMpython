import numpy as np
from scipy.linalg import inv
import matplotlib.pyplot as plt
from loguru import logger

from elements.TriangularElement2D import TriangularPlaneStressElement as Element


#TODO:
# - Set displacements to element object
# - Set stresses to element object
# - Set strains to element object

def set_geo_data():
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

    return {"nodes": nodes, "elements": elements, "rdofs": restrained_dofs, "forces": forces, "stiff": E, "poisson": nu, "thick": t, "ndofs":ndofs}

def flatten_dofs(dofs:dict, restrained: bool):
    out = []
    for key in dofs.keys():
        if dofs[key][0] == restrained:
            out.append(key*2)
        if dofs[key][1] == restrained:
            out.append(key*2 + 1)
    return out

def plot_nodes(nodes):
    fig, ax = plt.subplots()
    x = [val[0] for val in nodes.values()]
    y = [val[1] for val in nodes.values()]
    # Plot nodes and node ids
    ax.scatter(x, y, s=100, color='darkorange', zorder=4)
    for i, pos in enumerate(zip(x,y)):
        ax.text(pos[0]-0.04, pos[1]-0.05, str(i), color='black', fontsize=8, zorder=5)

def create_elements(model):
    # el = Element(id, points, thickness, stiffness, poisson)
    elements = model['elements']
    nodes = model['nodes']
    thick = model['thick']
    stiff = model['stiff']
    poisson = model['poisson']

    new_elements = {}
    for id, element in enumerate(elements.values()):
        points = np.array([nodes[element[0]], nodes[element[1]], nodes[element[2]]], dtype='float')
        el = Element(id, points, element, thick, stiff, poisson)
        new_elements[id] = el
    
    model['elements'] = new_elements
    return model

def create_global_stiffness_matrix(nodes: dict, elements: Element, ndofs: int):
    logger.info('Creating global stiffness matrix')
    K = np.zeros([ndofs, ndofs])
    for element in elements.values():
        logger.info(f'Processing Element: {element.id}')
        ke = element.stiffness_matrix
        dofs = element.dofs
        for i in range(ke.shape[0]):
            for j in range(ke.shape[1]):
                K[dofs[i], dofs[j]] = K[dofs[i], dofs[j]] + ke[i, j]
    return K

def create_force_vector(forces, ndofs):
    F = []
    for force in forces.values():
        F.extend(force)
    return np.array(F)

def apply_boundary_conditions(K, F, restrained_dofs):
    logger.info('Applying boundary conditions')
    rdofs = flatten_dofs(restrained_dofs, True)
    rdofs = np.array(rdofs)
    K = np.delete(K, rdofs, axis=0)
    K = np.delete(K, rdofs, axis=1)
    F = np.delete(F,rdofs)

    return K, F

def calculate_displacements(K, F):
    logger.info('Calculating displacements')
    return inv(K).dot(F)

def create_global_displacement_matrix(restrained_dofs, ndofs, U):
    fdofs = flatten_dofs(restrained_dofs, False)
    Ug = np.zeros(ndofs)
    for i in range(len(fdofs)):
        Ug[fdofs[i]] = U[i]
    return Ug

def calculate_stress(Ug, elements, ndofs):
    dofs_list = [dof for dof in range(ndofs)]
    
    for element in elements.values():
        Ug_ = np.ma.array(Ug, mask=False)
        mask = list(set(dofs_list).difference(element.dofs))
        Ug_.mask[mask] = True
        print(mask)
        print(Ug_)
        stress = element.elasticity_matrix.dot(element.b_matrix).dot(Ug_.compressed())
        print(stress)

def print_results(U):
    print('\nResults:')
    with np.printoptions(precision=5, suppress=True):
        print('Nodal Displacements:\n', U)


def main():
    # Assume plane stress problem (t << w, h)
    model = set_geo_data()
    plot_nodes(model['nodes'])
    create_elements(model)

    K = create_global_stiffness_matrix(model['nodes'], model['elements'], model['ndofs'])
    F = create_force_vector(model['forces'], model['ndofs'])
    Kc, Fc = apply_boundary_conditions(K, F, model['rdofs'])
    U = calculate_displacements(Kc, Fc)
    Ug = create_global_displacement_matrix(model['rdofs'], model['ndofs'], U)
    
    calculate_stress(Ug, model['elements'], model['ndofs'])

    print_results(Ug)

    print_results(U)
    plt.show()
    
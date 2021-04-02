import numpy as np
from scipy.linalg import inv
import matplotlib.pyplot as plt
from loguru import logger

from elements.TriangularElement2D import TriangularPlaneStressElement as Element


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
    try:
        out = []
        for key in dofs.keys():
            if dofs[key][0] == restrained:
                out.append(key*2)
            if dofs[key][1] == restrained:
                out.append(key*2 + 1)
        return out
    except Exception as e:
        return None

def plot_nodes(nodes, elements):
    try:
        fig, ax = plt.subplots()
        x = [val[0] for val in nodes.values()]
        y = [val[1] for val in nodes.values()]
        tri = [el.nodes for el in elements.values()]
        # Plot nodes and node ids
        ax.scatter(x, y, s=100, color='darkorange', zorder=4)
        for i, pos in enumerate(zip(x,y)):
            ax.text(pos[0], pos[1], str(i), color='black', ha='center', va='center', fontsize=8, zorder=5)
        ax.triplot(x, y, triangles=tri, color='slategray')
        for el in elements.values():
            ax.text(el.centroid[0], el.centroid[1], str(el.id), color='black', ha='center', va='center', fontsize=10, zorder=5)
    except Exception as e:
        logger.error(e)
        return None

def create_elements(model):
    try:
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
    except Exception as e:
        logger.error(e)
        return None

def create_global_stiffness_matrix(nodes: dict, elements: Element, ndofs: int):
    try:
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
    except Exception as e:
        logger.error(e)
        return None

def create_force_vector(forces, ndofs):
    try:
        F = []
        for force in forces.values():
            F.extend(force)
        return np.array(F)
    except Exception as e:
        logger.error(e)
        return None

def apply_boundary_conditions(K, F, restrained_dofs):
    try:
        logger.info('Applying boundary conditions')
        rdofs = flatten_dofs(restrained_dofs, True)
        rdofs = np.array(rdofs)
        K = np.delete(K, rdofs, axis=0)
        K = np.delete(K, rdofs, axis=1)
        F = np.delete(F,rdofs)
        return K, F
    except Exception as e:
        logger.error(e)
        return None

def calculate_displacements(K, F):
    try:
        logger.info('Calculating displacements')
        return inv(K).dot(F)
    except Exception as e:
        logger.error(e)
        print(None)

def create_global_displacement_matrix(restrained_dofs, ndofs, U):
    try:
        fdofs = flatten_dofs(restrained_dofs, False)
        Ug = np.zeros(ndofs)
        for i in range(len(fdofs)):
            Ug[fdofs[i]] = U[i]
        return Ug
    except Exception as e:
        logger.error(e)
        return None

def store_displ_on_element(Ug, elements):
    try:
        for element in elements.values():
            Ue = Ug[element.dofs]
            element.displacements = Ue
    except Exception as e:
        logger.error(e)

def calculate_stress(elements):
    try:
        logger.info('Calculating stress')
        for element in elements.values():
            stress = element.elasticity_matrix.dot(element.b_matrix).dot(element.displacements)
            element.stress = stress
    except Exception as e:
        logger.error(e)

def calculate_strain(elements):
    try:
        logger.info('Calculating strain')
        for element in elements.values():
            strain = inv(element.elasticity_matrix).dot(element.stress)
            element.strain = strain
    except Exception as e:
        logger.error(e)

def print_results(Ug, nodes, ndofs, elements):
    print('\nResults:')
    print('Nodal Displacements:')
    for node in range(len(nodes)):
        print(f' - Node: {node}')
        print(f'\t u_{node}: {Ug[node*2]:+.5f}'.replace('+', ' '))
        print(f'\t v_{node}: {Ug[node*2+1]:+.5f}'.replace('+', ' '))
    print('\nElement Stress:')
    for element in elements.values():
        print(f' - Element: {element.id}')
        print(f'\t σ_x : {element.stress[0]:+.3e}'.replace('+', ' '))
        print(f'\t σ_y : {element.stress[1]:+.3e}'.replace('+', ' '))
        print(f'\t τ_xy: {element.stress[2]:+.3e}'.replace('+', ' '))
    print('\nElement Strain:')
    for element in elements.values():
        print(f' - Element: {element.id}')
        print(f'\t ε_x : {element.strain[0]:+.3e}'.replace('+', ' '))
        print(f'\t ε_y : {element.strain[1]:+.3e}'.replace('+', ' '))
        print(f'\t γ_xy: {element.strain[2]:+.3e}'.replace('+', ' '))


def main():
    # Assume plane stress problem (t << w, h)
    model = set_geo_data()
    
    create_elements(model)
    plot_nodes(model['nodes'], model['elements'])

    K = create_global_stiffness_matrix(model['nodes'], model['elements'], model['ndofs'])
    F = create_force_vector(model['forces'], model['ndofs'])
    Kc, Fc = apply_boundary_conditions(K, F, model['rdofs'])
    U = calculate_displacements(Kc, Fc)
    Ug = create_global_displacement_matrix(model['rdofs'], model['ndofs'], U)
    
    store_displ_on_element(Ug, model['elements'])
    calculate_stress(model['elements'])
    calculate_strain(model['elements'])

    print_results(Ug, model['nodes'], model['ndofs'], model['elements'])

    plt.show()
    
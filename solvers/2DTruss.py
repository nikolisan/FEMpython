import numpy as np
from numpy.linalg import norm
from scipy.linalg import inv
import matplotlib.pyplot as plt
from loguru import logger

from elements.BarElement2D import BarElement2D
    

def set_geo_data():
    nodes = {0: [0,0], 1: [10,5], 2: [0,10]}
    elements = {0: [0, 1], 1: [1, 2]}
    ext_forces = {0: [0, 0], 1: [0, -200], 2: [0, 0]}
    # For each node we have u, v dof
    dofs = {0: [1, 1], 1: [1, 1], 2: [1, 1]}
    # The restrained dofs are the 1 and the free dofs are 0 (u,v)
    restrained_dofs = {0: [1, 1], 1: [0, 0], 2: [1, 1]}
    areas = {0: 1, 1: 2}
    stiffs = {0: 30e6, 1: 30e6}
    ndofs = len(nodes)*2
    return {
        'nodes': nodes,
        'elements': elements,
        'ext_forces': ext_forces,
        'dofs': dofs,
        'restrained_dofs': restrained_dofs,
        'areas': areas,
        'stiffs': stiffs,
        'ndofs': ndofs
    }

def flatten_dofs(dofs:dict, restrained: bool):
    out = []
    for key in dofs.keys():
        if dofs[key][0] == restrained:
            out.append(key*2)
        if dofs[key][1] == restrained:
            out.append(key*2 + 1)
    return out

def create_elements(model):
    elements = model['elements']
    nodes = model['nodes']
    areas = model['areas']
    stiffs = model['stiffs']
    newElements = {}
    for id_ in elements.keys():
        newElement = BarElement2D(id_, elements[id_][0], elements[id_][1], areas[id_], stiffs[id_], nodes)
        newElements[id_] = newElement
    return newElements

def plot_system(nodes, elements, forces):
    fig, ax = plt.subplots()
    x = [val[0] for val in nodes.values()]
    y = [val[1] for val in nodes.values()]
    # Plot nodes and node ids
    ax.scatter(x, y, s=100, color='darkorange', zorder=4)
    for i, pos in enumerate(zip(x,y)):
        ax.text(pos[0]-0.04, pos[1]-0.05, str(i), color='black', fontsize=8, zorder=5)
    # Plot Elements and element ids
    for e in elements.values():
        ax.plot((e.fromPoint[0], e.toPoint[0]), (e.fromPoint[1], e.toPoint[1]), '-', color='slategray', zorder=2)
        idpos = abs((e.fromPoint - e.toPoint) / 2) + 0.04
        ax.text(idpos[0], idpos[1], str(e.elid), color='black', fontsize=8, zorder=5)
    # Plot forces
    # maxforce = max(forces[max(abs(forces, key=lambda key: forces[key]))][0], forces[max(abs(forces, key=lambda key: forces[key]))][1])
    maxforce = 200
    for fid in forces.keys():
        node = nodes[fid]
        fx = forces[fid][0]
        fy = forces[fid][1]
        if fx != 0 or fy != 0:
            ax.arrow(node[0], node[1], fx/round(maxforce*10/ax.get_xlim()[1]), fy/round(maxforce*10/ax.get_ylim()[1]), head_width=0.08, zorder=3)

def create_transformation_matrix(element):
    # Find the angle of the element's axis and the global x-axis
    element_vector = element.toPoint - element.fromPoint
    x_axis = np.array([1,0])
    y_axis = np.array([0,1])
    cosphi = np.dot(element_vector, x_axis) / norm(element_vector)
    sinphi = np.dot(element_vector, y_axis) / norm(element_vector)
    # # Create the transformation matrix
    transformMatrix = np.array([[cosphi, sinphi, 0, 0], [-sinphi, cosphi, 0, 0], [0, 0, cosphi, sinphi], [0, 0, -sinphi, cosphi]])
    return transformMatrix

def create_global_stiffness_matrix(nodes: dict, elements: BarElement2D, forces: dict, areas: dict, ndofs: int):
    logger.info('Creating global stiffness matrix')
    K = np.zeros([ndofs, ndofs])
    for element in elements.values():
        logger.info(f'Processing: {element}')
        L = element.length
        A = element.area
        E = element.stiff
        dofs = element.dofs
        # Construct the local stiffness matrix
        k_local = ( E * A / L ) * np.array([[1, 0, -1, 0], [0, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 0]])
        # Calculate the transformation matrix
        transformMatrix = create_transformation_matrix(element)
        # Derive the element's stiffness matrix to the global coordinate system
        k_global = np.transpose(transformMatrix).dot(k_local).dot(transformMatrix)
        
        # Assemble the global stiffness matrix by position the element's k matrix
        for i in range(k_global.shape[0]):
            for j in range(k_global.shape[1]):
                K[dofs[i], dofs[j]] = K[dofs[i], dofs[j]] + k_global[i, j]
    return K

def create_ext_force_vector(forces, ndofs):
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

def construct_global_displ_matrix(restrained_dofs, ndofs, U):
    fdofs = flatten_dofs(restrained_dofs, False)
    Ug = np.zeros(ndofs)
    for i in range(len(fdofs)):
        Ug[fdofs[i]] = U[i]
    return Ug

def calculate_reaction_forces(restrained_dofs, ndofs, U, K):
    logger.info('Calculating the forces on the structure (External, Reactions)')
    Ug = construct_global_displ_matrix(restrained_dofs, ndofs, U)
    Fr = K.dot(Ug)
    return Fr

def calculate_strain_stress(elements, Ug):
    logger.info('Calculating Strain and Stress')
    stresses = {}
    internalForces = {}
    strains = {}
    for element in elements.values():
        logger.info(f'Processing: Element {element.elid}')
        E = element.stiff
        L = element.length
        A = element.area
        k_local = ( E * A / L ) * np.array([[1, 0, -1, 0], [0, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 0]])
        transMatrix = create_transformation_matrix(element)
        # Keep from the global displacement matrix only the rows with the element's dofs
        Ueg = Ug[ element.dofs ]
        q = transMatrix.dot(Ueg)
        # calculate the local forces
        internalF = k_local.dot(q)
        # Remove rows with local v (as we have displacement only in local u)
        q = np.delete(q, [1, 3])
        internalF = np.delete(internalF, [1,3])
        strain = (q[1] - q[0]) / element.length
        stress = E*strain

        strains[element.elid] = strain
        stresses[element.elid] = stress
        internalForces[element.elid] = internalF
    
    return stresses, strains, internalForces
        
def print_results(U, F, stresess, strains):
    print('\nResults:')
    with np.printoptions(precision=4):
        print('Nodal Displacements:\n', U)
    with np.printoptions(precision=2, suppress=True):
        print('\nExternal and Reactions Forces:\n', F)
    print('\nStresses:')
    for key in stresses:
        print(f'\tElement {key}: {stresses[key]:.2f}')
    print('\nStrains:')
    for key in strains:
        print(f'\tElement {key}: {strains[key]:.4E}')


################
##  MAIN RUN  ##
################

model = set_geo_data()
model['elements'] = create_elements(model)
nodes = model['nodes']
elements = model['elements']
ext_forces = model['ext_forces']
plot_system(nodes, elements, ext_forces)

K = create_global_stiffness_matrix(nodes, elements, model['ext_forces'], model['areas'], model['ndofs'])
Fext = create_ext_force_vector(model['ext_forces'], model['ndofs'])
Kcondensed, Fcondensed = apply_boundary_conditions(K, Fext, model['restrained_dofs'])
U = calculate_displacements(Kcondensed, Fcondensed)
# print('U:\n',U)
F = calculate_reaction_forces(model['restrained_dofs'], model['ndofs'], U, K)
# print('F\n',F)
Ug = construct_global_displ_matrix(model['restrained_dofs'], model['ndofs'], U)
# print('UG:\n',Ug)
stresses, strains, internalForces = calculate_strain_stress(elements, Ug)

print_results(Ug, F, stresses, strains)

plt.show()
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import loguru
from BarElement2D import BarElement2D
    

def set_geo_data():
    nodes = {0: [0,0], 1: [5,0], 2: [0,5]}
    elements = {0: [0, 1], 1: [0, 2], 2: [1, 2]}
    ext_forces = {0: [10, 10], 1: [20, 30], 2: [0, 0]}
    # For each node we have u, v dof
    dofs = {0: [1, 1], 1: [1, 1], 2: [1, 1]}
    # The restrained dofs are the 1 and the free dofs are 0 (u,v)
    restrained_dofs = {0: [1, 1], 1: [0, 0], 2: [1, 1]}
    areas = {0: 1, 1: 1, 2: 1}
    stiffs = {0: 30e6, 1: 30e6, 2: 30e6}
    return {
        'nodes': nodes,
        'elements': elements,
        'ext_forces': ext_forces,
        'dofs': dofs,
        'restrained_dofs': restrained_dofs,
        'areas': areas,
        'stiffs': stiffs
    }
        
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
    ax.scatter(x, y, s=100, color='darkorange', zorder=4)
    for i, pos in enumerate(zip(x,y)):
        ax.text(pos[0]-0.04, pos[1]-0.05, str(i), color='black', fontsize=8, zorder=5)
    

    for e in elements.values():
        ax.plot((e.fromPoint[0], e.toPoint[0]), (e.fromPoint[1], e.toPoint[1]), '-', color='slategray', zorder=2)
        idpos = abs((e.fromPoint - e.toPoint) / 2) + 0.04
        ax.text(idpos[0], idpos[1], str(e.elid), color='black', fontsize=8, zorder=5)
    
    maxforce = max(forces[max(forces, key=lambda key: forces[key])][0], forces[max(forces, key=lambda key: forces[key])][1])
    for fid in forces.keys():
        node = nodes[fid]
        fx = forces[fid][0]
        fy = forces[fid][1]
        if fx != 0 or fy != 0:
            ax.arrow(node[0], node[1], fx/round(maxforce*10/ax.get_xlim()[1]), fy/round(maxforce*10/ax.get_ylim()[1]), head_width=0.08, zorder=3)
    
    return True


if __name__ == '__main__':
    model = set_geo_data()
    model['elements'] = create_elements(model)
    nodes = model['nodes']
    elements = model['elements']
    ext_forces = model['ext_forces']
    plot_system(nodes, elements, ext_forces)
    plt.show()
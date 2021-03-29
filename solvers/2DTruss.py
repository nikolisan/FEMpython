import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import loguru
from BarElement2D import BarElement2D
    

def set_geo_data():
    nodes = {0: [0,0], 1: [5,0], 2: [0,5]}
    elements = {0: [0, 1], 1: [0, 2], 2: [1, 2]}
    ext_forces = {0: [0, 0], 1: [200, -200], 2: [0, 0]}
    dofs = {0: [0, 1], 1: [0, 1], 2: [0, 1]}
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
    x = [val[0] for val in nodes.values()]
    y = [val[1] for val in nodes.values()]
    plt.scatter(x, y, s=100, color='darkorange', zorder=3)
    for i, pos in enumerate(zip(x,y)):
        plt.text(pos[0]-0.04, pos[1]-0.05, str(i), color='black', fontsize=8, zorder=4)

    for e in elements.values():
        plt.plot((e.fromPoint[0], e.toPoint[0]), (e.fromPoint[1], e.toPoint[1]), '-', color='slategray', zorder=1)
        idpos = abs((e.fromPoint - e.toPoint) / 2) + 0.04
        plt.text(idpos[0], idpos[1], str(e.elid), color='black', fontsize=8, zorder=4)
    return True


if __name__ == '__main__':
    model = set_geo_data()
    model['elements'] = create_elements(model)
    nodes = model['nodes']
    elements = model['elements']
    ext_forces = model['ext_forces']
    plot_system(nodes, elements, ext_forces)
    plt.show()
import numpy as np
from scipy.linalg import inv
from matplotlib import colors
import matplotlib.pyplot as plt
import traceback
from loguru import logger
import random

from ..elements.TriangularElement2D import TriangularPlaneStressElement as StressElement
from ..elements.TriangularElement2D import TriangularPlaneStrainElement as StrainElement

class PlaneTriangular2D():
    def __init__(self, model, element_type='stress'):
        self.model = model
        self._element_type = element_type

    # Helper functions
    def flatten_dofs(self, dofs:dict, restrained: bool):
        try:
            out = []
            for key in dofs.keys():
                if dofs[key][0] == restrained:
                    out.append(key*2)
                if dofs[key][1] == restrained:
                    out.append(key*2 + 1)
            return out
        except Exception as e:
            logger.exception(e)
            return None

    # Main solver functions
    def create_elements(self):
        logger.info('Creating Elements')
        try:
            elements = self.model['elements']
            nodes = self.model['nodes']
            thick = self.model['thick']
            stiff = self.model['stiff']
            poisson = self.model['poisson']

            new_elements = {}
            for id, element in enumerate(elements.values()):
                points = np.array([nodes[element[0]], nodes[element[1]], nodes[element[2]]], dtype='float')
                if self._element_type.lower() == 'stress':
                    el = StressElement(id, points, element, thick, stiff, poisson)
                elif self._element_type.lower() == 'strain':
                    el = StrainElement(id, points, element, thick, stiff, poisson)
                else:
                    logger.error('Valid element types are: "stress" or "strain"')
                new_elements[id] = el
            
            self.model['elements'] = new_elements
            return True
        except Exception as e:
            logger.exception(e)
            return False

    def create_global_stiffness_matrix(self):
        try:
            elements = self.model['elements']
            ndofs = self.model['ndofs']
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
            logger.exception(e)
            return None

    def create_force_vector(self):
        try:
            forces = self.model['forces']
            ndofs = self.model['ndofs']
            F = []
            for force in forces.values():
                F.extend(force)
            return np.array(F)
        except Exception as e:
            logger.exception(e)
            return None

    def apply_boundary_conditions(self):
        try:
            K = self.K
            F = self.F
            restrained_dofs = self.model['rdofs']
            logger.info('Applying boundary conditions')
            rdofs = self.flatten_dofs(restrained_dofs, True)
            rdofs = np.array(rdofs)
            K = np.delete(K, rdofs, axis=0)
            K = np.delete(K, rdofs, axis=1)
            F = np.delete(F,rdofs)
            return K, F
        except Exception as e:
            logger.exception(e)
            return None

    def calculate_displacements(self, K, F):
        try:
            logger.info('Calculating displacements')
            return inv(K).dot(F)
        except Exception as e:
            logger.exception(e)
            print(None)

    def create_global_displacement_matrix(self, U):
        try:
            restrained_dofs = self.model['rdofs']
            ndofs = self.model['ndofs']
            fdofs = self.flatten_dofs(restrained_dofs, False)
            Ug = np.zeros(ndofs)
            for i in range(len(fdofs)):
                Ug[fdofs[i]] = U[i]
            return Ug
        except Exception as e:
            logger.exception(e)
            return None
    
    def store_displ_on_element(self):
        try:
            Ug = self.Ug
            elements = self.model['elements']
            for element in elements.values():
                Ue = Ug[element.dofs]
                element.displacements = Ue
        except Exception as e:
            logger.exception(e)

    def calculate_stress(self):
        try:
            elements = self.model['elements']
            logger.info('Calculating stress')
            for element in elements.values():
                stress = element.elasticity_matrix.dot(element.b_matrix).dot(element.displacements)
                element.stress = stress
        except Exception as e:
            logger.exception(e)

    def calculate_strain(self):
        try:
            elements = self.model['elements']
            logger.info('Calculating strain')
            for element in elements.values():
                strain = inv(element.elasticity_matrix).dot(element.stress)
                element.strain = strain
        except Exception as e:
            logger.exception(e)

    # Solve the system
    def solve(self):
        # Create elements
        if isinstance(self.model['elements'][0], list):
            self.create_elements()
        # Create global K
        self.K = self.create_global_stiffness_matrix()
        if self.K.all() == None:
            logger.error('Error in creating K matrix')
            return False
        # Create F vector
        self.F = self.create_force_vector()
        if self.F.all() == None:
            logger.error('Error in creating F matrix')
            return False
        # Apply boundaries
        Kc, Fc = self.apply_boundary_conditions()
        # Calculate Displs
        U = self.calculate_displacements(Kc, Fc)
        if U.all() == None:
            logger.error('Error in calculating Displacements')
            return False
        # Create Ug
        self.Ug = self.create_global_displacement_matrix(U)
        if self.Ug.all() == None:
            logger.error('Error in calculating Global Displacements')
            return False
        try:
            # Store Displ
            self.store_displ_on_element()
            # Calculate stress
            self.calculate_stress()
            # Calculate strain
            self.calculate_strain()
            return True
        except Exception as e:
            logger.exception(e)
            return False

    # Results
    def print_results(self):
        Ug = self.Ug
        nodes = self.model['nodes']
        elements = self.model['elements']

        print('\nResults:')
        print('Nodal Displacements:')
        for node in range(len(nodes)):
            print(f' - Node: {node}')
            print(f'\t u_{node}: {Ug[node*2]:+.5f}'.replace('+', ' '))
            print(f'\t v_{node}: {Ug[node*2+1]:+.5f}'.replace('+', ' '))
        print('\nElement Stress:')
        for element in elements.values():
            print(f' - Element: {element.id}')
            print(f'\t σ_x : {element.stress[0]:+.3e}')
            print(f'\t σ_y : {element.stress[1]:+.3e}')
            print(f'\t τ_xy: {element.stress[2]:+.3e}')
        print('\nElement Strain:')
        for element in elements.values():
            print(f' - Element: {element.id}')
            print(f'\t ε_x : {element.strain[0]:+.3e}')
            print(f'\t ε_y : {element.strain[1]:+.3e}')
            print(f'\t γ_xy: {element.strain[2]:+.3e}')

    # Plot the model geometry
    def plot_model(self, show_nodes=True, show_text=True, blocking=True):
        name = self.model['name']
        nodes = self.model['nodes']
        forces = self.model['forces']

        if isinstance(self.model['elements'][0], list):
            self.create_elements()
        
        elements = self.model['elements']

        try:
            fig, ax = plt.subplots(num=f'PlaneTriangular2D: Plane {self._element_type} Elements - Model: {name}')
            ax.set_aspect('equal')
            ax.set_title(f'Model: {name}')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            x = [val[0] for val in nodes.values()]
            y = [val[1] for val in nodes.values()]
            tri = [el.nodes for el in elements.values()]
            
            # Plot nodes and node ids
            if show_nodes:
                ax.scatter(x, y, s=100, color='darkorange', zorder=4)           
            ax.triplot(x, y, triangles=tri, color='slategray')

            if show_text:
                for i, pos in enumerate(zip(x,y)):
                    ax.text(pos[0], pos[1], str(i), color='black', ha='center', va='center', fontsize=8, zorder=5)
                for el in elements.values():
                    ax.text(el.centroid[0], el.centroid[1], str(el.id), color='black', ha='center', va='center', fontsize=10, zorder=5)
            
            for i, force in enumerate(forces.values()):
                if force[0] or force[1]:
                    origin = np.asarray(nodes[i])
                    Fx = np.asarray(force[0])
                    Fy = np.asarray(force[1])
                    print(origin, Fx, Fy)
                    ax.quiver(*origin, Fx, Fy, pivot='tip', zorder=5)

            plt.show(block=blocking)
        except Exception as e:
            logger.exception(e)
            return None

    def plot_stresses(self, show_nodes=True, show_text=True):
        try:
            name = self.model['name']
            nodes = self.model['nodes']
            elements = self.model['elements']
            # Get values from the model and the elements
            x = np.asarray([val[0] for val in nodes.values()])
            y = np.asarray([val[1] for val in nodes.values()])
            tri = np.asarray([el.nodes for el in elements.values()])
            sx = np.asarray([el.stress[0] for el in elements.values()])
            sy = np.asarray([el.stress[1] for el in elements.values()])
            txy = np.asarray([el.stress[2] for el in elements.values()])
            stress = [sx, sy, txy]
            
            for i in range(3):
                stress_name = ['σ_x', 'σ_y', 'τ_xy']
                fig , ax = plt.subplots(num=f'PlaneTriangular2D: Stress {stress_name[i]} - Model: {name}')
                ax.set_aspect('equal')
                ax.set_xlabel('X')
                ax.set_ylabel('Y')
                # Plot nodes as black dots
                if show_nodes:
                    ax.scatter(x, y, color='k', zorder=4)
                ax.set_title(f'Model: {name} - Stress: {stress_name[i]}')
                tpc = ax.tripcolor(x, y, tri, facecolors=stress[i], edgecolor='k', norm=colors.CenteredNorm(), cmap='bwr')
                fig.colorbar(tpc)
                if show_text:
                    for el in elements.values():
                        ax.text(el.centroid[0], el.centroid[1], str(f'{el.stress[i]:+.3e}'), color='black', ha='center', va='center', fontsize=10, zorder=5)
            plt.show()
        except Exception as e:
            logger.exception(e)
            return None


import calfem.utils as cfu
import calfem.core as cfc
import numpy as np
import scipy.sparse as sp
from plantml import plantml
from constants import top_element, half_circle, circle1, circle2, circle3, T_out, T_in, T_inf, boundary_markers, thickness, alphan, alphac, cp, rho, k

class pre_processor():
    def __init__(self, mesh):
        self.coords, self.edof, self.dofs, self.bdofs, self.elementmarkers, self.boundary_elements = mesh.create()
        self.num_dofs = len(self.dofs)

    def stationary_boundary_conditions(self):
        bc = np.array([], 'i')
        bcVal = np.array([], 'f')

        # Apply the boundary conditions for the stationary problem
        bc, bcVal = cfu.apply_bc(self.bdofs, bc, bcVal, half_circle, T_out, 0)
        bc, bcVal = cfu.apply_bc(self.bdofs, bc, bcVal, circle1, T_in, 0)
        bc, bcVal = cfu.apply_bc(self.bdofs, bc, bcVal, circle2, T_out, 0)
        bc, bcVal = cfu.apply_bc(self.bdofs, bc, bcVal, circle3, T_in, 0)
        bc, bcVal = cfu.apply_bc(self.bdofs, bc, bcVal, top_element, T_inf, 2)
        return bc, bcVal

    def convection_boundary_conditions(self):
        boundary_node_list = []
        K_c = np.zeros([self.num_dofs, self.num_dofs])
        F_c = np.zeros([self.num_dofs, 1])

        # Retrieve the boundary nodes
        for marker in boundary_markers:
            for i in range(len(self.boundary_elements[marker])):
                nodes = self.boundary_elements[marker][i]['node-number-list']
                boundary_node_list.append([nodes, marker])

        # Compute and apply the boundary conditions for the convection problem
        for nodes in boundary_node_list:
            nx0, ny0 = self.coords[nodes[0][0] - 1]
            nx1, ny1 = self.coords[nodes[0][1] - 1]

            L = np.linalg.norm(np.array([nx1-nx0, ny1-ny0]))

            if nodes[1] == top_element:
              alpha = alphan
              temp = T_inf

            else:
                alpha = alphac 

                if nodes[1] == half_circle or nodes[1] == circle2:
                    temp = T_out

                else:
                    temp = T_in

            K_c[nodes[0][0], nodes[0][0]] += alpha*thickness*L/3 
            K_c[nodes[0][1], nodes[0][1]] += alpha*thickness*L/3 
            K_c[nodes[0][0], nodes[0][1]] += alpha*thickness*L/6
            K_c[nodes[0][1], nodes[0][0]] += alpha*thickness*L/6

            F_c[nodes[0][0]] += temp*alpha*thickness*L/2
            F_c[nodes[0][1]] += temp*alpha*thickness*L/2
        
        K_c = sp.csr_matrix(K_c)
        F_c = sp.csr_matrix(F_c)

        return K_c, F_c

    def transient_pre_processor(self):
        K = np.zeros([self.num_dofs, self.num_dofs])
        F = np.zeros([self.num_dofs, 1])
        V = np.zeros([self.num_dofs , 1])
        C = np.zeros([self.num_dofs, self.num_dofs])
        ex, ey = cfc.coordxtr(self.edof, self.coords, self.dofs)
        
        ptype = 1
        ep = [ptype, thickness]
        D = k * np.eye(2)

        # Compute the K and C matrix of the problem
        for eltopo, elx, ely in zip(self.edof, ex, ey):
            Ke = cfc.flw2te(elx, ely, ep, D)
            Ce = plantml(elx, ely, cp*rho)
            cfc.assem(eltopo, C, Ce)
            cfc.assem(eltopo, K, Ke)

        # Assign the same boundary conditions as in the convection problem
        K_c, F_c = self.convection_boundary_conditions()    

        K += K_c
        F += F_c

        # Compute the surface integral of each element
        for e in self.edof:
            L1 = np.abs(np.subtract(self.coords[e[0] - 1], self.coords[e[1] - 1]))
            L2 = np.abs(np.subtract(self.coords[e[0] - 1], self.coords[e[2] - 1]))

            A = np.abs(L1[0] * L2[1] - L1[1] * L2[0])
            Ve = A*thickness/3

            V[e[0] - 1] += Ve
            V[e[1] - 1] += Ve
            V[e[2] - 1] += Ve   

        return K, F, V, C
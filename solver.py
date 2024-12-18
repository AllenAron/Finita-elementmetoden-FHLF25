import calfem.core as cfc
import calfem.utils as cfu
from pre_processor import pre_processor
import numpy as np
import scipy.sparse as sp
import colorama

from constants import bottomright_element, left_element, T_initial, thickness, k, E, nu, alpha
from strategies import f1, f2


class solver():

    def __init__(self, mesh):
        self.coords, self.edof, self.dofs, self.bdofs, self.elementmarkers, self.boundary_elements = mesh.create()
        self.num_dofs = np.size(self.dofs)
        self.prep = pre_processor(mesh)

    def stationary_solver(self):
        K = np.zeros([self.num_dofs, self.num_dofs])
        F = np.zeros([self.num_dofs, 1])
        ex, ey = cfc.coordxtr(self.edof, self.coords, self.dofs)
        ptype = 1
        ep = [ptype, thickness]
        D = k * np.eye(2)

        # Calculate the stiffness matrix for the problem
        for eltopo, elx, ely in zip(self.edof, ex, ey):
            Ke = cfc.flw2te(elx, ely, ep, D)
            cfc.assem(eltopo, K, Ke)

        bc, bcVal = self.prep.stationary_boundary_conditions()

        # Solve the equation system with stationary boundary conditions
        a, r = cfc.solveq(K, F, bc, bcVal)
        return a

    def convection_solver(self):
        K = np.zeros([self.num_dofs, self.num_dofs])
        F = np.zeros([self.num_dofs, 1])
        ex, ey = cfc.coordxtr(self.edof, self.coords, self.dofs)
        ptype = 1
        ep = [ptype, thickness]
        D = k * np.eye(2)
        
        # Assemble the stiffness matrix for the problem
        print('Assembling stiffness matrix...')
        idx = 0
        for eltopo, elx, ely in zip(self.edof, ex, ey):
            Ke = cfc.flw2te(elx, ely, ep, D)
            cfc.assem(eltopo, K, Ke)
            idx += 1
            self.progress_bar(idx, len(self.edof))

        # Assign the convection boundary conditions
        K_c, F_c = self.prep.convection_boundary_conditions()

        K += K_c
        F += F_c
        K = sp.csr_matrix(K)
        F = sp.csr_matrix(F)
        a = sp.linalg.spsolve(K, F) # Solve the equation system

        return a

    def transient_solver(self, t0, t1, N):
        print('Performing pre-processing...')
        K_b, F_b, V, C = self.prep.transient_pre_processor()
        A = np.zeros((N, self.num_dofs))
        h = (t1-t0) / N
        tgrid = np.linspace(t0, t1, N)[1:]

        K = K_b + C/h
        K = sp.csr_matrix(K)
        F = sp.csr_matrix(np.zeros([self.num_dofs, 1]))
        
        # Assign the initial temperatures
        for i in range(self.num_dofs):
            A[0][i] = T_initial

        # Solve the transient equation system
        print('Solving transient equation system...')
        for idx, t in enumerate(tgrid):
            F = F_b + V * f1(t) + np.array([np.matmul(C/h, A[idx])]).T
            a = sp.linalg.spsolve(K, F)
            A[idx + 1] = a
            self.progress_bar(t, t1)


        # Find the time-frame with maximum and minimum temperature
        tmax = 0
        tmin = 0
        B = np.abs(A-T_initial)
        prevmax = np.max(B[0])
        prevmin = np.min(A[0])
        for i in range(N):
            nextmax = np.max(B[i])
            nextmin = np.min(A[i])
            if nextmax > prevmax:
                prevmax = nextmax
                tmax = i
            if nextmin < prevmin:
                prevmin = nextmin
                tmin = i

        return A, tgrid, tmax, tmin

    def von_mises_solver(self, temperature):
        ptype = 2
        ep = [ptype, thickness]
        D = cfc.hooke(ptype, E, nu)
        ex, ey = cfc.coordxtr(self.edof, self.coords, self.dofs)
        K = np.zeros([self.num_dofs, self.num_dofs])
        F = np.zeros([self.num_dofs, 1])

        # Assemble the stiffness matrix for the stress problem
        print('Assembling stiffness matrix...')
        idx = 0
        for eltopo, elx, ely in zip(self.edof, ex, ey):
            Ke = cfc.plante(elx, ely, ep, D)
            cfc.assem(eltopo, K, Ke)
            idx += 1
            self.progress_bar(idx, len(self.edof))


        # Calculate the thermal stresses
        delta_T_list = []
        idx = 0
        print('Calculating thermal stresses...')
        for element, elx, ely in zip(self.edof, ex, ey):
            delta_T = 0
            for dof in element:
                delta_T += (temperature[int(np.round(dof/2)) - 1] - T_initial) / 6

            delta_T_list.append(delta_T)
            eps0 = (1)*alpha*delta_T*np.array([1, 1, 0])
            xi, xj, xk = elx[0], elx[1], elx[2]
            yi, yj, yk = ely[0], ely[1], ely[2]
            
            B = 1/2*np.matrix(
                [[yj-yk, 0, yk-yi, 0, yi-yj, 0],
                 [0, xk-xj, 0, xi-xk, 0, xj-xi],
                 [xk-xj, yj-yk, xi-xk, yk-yi, xj-xi, yi-yj]]
            )
            D0 = np.matmul(np.delete(np.delete(D, 2, axis=0), 2, axis=1), eps0)
            f0 = np.matmul(D0, B) * thickness

            for dof, f in zip(element, f0.T):
                F[dof - 1] += f[0,0]
            idx += 1
            self.progress_bar(idx, len(self.edof))

        # Assign the boundary conditions for the stress
        bc = np.array([],'i')
        bcVal = np.array([],'f')
        bc, bcVal = cfu.applybc(self.bdofs, bc, bcVal, bottomright_element, 0, 0)
        bc, bcVal = cfu.applybc(self.bdofs, bc, bcVal, left_element, 0, 1)

        print('Solving for boundary conditions...')
        a,r = cfc.solveq(K,F,bc,bcVal)


        ed = cfc.extract_eldisp(self.edof, a)
        vonMises = np.zeros(self.edof.shape[0])

        print('Calculating von Mises stresses...')
        for i in range(self.edof.shape[0]):
            # Determine element stresses and strains in the element.
            es, et = cfc.plants(ex[i,:], ey[i,:], ep, D, ed[i,:])
            
            # Calculate and append effective stress to list.
            es = es[0]
            sigeff = np.sqrt(np.power(es[0],2) + np.power(es[1],2) + np.power(es[2],2) - es[0]*es[1] - es[0]*es[2] - es[1]*es[2] + 3*np.power(es[3],2))
            vonMises[i] = sigeff
            ## es: [sigx sigy tauxy]
            ## et: [epsx epsy gamxy]
            self.progress_bar(i, len(self.edof) - 1)
            
        return vonMises, a
    
    # Method of displaying the progress in the terminal
    def progress_bar(self, progress, total):
        percent = 100 * (progress / float(total))
        bar = '█' * int(percent) + '-' * (100 - int(percent))
        color = colorama.Fore.RED
        if percent > 33:
            color = colorama.Fore.YELLOW
        if percent > 66:
            color = colorama.Fore.GREEN
        print(color + f"\r|{bar}| {percent:.2f}%", end='\r')
        if progress == total:
            print(colorama.Fore.MAGENTA + f"\r|{bar}| {percent:.2f}%", end='\r')
            print(colorama.Fore.RESET)
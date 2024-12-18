import calfem.vis_mpl as cfv
import numpy as np
import matplotlib.pyplot as plt

from constants import T_initial

class visualizer():
    def __init__(self, geometry):
        self.geometry = geometry

    def visualize_temperature_solution(self, data, mesh, vmin, vmax, draw_mesh = True, draw_geometry = True, cmap = 'turbo'):
        coords, edof, dofs, bdofs, elementmarkers, boundary_elements = mesh.create()
        left_coords = np.array([[-c[0], c[1]] for c in coords])
        
        cfv.figure()

        # Apply the colormap
        cfv.draw_nodal_values_shaded(data, coords, edof, colormap = cmap, title='', vmin=vmin, vmax=vmax)
        cfv.draw_nodal_values_shaded(data, left_coords, edof, colormap = cmap, title='', vmin=vmin, vmax=vmax)
        cfv.colorbar('horizontal')

        # Draw the mesh
        if draw_mesh:
            cfv.drawMesh(
            coords=coords,
            edof=edof,
            dofs_per_node=mesh.dofsPerNode,
            el_type=mesh.elType,
            filled=False,
            title=""
            )

        # Showing geometry
        if draw_geometry:
            cfv.draw_geometry(self.geometry)

        # Set axis and colorbar label
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        for c in cfv.gca().collections:
            if c.colorbar:
                c.colorbar.set_label('Temperature [K]')

    def visualize_stress_solution(self, data, disp, mesh, dofsPerNode, elType, draw_values = True, draw_displacements = False):
        coords, edof, dofs, bdofs, elementmarkers, boundary_elements = mesh.create()
        left_coords = np.array([[-c[0], c[1]] for c in coords])
        cfv.figure()
        if draw_values:
            cfv.draw_element_values(data, coords, edof, dofsPerNode, elType, None, draw_elements=False, draw_undisplaced_mesh=False, title="von Mises stress", cmap='turbo')
            cfv.draw_element_values(data, left_coords, edof, dofsPerNode, elType, None, draw_elements=False, draw_undisplaced_mesh=False, title="von Mises stress", cmap='turbo')
            cfv.colorbar('horizontal')
            
            plt.xlabel('x [m]')
            plt.ylabel('y [m]')
            for c in cfv.gca().collections:
                if c.colorbar:
                    c.colorbar.set_label('von Mises stress [Pa]')
        if draw_displacements:
            cfv.figure()
            #cfv.draw_mesh(coords, edof, dofsPerNode, elType, title='Mesh')
            cfv.draw_displacements(disp, coords, edof, dofsPerNode, elType, draw_undisplaced_mesh=True, title='Displacements', magnfac=25)
            cfv.draw_displacements(disp, left_coords, edof, dofsPerNode, elType, draw_undisplaced_mesh=True, title='Displacements', magnfac=25)
            plt.xlabel('x [m]')
            plt.ylabel('y [m]')

            cfv.figure()
            cfv.draw_element_values(data, coords, edof, dofsPerNode, elType, disp, draw_elements=False, draw_undisplaced_mesh=True, title="Displacements", cmap='turbo', magnfac=25)
            cfv.draw_element_values(data, left_coords, edof, dofsPerNode, elType, disp, draw_elements=False, draw_undisplaced_mesh=True, title="Displacements", cmap='turbo', magnfac=25)
            cfv.colorbar('horizontal')
            plt.xlabel('x [m]')
            plt.ylabel('y [m]')
        
    def draw_mesh(self, mesh):
        coords, edof, dofs, bdofs, elementmarkers, boundary_elements = mesh.create()
        left_coords = np.array([[-c[0], c[1]] for c in coords])
        cfv.figure()
        cfv.drawMesh(coords=coords, edof=edof, dofs_per_node=mesh.dofsPerNode, el_type=mesh.elType, filled=False, title="Mesh")
        cfv.drawMesh(coords=left_coords, edof=edof, dofs_per_node=mesh.dofsPerNode, el_type=mesh.elType, filled=False, title="Mesh")

    def show(self):
        cfv.show_and_wait()

    def plotminmax(self, data, tgrid, N):
        a = np.zeros(N-1)
        for i in range(N- 1):
            a[i] = np.max(data[i])
        plt.figure()
        plt.plot(tgrid,a)
        plt.title('Maximum temperature')
        plt.xlabel('Time [s]')
        plt.ylabel('Temperature [K]')
        plt.show()

        for i in range(N-1):
            a[i] = np.min(data[i])
        plt.figure()
        plt.plot(tgrid,a)
        plt.title('Minimum temperature')
        plt.xlabel('Time [s]')
        plt.ylabel('Temperature [K]')
        plt.show()

        for i in range(N - 1):
            a[i] = np.max(np.abs(data[i]-T_initial))
        plt.figure()
        plt.plot(tgrid,a)
        plt.title('Maximum temperature deviation')
        plt.xlabel('Time [s]')
        plt.ylabel('Temperature deviation [K]')
        plt.show()